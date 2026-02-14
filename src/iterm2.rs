/// iTerm2 Graphics Protocol encoder for headless terminal rendering.
///
/// Encodes RGBA pixel buffers into iTerm2 inline image escape sequences
/// (PNG + base64) for display in compatible terminals (WezTerm, iTerm2).
pub struct Iterm2Encoder {
    png_buf: Vec<u8>,
    b64_buf: String,
    seq_buf: Vec<u8>,
}

impl Iterm2Encoder {
    pub fn new() -> Self {
        Self {
            png_buf: Vec::new(),
            b64_buf: String::new(),
            seq_buf: Vec::new(),
        }
    }

    /// Encode an RGBA pixel buffer into an iTerm2 inline image escape sequence.
    /// Returns the complete escape sequence bytes ready for stdout.
    pub fn encode(&mut self, rgba: &[u8], width: usize, height: usize) -> &[u8] {
        use base64::Engine;
        use std::io::Write;

        // 1. Encode RGBA directly to PNG (RGB strip happens inside the encoder)
        self.png_buf.clear();
        {
            let mut encoder = png::Encoder::new(&mut self.png_buf, width as u32, height as u32);
            encoder.set_color(png::ColorType::Rgb);
            encoder.set_depth(png::BitDepth::Eight);
            encoder.set_compression(png::Compression::Fast);
            let mut writer = encoder.write_header().expect("PNG header write");
            // Strip alpha: RGBA → RGB
            let rgb: Vec<u8> = rgba
                .chunks_exact(4)
                .flat_map(|px| [px[0], px[1], px[2]])
                .collect();
            writer.write_image_data(&rgb).expect("PNG data write");
        }

        // 2. Base64 encode the PNG bytes
        self.b64_buf.clear();
        base64::engine::general_purpose::STANDARD.encode_string(&self.png_buf, &mut self.b64_buf);

        // 3. Build escape sequence
        self.seq_buf.clear();
        write!(
            self.seq_buf,
            "\x1b]1337;File=inline=1;size={};width={}px;height={}px;preserveAspectRatio=0:{}\x07",
            self.png_buf.len(),
            width,
            height,
            self.b64_buf,
        )
        .expect("escape sequence write");

        &self.seq_buf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_rgba(width: usize, height: usize) -> Vec<u8> {
        let mut rgba = vec![0u8; width * height * 4];
        for y in 0..height {
            for x in 0..width {
                let i = (y * width + x) * 4;
                rgba[i] = (x * 255 / width.max(1)) as u8; // R
                rgba[i + 1] = (y * 255 / height.max(1)) as u8; // G
                rgba[i + 2] = 128; // B
                rgba[i + 3] = 255; // A
            }
        }
        rgba
    }

    #[test]
    fn test_encode_produces_valid_escape_sequence() {
        let mut encoder = Iterm2Encoder::new();
        let rgba = make_test_rgba(32, 16);
        let result = encoder.encode(&rgba, 32, 16);

        assert!(!result.is_empty(), "Encoded output should not be empty");

        // Must start with ESC ]1337;File=
        let prefix = b"\x1b]1337;File=";
        assert!(
            result.starts_with(prefix),
            "Should start with iTerm2 escape prefix, got {:?}",
            &result[..result.len().min(20)]
        );

        // Must end with BEL (\x07)
        assert_eq!(
            result.last(),
            Some(&0x07),
            "Should end with BEL character"
        );

        // Must contain inline=1
        let as_str = std::str::from_utf8(result).expect("Should be valid UTF-8");
        assert!(as_str.contains("inline=1"), "Should contain inline=1");
    }

    #[test]
    fn test_encode_contains_base64_png() {
        let mut encoder = Iterm2Encoder::new();
        let rgba = make_test_rgba(16, 8);
        let result = encoder.encode(&rgba, 16, 8);
        let as_str = std::str::from_utf8(result).expect("Should be valid UTF-8");

        // Extract base64 payload after the colon separator
        let colon_pos = as_str.rfind(':').expect("Should contain colon separator");
        let b64_data = &as_str[colon_pos + 1..as_str.len() - 1]; // strip trailing BEL

        // Decode base64
        use base64::Engine;
        let png_bytes = base64::engine::general_purpose::STANDARD
            .decode(b64_data)
            .expect("Base64 payload should decode");

        // Verify PNG header
        assert!(
            png_bytes.starts_with(&[0x89, b'P', b'N', b'G']),
            "Decoded payload should be a valid PNG (starts with PNG header)"
        );
    }

    #[test]
    fn test_encode_reuse_buffers() {
        let mut encoder = Iterm2Encoder::new();

        // First encode
        let rgba1 = make_test_rgba(8, 6);
        let result1 = encoder.encode(&rgba1, 8, 6);
        assert!(!result1.is_empty());
        let len1 = result1.len();

        // Second encode with different dimensions
        let rgba2 = make_test_rgba(16, 12);
        let result2 = encoder.encode(&rgba2, 16, 12);
        assert!(!result2.is_empty());

        // Second should be larger (more pixels → more data)
        assert!(
            result2.len() > len1,
            "Larger image should produce larger output"
        );

        // Both should be valid escape sequences
        assert!(result2.starts_with(b"\x1b]1337;File="));
        assert_eq!(result2.last(), Some(&0x07));
    }

    #[test]
    fn test_encode_small_image() {
        let mut encoder = Iterm2Encoder::new();
        let rgba = make_test_rgba(8, 6);
        let result = encoder.encode(&rgba, 8, 6);

        assert!(!result.is_empty(), "Should encode even small images");
        assert!(result.starts_with(b"\x1b]1337;File="));
        assert_eq!(result.last(), Some(&0x07));

        // Verify dimensions are in the header
        let as_str = std::str::from_utf8(result).expect("Valid UTF-8");
        assert!(as_str.contains("width=8px"), "Should specify width");
        assert!(as_str.contains("height=6px"), "Should specify height");
    }
}
