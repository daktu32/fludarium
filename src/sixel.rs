use std::io::{self, Write};

/// Encode RGBA buffer to Sixel format using icy_sixel.
pub fn encode_sixel(rgba: &[u8], width: usize, height: usize) -> Result<Vec<u8>, String> {
    // icy_sixel expects a flat pixel buffer
    let sixel_output = icy_sixel::sixel_string(
        rgba,
        width as i32,
        height as i32,
        icy_sixel::PixelFormat::RGBA8888,
        icy_sixel::DiffusionMethod::Atkinson,
        icy_sixel::MethodForLargest::Auto,
        icy_sixel::MethodForRep::Auto,
        icy_sixel::Quality::AUTO,
    )
    .map_err(|e| format!("Sixel encoding error: {}", e))?;

    Ok(sixel_output.into_bytes())
}

/// Output a Sixel-encoded frame to stdout.
/// Uses synchronized output (DEC 2026) and single write to minimize flicker.
pub fn output_frame(sixel_data: &[u8]) -> io::Result<()> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    // Build entire frame in one buffer: BSU + cursor home + sixel + ESU
    let mut buf = Vec::with_capacity(10 + 3 + sixel_data.len() + 10);
    buf.extend_from_slice(b"\x1b[?2026h"); // begin synchronized update
    buf.extend_from_slice(b"\x1b[H"); // cursor home
    buf.extend_from_slice(sixel_data);
    buf.extend_from_slice(b"\x1b[?2026l"); // end synchronized update
    handle.write_all(&buf)?;
    handle.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_sixel_not_empty() {
        // Create a small 8x6 red image
        let width = 8;
        let height = 6;
        let mut rgba = vec![0u8; width * height * 4];
        for i in 0..(width * height) {
            rgba[i * 4] = 255; // R
            rgba[i * 4 + 1] = 0; // G
            rgba[i * 4 + 2] = 0; // B
            rgba[i * 4 + 3] = 255; // A
        }

        let result = encode_sixel(&rgba, width, height);
        assert!(result.is_ok(), "Encoding should succeed");
        let data = result.unwrap();
        assert!(!data.is_empty(), "Output should not be empty");
    }

    #[test]
    fn test_sixel_dcs_header() {
        let width = 8;
        let height = 6;
        let rgba = vec![128u8; width * height * 4];

        let data = encode_sixel(&rgba, width, height).unwrap();
        let s = String::from_utf8_lossy(&data);

        // Sixel data should start with DCS (ESC P or \x90)
        assert!(
            s.starts_with("\x1bP") || s.starts_with("\u{0090}") || data[0] == 0x90,
            "Should start with DCS header, got: {:?}",
            &data[..data.len().min(10)]
        );
    }

    #[test]
    fn test_sixel_st_footer() {
        let width = 8;
        let height = 6;
        let rgba = vec![128u8; width * height * 4];

        let data = encode_sixel(&rgba, width, height).unwrap();
        let s = String::from_utf8_lossy(&data);

        // Sixel data should end with ST (ESC \ or \x9c)
        assert!(
            s.ends_with("\x1b\\") || s.ends_with("\u{009c}") || *data.last().unwrap() == 0x9c,
            "Should end with ST footer, got: {:?}",
            &data[data.len().saturating_sub(10)..]
        );
    }
}
