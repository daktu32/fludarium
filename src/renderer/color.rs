/// Tokyo Night-inspired color stops for field mapping.
/// Deep navy -> blue -> purple -> pink -> orange
pub(crate) const COLOR_STOPS: [(f64, f64, f64); 5] = [
    (26.0, 27.0, 38.0),    // #1a1b26 navy         (0.00)
    (122.0, 162.0, 247.0), // #7aa2f7 blue         (0.25)
    (187.0, 154.0, 247.0), // #bb9af7 purple       (0.50)
    (247.0, 118.0, 142.0), // #f7768e pink         (0.75)
    (255.0, 158.0, 100.0), // #ff9e64 orange       (1.00)
];

/// Convert temperature [0.0, 1.0] to RGBA color (Tokyo Night palette).
pub fn temperature_to_rgba(t: f64) -> [u8; 4] {
    let t = t.clamp(0.0, 1.0);
    let seg = t * 4.0;
    let i = (seg as usize).min(3);
    let s = seg - i as f64;

    let (r0, g0, b0) = COLOR_STOPS[i];
    let (r1, g1, b1) = COLOR_STOPS[i + 1];

    [
        (r0 + s * (r1 - r0)) as u8,
        (g0 + s * (g1 - g0)) as u8,
        (b0 + s * (b1 - b0)) as u8,
        255,
    ]
}

/// Color bar layout constants.
pub(crate) const BAR_GAP: usize = 6;
pub(crate) const BAR_WIDTH: usize = 20;
pub(crate) const TICK_LEN: usize = 4;
pub(crate) const LABEL_GAP: usize = 2;
pub(crate) const LABEL_WIDTH: usize = 24;
pub(crate) const BAR_TOTAL: usize = BAR_GAP + BAR_WIDTH + TICK_LEN + LABEL_GAP + LABEL_WIDTH;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_color_cold_is_navy() {
        let rgba = temperature_to_rgba(0.0);
        assert_eq!(rgba[0], 26, "R should be 26");
        assert_eq!(rgba[1], 27, "G should be 27");
        assert_eq!(rgba[2], 38, "B should be 38");
        assert_eq!(rgba[3], 255, "A should be 255");
    }

    #[test]
    fn test_color_hot_is_orange() {
        let rgba = temperature_to_rgba(1.0);
        assert_eq!(rgba[0], 255, "R should be 255");
        assert_eq!(rgba[1], 158, "G should be 158");
        assert_eq!(rgba[2], 100, "B should be 100");
    }

    #[test]
    fn test_color_mid_is_purple() {
        let rgba = temperature_to_rgba(0.5);
        assert_eq!(rgba[0], 187, "R should be 187");
        assert_eq!(rgba[1], 154, "G should be 154");
        assert_eq!(rgba[2], 247, "B should be 247");
    }

    #[test]
    fn test_color_clamp() {
        let lo = temperature_to_rgba(-1.0);
        let hi = temperature_to_rgba(2.0);
        assert_eq!(lo, temperature_to_rgba(0.0));
        assert_eq!(hi, temperature_to_rgba(1.0));
    }

    #[test]
    fn test_gradient_continuity() {
        let steps = 256;
        for i in 1..steps {
            let t0 = (i - 1) as f64 / (steps - 1) as f64;
            let t1 = i as f64 / (steps - 1) as f64;
            let c0 = temperature_to_rgba(t0);
            let c1 = temperature_to_rgba(t1);
            for ch in 0..3 {
                let diff = (c1[ch] as i32 - c0[ch] as i32).abs();
                assert!(
                    diff <= 5,
                    "Color channel {} jumped by {} between t={} and t={}",
                    ch, diff, t0, t1
                );
            }
        }
    }
}
