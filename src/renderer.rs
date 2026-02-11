use crate::state::{idx, SimState, N};

/// Convert temperature [0.0, 1.0] to RGBA color.
/// Blue(0.0) -> Cyan -> Green -> Yellow -> Red(1.0)
pub fn temperature_to_rgba(t: f64) -> [u8; 4] {
    let t = t.clamp(0.0, 1.0);

    let (r, g, b) = if t < 0.25 {
        // Blue -> Cyan
        let s = t / 0.25;
        (0.0, s, 1.0)
    } else if t < 0.5 {
        // Cyan -> Green
        let s = (t - 0.25) / 0.25;
        (0.0, 1.0, 1.0 - s)
    } else if t < 0.75 {
        // Green -> Yellow
        let s = (t - 0.5) / 0.25;
        (s, 1.0, 0.0)
    } else {
        // Yellow -> Red
        let s = (t - 0.75) / 0.25;
        (1.0, 1.0 - s, 0.0)
    };

    [
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
        255,
    ]
}

/// Display scale factor (nearest-neighbor upscale).
pub const SCALE: usize = 2;
pub const DISPLAY_SIZE: usize = N * SCALE;

/// Color bar layout constants.
const BAR_GAP: usize = 6;
const BAR_WIDTH: usize = 20;
const TICK_LEN: usize = 4;
/// Total output width including the color bar.
pub const FRAME_WIDTH: usize = DISPLAY_SIZE + BAR_GAP + BAR_WIDTH + TICK_LEN;
pub const FRAME_HEIGHT: usize = DISPLAY_SIZE;

/// Render temperature field + color bar to RGBA buffer.
/// Output size: FRAME_WIDTH x FRAME_HEIGHT x 4 bytes.
pub fn render(state: &SimState) -> Vec<u8> {
    let mut buf = vec![0u8; FRAME_WIDTH * FRAME_HEIGHT * 4];

    // Draw simulation (y-flipped)
    for screen_y in 0..DISPLAY_SIZE {
        let sim_y = (N - 1) - screen_y / SCALE;
        for screen_x in 0..DISPLAY_SIZE {
            let sim_x = screen_x / SCALE;
            let t = state.temperature[idx(sim_x as i32, sim_y as i32)];
            let rgba = temperature_to_rgba(t);
            let offset = (screen_y * FRAME_WIDTH + screen_x) * 4;
            buf[offset] = rgba[0];
            buf[offset + 1] = rgba[1];
            buf[offset + 2] = rgba[2];
            buf[offset + 3] = rgba[3];
        }
    }

    // Draw color bar: top=cold(0.0), bottom=hot(1.0)
    let bar_x = DISPLAY_SIZE + BAR_GAP;
    for y in 0..FRAME_HEIGHT {
        let t = y as f64 / (FRAME_HEIGHT - 1) as f64;
        let rgba = temperature_to_rgba(t);
        for bx in 0..BAR_WIDTH {
            let offset = (y * FRAME_WIDTH + bar_x + bx) * 4;
            buf[offset] = rgba[0];
            buf[offset + 1] = rgba[1];
            buf[offset + 2] = rgba[2];
            buf[offset + 3] = rgba[3];
        }
    }

    // Draw tick marks at 0%, 25%, 50%, 75%, 100%
    let tick_x = bar_x + BAR_WIDTH;
    for tick in 0..5u32 {
        let y = (tick as usize) * (FRAME_HEIGHT - 1) / 4;
        // White tick extending right from bar
        for dy in 0..2usize {
            let yy = (y + dy).min(FRAME_HEIGHT - 1);
            for tx in 0..TICK_LEN {
                let offset = (yy * FRAME_WIDTH + tick_x + tx) * 4;
                buf[offset] = 255;
                buf[offset + 1] = 255;
                buf[offset + 2] = 255;
                buf[offset + 3] = 255;
            }
        }
    }

    buf
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::SimState;

    #[test]
    fn test_color_cold_is_blue() {
        let rgba = temperature_to_rgba(0.0);
        assert_eq!(rgba[0], 0, "R should be 0");
        assert_eq!(rgba[1], 0, "G should be 0");
        assert_eq!(rgba[2], 255, "B should be 255");
        assert_eq!(rgba[3], 255, "A should be 255");
    }

    #[test]
    fn test_color_hot_is_red() {
        let rgba = temperature_to_rgba(1.0);
        assert_eq!(rgba[0], 255, "R should be 255");
        assert_eq!(rgba[1], 0, "G should be 0");
        assert_eq!(rgba[2], 0, "B should be 0");
    }

    #[test]
    fn test_color_mid_is_green() {
        let rgba = temperature_to_rgba(0.5);
        assert_eq!(rgba[0], 0, "R should be 0 at midpoint");
        assert_eq!(rgba[1], 255, "G should be 255 at midpoint");
        assert_eq!(rgba[2], 0, "B should be 0 at midpoint");
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
        // Check that adjacent colors don't jump too much
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
                    ch,
                    diff,
                    t0,
                    t1
                );
            }
        }
    }

    #[test]
    fn test_render_buffer_size() {
        let state = SimState::new();
        let buf = render(&state);
        assert_eq!(buf.len(), FRAME_WIDTH * FRAME_HEIGHT * 4);
    }

    #[test]
    fn test_render_y_flip() {
        let state = SimState::new();
        let buf = render(&state);

        // Screen top (y=0) should be cold (blue), screen bottom should be hot (red)
        let top = &buf[0..4];
        let bottom_row = DISPLAY_SIZE - 1;
        let bottom = &buf[(bottom_row * FRAME_WIDTH) * 4..(bottom_row * FRAME_WIDTH) * 4 + 4];

        // Top should be blueish (cold = 0.0 -> blue)
        assert!(top[2] > top[0], "Top should be more blue than red");
        // Bottom should be reddish (hot = 1.0 -> red)
        assert!(bottom[0] > bottom[2], "Bottom should be more red than blue");
    }

    #[test]
    fn test_color_bar_gradient() {
        let state = SimState::new();
        let buf = render(&state);

        // Color bar top pixel should be blue (cold), bottom should be red (hot)
        let bar_x = DISPLAY_SIZE + BAR_GAP + BAR_WIDTH / 2;
        let top_offset = (0 * FRAME_WIDTH + bar_x) * 4;
        let bot_offset = ((FRAME_HEIGHT - 1) * FRAME_WIDTH + bar_x) * 4;
        // Top of bar = T=0.0 → blue
        assert!(buf[top_offset + 2] > buf[top_offset], "Bar top should be blue");
        // Bottom of bar = T=1.0 → red
        assert!(buf[bot_offset] > buf[bot_offset + 2], "Bar bottom should be red");
    }
}
