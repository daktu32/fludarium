// 1D line plot renderer for fludarium playback.
//
// Renders a single field as a line chart with axes, grid lines,
// and a field-name badge. Used when jm == 1 (1D periodic data).

use super::color::{self, ColorMap};
use super::font::{self, FONT_HEIGHT, FONT_WIDTH, STATUS_BAR_HEIGHT};

/// Layout configuration for the line plot.
pub struct LinePlotConfig {
    /// Total frame width in pixels.
    pub frame_width: usize,
    /// Total frame height in pixels (including status bar).
    pub frame_height: usize,
    /// Plot area left margin (for y-axis labels).
    pub margin_left: usize,
    /// Plot area right margin.
    pub margin_right: usize,
    /// Plot area top margin (for field badge).
    pub margin_top: usize,
    /// Plot area bottom margin (for x-axis labels, above status bar).
    pub margin_bottom: usize,
}

impl LinePlotConfig {
    pub fn new(win_width: usize, win_height: usize) -> Self {
        Self {
            frame_width: win_width,
            frame_height: win_height,
            margin_left: 55,
            margin_right: 15,
            margin_top: 35,
            margin_bottom: 25,
        }
    }

    /// Width of the plot area in pixels.
    pub fn plot_width(&self) -> usize {
        self.frame_width
            .saturating_sub(self.margin_left + self.margin_right)
            .max(4)
    }

    /// Height of the plot area in pixels.
    pub fn plot_height(&self) -> usize {
        self.frame_height
            .saturating_sub(self.margin_top + self.margin_bottom + STATUS_BAR_HEIGHT)
            .max(4)
    }
}

/// Render a 1D line plot into an RGBA buffer.
pub fn render_lineplot(
    buf: &mut Vec<u8>,
    cfg: &LinePlotConfig,
    data: &[f64],
    domain_length: f64,
    y_range: (f64, f64),
    field_name: &str,
    field_index: usize,
    field_count: usize,
    colormap: ColorMap,
) {
    let fw = cfg.frame_width;
    let fh = cfg.frame_height;
    let total = fw * fh * 4;
    buf.resize(total, 0);

    // Dark background (#0A0C14)
    for i in 0..fw * fh {
        let off = i * 4;
        buf[off] = 0x0A;
        buf[off + 1] = 0x0C;
        buf[off + 2] = 0x14;
        buf[off + 3] = 0xFF;
    }

    let pw = cfg.plot_width();
    let ph = cfg.plot_height();
    let x0 = cfg.margin_left;
    let y0 = cfg.margin_top;
    let (vmin, vmax) = y_range;
    let vrange = (vmax - vmin).max(1e-15);

    // --- Grid lines (dashed, dim) ---
    let grid_color: [u8; 4] = [0x22, 0x28, 0x35, 0xFF];
    let dash_len = 4;
    let gap_len = 4;

    // Horizontal grid lines (5 lines: 0%, 25%, 50%, 75%, 100%)
    for i in 0..=4 {
        let py = y0 + ph * i / 4;
        for px in x0..x0 + pw {
            let phase = (px - x0) % (dash_len + gap_len);
            if phase < dash_len {
                set_pixel(buf, fw, px, py, grid_color);
            }
        }
    }

    // Vertical grid lines (5 lines)
    for i in 0..=4 {
        let px = x0 + pw * i / 4;
        for py in y0..y0 + ph {
            let phase = (py - y0) % (dash_len + gap_len);
            if phase < dash_len {
                set_pixel(buf, fw, px, py, grid_color);
            }
        }
    }

    // --- Axis frame (solid, brighter) ---
    let axis_color: [u8; 4] = [0x44, 0x4C, 0x5C, 0xFF];
    // Left edge
    for py in y0..=y0 + ph {
        set_pixel(buf, fw, x0, py, axis_color);
    }
    // Bottom edge
    for px in x0..=x0 + pw {
        set_pixel(buf, fw, px, y0 + ph, axis_color);
    }

    // --- Y-axis labels ---
    let label_color: [u8; 3] = [0x66, 0x6C, 0x7C];
    for i in 0..=4 {
        let frac = 1.0 - i as f64 / 4.0;
        let val = vmin + frac * vrange;
        let py = y0 + ph * i / 4;
        let label = format_axis_value(val);
        let label_w = label.len() * (FONT_WIDTH + 1);
        let lx = if cfg.margin_left > label_w + 4 {
            cfg.margin_left - label_w - 4
        } else {
            1
        };
        let ly = if py >= FONT_HEIGHT / 2 {
            py - FONT_HEIGHT / 2
        } else {
            0
        };
        font::draw_text(buf, fw, lx, ly, &label, label_color);
    }

    // --- X-axis labels ---
    let x_labels = ["0", "", "", "", ""];
    let x_label_end = format_axis_value(domain_length);
    let x_labels_text = [x_labels[0], "", "", "", &x_label_end];
    for (i, label) in x_labels_text.iter().enumerate() {
        if label.is_empty() {
            continue;
        }
        let px = x0 + pw * i / 4;
        let label_w = label.len() * (FONT_WIDTH + 1);
        let lx = px.saturating_sub(label_w / 2);
        let ly = y0 + ph + 6;
        font::draw_text(buf, fw, lx, ly, label, label_color);
    }

    // --- Line plot (Bresenham) ---
    let line_rgba = color::map_to_rgba(0.7, colormap);
    let line_color: [u8; 4] = [line_rgba[0], line_rgba[1], line_rgba[2], 0xFF];

    if data.len() >= 2 {
        let n = data.len();
        // Convert data points to screen coordinates
        let to_screen = |i: usize| -> (isize, isize) {
            let sx = x0 as f64 + (i as f64 / (n - 1) as f64) * pw as f64;
            let t = ((data[i] - vmin) / vrange).clamp(0.0, 1.0);
            let sy = y0 as f64 + (1.0 - t) * ph as f64;
            (sx.round() as isize, sy.round() as isize)
        };

        let (mut px0, mut py0) = to_screen(0);
        for i in 1..n {
            let (px1, py1) = to_screen(i);
            draw_line(buf, fw, px0, py0, px1, py1, x0, y0, pw, ph, line_color);
            px0 = px1;
            py0 = py1;
        }
    }

    // --- Field badge (top-left, matching spherical renderer style) ---
    font::render_field_badge(
        buf,
        fw,
        y0 + ph, // display_height for badge clipping
        field_name,
        field_index,
        field_count,
    );
}

/// Format an axis value compactly.
fn format_axis_value(val: f64) -> String {
    let abs = val.abs();
    if abs == 0.0 {
        "0".to_string()
    } else if abs >= 1000.0 || abs < 0.01 {
        format!("{:.1e}", val)
    } else if abs >= 10.0 {
        format!("{:.1}", val)
    } else if abs >= 1.0 {
        format!("{:.2}", val)
    } else {
        format!("{:.3}", val)
    }
}

/// Set a single pixel in the RGBA buffer.
#[inline]
fn set_pixel(buf: &mut [u8], fw: usize, x: usize, y: usize, rgba: [u8; 4]) {
    let off = (y * fw + x) * 4;
    if off + 3 < buf.len() {
        buf[off] = rgba[0];
        buf[off + 1] = rgba[1];
        buf[off + 2] = rgba[2];
        buf[off + 3] = rgba[3];
    }
}

/// Bresenham line drawing clipped to the plot area.
fn draw_line(
    buf: &mut [u8],
    fw: usize,
    x0: isize,
    y0: isize,
    x1: isize,
    y1: isize,
    clip_x: usize,
    clip_y: usize,
    clip_w: usize,
    clip_h: usize,
    color: [u8; 4],
) {
    let mut cx = x0;
    let mut cy = y0;
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx: isize = if x0 < x1 { 1 } else { -1 };
    let sy: isize = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;

    loop {
        if cx >= clip_x as isize
            && (cx as usize) < clip_x + clip_w
            && cy >= clip_y as isize
            && (cy as usize) <= clip_y + clip_h
        {
            set_pixel(buf, fw, cx as usize, cy as usize, color);
        }
        if cx == x1 && cy == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            cx += sx;
        }
        if e2 <= dx {
            err += dx;
            cy += sy;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lineplot_config_dimensions() {
        let cfg = LinePlotConfig::new(800, 600);
        assert_eq!(cfg.frame_width, 800);
        assert_eq!(cfg.frame_height, 600);
        assert!(cfg.plot_width() > 0);
        assert!(cfg.plot_height() > 0);
        assert!(cfg.plot_width() < 800);
        assert!(cfg.plot_height() < 600);
    }

    #[test]
    fn test_render_lineplot_buffer_size() {
        let cfg = LinePlotConfig::new(640, 480);
        let mut buf = Vec::new();
        let data: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();
        render_lineplot(
            &mut buf,
            &cfg,
            &data,
            10.0,
            (-1.0, 1.0),
            "u",
            0,
            1,
            ColorMap::BlueWhiteRed,
        );
        assert_eq!(buf.len(), 640 * 480 * 4);
    }

    #[test]
    fn test_render_lineplot_draws_nonblack_pixels() {
        let cfg = LinePlotConfig::new(640, 480);
        let mut buf = Vec::new();
        let data: Vec<f64> = (0..512).map(|i| (i as f64 * 0.05).sin() * 2.0).collect();
        render_lineplot(
            &mut buf,
            &cfg,
            &data,
            62.83,
            (-2.0, 2.0),
            "u",
            0,
            1,
            ColorMap::SolarWind,
        );
        // The line should have drawn bright-colored pixels in the plot area
        let plot_x0 = cfg.margin_left;
        let plot_y0 = cfg.margin_top;
        let pw = cfg.plot_width();
        let ph = cfg.plot_height();
        let mut bright_count = 0;
        for y in plot_y0..plot_y0 + ph {
            for x in plot_x0..plot_x0 + pw {
                let off = (y * cfg.frame_width + x) * 4;
                let lum = buf[off] as u32 + buf[off + 1] as u32 + buf[off + 2] as u32;
                if lum > 100 {
                    bright_count += 1;
                }
            }
        }
        assert!(
            bright_count > 100,
            "Line plot should draw bright pixels, got {bright_count}"
        );
    }

    #[test]
    fn test_format_axis_value() {
        assert_eq!(format_axis_value(0.0), "0");
        assert_eq!(format_axis_value(1.5), "1.50");
        assert_eq!(format_axis_value(42.3), "42.3");
        assert_eq!(format_axis_value(0.005), "5.0e-3");
        assert_eq!(format_axis_value(12345.0), "1.2e4");
    }
}
