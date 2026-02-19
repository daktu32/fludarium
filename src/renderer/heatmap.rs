// 2D heatmap renderer for channel grid playback.
//
// Renders a 2D field (nz x nx) as a colormapped image with axes.
// Used for PeriodicChannel data (Boussinesq convection, etc).

use super::color::{self, ColorMap, BAR_GAP, BAR_TOTAL, BAR_WIDTH, TICK_LEN, LABEL_GAP};
use super::font::{self, FONT_HEIGHT, FONT_WIDTH, STATUS_BAR_HEIGHT};

/// Layout configuration for the heatmap.
pub struct HeatmapConfig {
    pub frame_width: usize,
    pub frame_height: usize,
    pub margin_left: usize,
    pub margin_right: usize,
    pub margin_top: usize,
    pub margin_bottom: usize,
}

impl HeatmapConfig {
    pub fn new(win_width: usize, win_height: usize) -> Self {
        Self {
            frame_width: win_width,
            frame_height: win_height,
            margin_left: 45,
            margin_right: BAR_TOTAL + 5,
            margin_top: 35,
            margin_bottom: 25,
        }
    }

    pub fn plot_width(&self) -> usize {
        self.frame_width
            .saturating_sub(self.margin_left + self.margin_right)
            .max(4)
    }

    pub fn plot_height(&self) -> usize {
        self.frame_height
            .saturating_sub(self.margin_top + self.margin_bottom + STATUS_BAR_HEIGHT)
            .max(4)
    }
}

/// Render a 2D heatmap into an RGBA buffer.
///
/// `data` is row-major: data[j * nx + i] where j=0..nz-1, i=0..nx-1.
/// z-axis is vertical (bottom=0, top=nz-1), x-axis is horizontal.
pub fn render_heatmap(
    buf: &mut Vec<u8>,
    cfg: &HeatmapConfig,
    data: &[f64],
    nx: usize,
    nz: usize,
    domain: (f64, f64),
    value_range: (f64, f64),
    field_name: &str,
    field_index: usize,
    field_count: usize,
    colormap: ColorMap,
) {
    let fw = cfg.frame_width;
    let fh = cfg.frame_height;
    let total = fw * fh * 4;
    buf.resize(total, 0);

    // Dark background
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
    let (vmin, vmax) = value_range;
    let vrange = (vmax - vmin).max(1e-15);

    if nx == 0 || nz == 0 {
        return;
    }

    // Pre-compute interpolated value grid at screen resolution (once)
    let mut vals = vec![0.0f64; pw * ph];
    for py in 0..ph {
        let fz = (1.0 - (py as f64 + 0.5) / ph as f64) * (nz - 1) as f64;
        let iz0 = (fz as usize).min(nz - 2);
        let iz1 = iz0 + 1;
        let sz = fz - iz0 as f64;
        let sz_inv = 1.0 - sz;

        for px in 0..pw {
            let fx = (px as f64 + 0.5) / pw as f64 * nx as f64;
            let ix0f = fx.floor();
            let ix0 = ((ix0f as isize) % nx as isize + nx as isize) as usize % nx;
            let ix1 = (ix0 + 1) % nx;
            let sx = fx - ix0f;
            let sx_inv = 1.0 - sx;

            let v = data[iz0 * nx + ix0] * sx_inv * sz_inv
                + data[iz0 * nx + ix1] * sx * sz_inv
                + data[iz1 * nx + ix0] * sx_inv * sz
                + data[iz1 * nx + ix1] * sx * sz;

            vals[py * pw + px] = v;
        }
    }

    // Draw heatmap from pre-computed values (with single-pass contour overlay)
    let num_contours = 16.0_f64;
    let contour_interval = vrange / num_contours;
    let contour_threshold = 0.12;

    for py in 0..ph {
        let row = py * pw;
        for px in 0..pw {
            let val = vals[row + px];
            let t = ((val - vmin) / vrange).clamp(0.0, 1.0);
            let mut rgba = color::map_to_rgba(t, colormap);

            // Single-pass contour detection via modular arithmetic
            // (same approach as the built-in RB renderer)
            if contour_interval > 1e-15 {
                let frac = ((val - vmin) / contour_interval).rem_euclid(1.0);
                let line_dist = frac.min(1.0 - frac) * 2.0; // 0=on line, 1=between
                if line_dist < contour_threshold {
                    let brightness = 1.0 - (line_dist / contour_threshold);
                    let alpha = brightness * 0.6;
                    let inv = 1.0 - alpha;
                    rgba[0] = (rgba[0] as f64 * inv + 255.0 * alpha) as u8;
                    rgba[1] = (rgba[1] as f64 * inv + 255.0 * alpha) as u8;
                    rgba[2] = (rgba[2] as f64 * inv + 255.0 * alpha) as u8;
                }
            }

            let off = ((y0 + py) * fw + (x0 + px)) * 4;
            if off + 3 < buf.len() {
                buf[off] = rgba[0];
                buf[off + 1] = rgba[1];
                buf[off + 2] = rgba[2];
                buf[off + 3] = rgba[3];
            }
        }
    }

    // Axis frame
    let axis_color: [u8; 3] = [0x44, 0x4C, 0x5C];
    // Left edge
    for py in y0..=y0 + ph {
        set_pixel_rgb(buf, fw, x0, py, axis_color);
    }
    // Bottom edge
    for px in x0..=x0 + pw {
        set_pixel_rgb(buf, fw, px, y0 + ph, axis_color);
    }
    // Right edge
    for py in y0..=y0 + ph {
        set_pixel_rgb(buf, fw, x0 + pw, py, axis_color);
    }
    // Top edge
    for px in x0..=x0 + pw {
        set_pixel_rgb(buf, fw, px, y0, axis_color);
    }

    // Axis tick labels
    let label_color: [u8; 3] = [0x66, 0x6C, 0x7C];
    let (lx_domain, lz_domain) = domain;

    // Z-axis labels (left side): 5 ticks from bottom (0) to top (lz)
    for i in 0..=4 {
        let frac = i as f64 / 4.0;
        let val = frac * lz_domain;
        let py = y0 + ph - (ph * i / 4); // bottom to top
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

    // X-axis labels (bottom): 5 ticks from left (0) to right (lx)
    for i in 0..=4 {
        let frac = i as f64 / 4.0;
        let val = frac * lx_domain;
        let px = x0 + pw * i / 4;
        let label = format_axis_value(val);
        let label_w = label.len() * (FONT_WIDTH + 1);
        let lx = px.saturating_sub(label_w / 2);
        let ly = y0 + ph + 6;
        font::draw_text(buf, fw, lx, ly, &label, label_color);
    }

    // Field badge
    font::render_field_badge(buf, fw, y0 + ph, field_name, field_index, field_count);

    // Color bar (right side of plot)
    let bar_x = x0 + pw + BAR_GAP;
    let bar_h = ph;

    // Bar gradient
    for py in 0..bar_h {
        let t = 1.0 - py as f64 / bar_h.max(1) as f64;
        let rgba = color::map_to_rgba(t, colormap);
        for bx in 0..BAR_WIDTH {
            let sx = bar_x + bx;
            let sy = y0 + py;
            let off = (sy * fw + sx) * 4;
            if off + 3 < buf.len() {
                buf[off] = rgba[0];
                buf[off + 1] = rgba[1];
                buf[off + 2] = rgba[2];
                buf[off + 3] = 255;
            }
        }
    }

    // Tick marks and labels
    let tick_x = bar_x + BAR_WIDTH;
    let label_x = tick_x + TICK_LEN + LABEL_GAP;
    let tick_label_color: [u8; 3] = [0x88, 0x88, 0x88];

    let tick_fracs = [1.0, 0.75, 0.5, 0.25, 0.0];
    for (i, &tf) in tick_fracs.iter().enumerate() {
        let py = y0 + i * bar_h.max(1) / 4;
        let actual_val = vmin + (vmax - vmin) * tf;
        let label = format_axis_value(actual_val);

        // Tick mark
        for dy in 0..2usize {
            let yy = (py + dy).min(y0 + bar_h - 1);
            for tx in 0..TICK_LEN {
                let off = (yy * fw + tick_x + tx) * 4;
                if off + 3 < buf.len() {
                    buf[off] = 0xFF;
                    buf[off + 1] = 0xFF;
                    buf[off + 2] = 0xFF;
                    buf[off + 3] = 0xFF;
                }
            }
        }

        // Label
        let ly = if py >= FONT_HEIGHT / 2 {
            py - FONT_HEIGHT / 2
        } else {
            0
        };
        font::draw_text(buf, fw, label_x, ly, &label, tick_label_color);
    }
}

// Meteor color palette (same as built-in model renderer).
const METEOR_COLORS: [[f64; 3]; 4] = [
    [255.0, 240.0, 200.0], // warm white (newest)
    [180.0, 220.0, 255.0], // light cyan
    [80.0, 140.0, 255.0],  // blue
    [30.0, 50.0, 120.0],   // deep blue (oldest)
];

#[inline]
fn meteor_color(t: f64) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let seg = t * (METEOR_COLORS.len() - 1) as f64;
    let i = (seg as usize).min(METEOR_COLORS.len() - 2);
    let f = seg - i as f64;
    let c0 = &METEOR_COLORS[METEOR_COLORS.len() - 1 - i];
    let c1 = &METEOR_COLORS[(METEOR_COLORS.len() - 2).saturating_sub(i)];
    [
        c0[0] + f * (c1[0] - c0[0]),
        c0[1] + f * (c1[1] - c0[1]),
        c0[2] + f * (c1[2] - c0[2]),
    ]
}

/// Render channel particles on top of the heatmap buffer.
///
/// Draws meteor-colored trail lines (oldest→newest) and bright heads with glow.
/// Coordinates: x in [0, lx) periodic, z in [0, lz] (z=0 at bottom of plot).
pub fn render_channel_particles(
    buf: &mut [u8],
    cfg: &HeatmapConfig,
    px: &[f64],
    pz: &[f64],
    trails: &[(&[f64], &[f64])],
    lx: f64,
    lz: f64,
) {
    let fw = cfg.frame_width;
    let pw = cfg.plot_width();
    let ph = cfg.plot_height();
    let x0 = cfg.margin_left;
    let y0 = cfg.margin_top;

    let to_screen = |sx: f64, sz: f64| -> (isize, isize) {
        let scr_x = x0 as f64 + (sx / lx) * pw as f64;
        let scr_y = y0 as f64 + (1.0 - sz / lz) * ph as f64;
        (scr_x as isize, scr_y as isize)
    };

    let in_plot = |sx: isize, sy: isize| -> bool {
        sx >= x0 as isize
            && sx < (x0 + pw) as isize
            && sy >= y0 as isize
            && sy < (y0 + ph) as isize
    };

    // Additive blend (like the built-in model's meteor mode)
    let additive = |buf: &mut [u8], off: usize, r: f64, g: f64, b: f64| {
        if off + 3 >= buf.len() {
            return;
        }
        buf[off] = (buf[off] as f64 + r).min(255.0) as u8;
        buf[off + 1] = (buf[off + 1] as f64 + g).min(255.0) as u8;
        buf[off + 2] = (buf[off + 2] as f64 + b).min(255.0) as u8;
    };

    let trail_count = trails.len();
    if trail_count < 2 {
        return;
    }

    // Draw trail segments: line from trails[i] → trails[i+1]
    for ti in 0..(trail_count - 1) {
        let frac = (ti + 1) as f64 / trail_count as f64;
        let alpha = frac.powf(2.5);
        let color = meteor_color(frac);

        let (xs0, zs0) = trails[ti];
        let (xs1, zs1) = trails[ti + 1];

        for p in 0..xs0.len() {
            let (sx0, sy0) = to_screen(xs0[p], zs0[p]);
            let (sx1, sy1) = to_screen(xs1[p], zs1[p]);

            // Bresenham line between consecutive trail points
            let dx = (sx1 - sx0).abs();
            let dy = -(sy1 - sy0).abs();
            let step_x: isize = if sx0 < sx1 { 1 } else { -1 };
            let step_y: isize = if sy0 < sy1 { 1 } else { -1 };
            let mut err = dx + dy;
            let mut cx = sx0;
            let mut cy = sy0;

            // Skip lines that wrap around the periodic boundary
            if dx > pw as isize / 2 {
                continue;
            }

            let max_steps = (dx.abs() + dy.abs() + 1) as usize;
            for _ in 0..max_steps {
                if in_plot(cx, cy) {
                    let off = (cy as usize * fw + cx as usize) * 4;
                    additive(
                        buf,
                        off,
                        color[0] * alpha,
                        color[1] * alpha,
                        color[2] * alpha,
                    );
                }
                if cx == sx1 && cy == sy1 {
                    break;
                }
                let e2 = 2 * err;
                if e2 >= dy {
                    err += dy;
                    cx += step_x;
                }
                if e2 <= dx {
                    err += dx;
                    cy += step_y;
                }
            }
        }
    }

    // Draw particle heads with 5x5 soft glow
    const GLOW: &[(isize, isize, f64)] = &[
        (0, 0, 1.0),
        (-1, 0, 0.6),
        (1, 0, 0.6),
        (0, -1, 0.6),
        (0, 1, 0.6),
        (-1, -1, 0.3),
        (1, -1, 0.3),
        (-1, 1, 0.3),
        (1, 1, 0.3),
        (-2, 0, 0.15),
        (2, 0, 0.15),
        (0, -2, 0.15),
        (0, 2, 0.15),
    ];

    for i in 0..px.len() {
        let (sx, sy) = to_screen(px[i], pz[i]);
        for &(dx, dy, intensity) in GLOW {
            let gx = sx + dx;
            let gy = sy + dy;
            if in_plot(gx, gy) {
                let off = (gy as usize * fw + gx as usize) * 4;
                additive(buf, off, 255.0 * intensity, 240.0 * intensity, 200.0 * intensity);
            }
        }
    }
}

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

#[inline]
fn set_pixel_rgb(buf: &mut [u8], fw: usize, x: usize, y: usize, rgb: [u8; 3]) {
    let off = (y * fw + x) * 4;
    if off + 3 < buf.len() {
        buf[off] = rgb[0];
        buf[off + 1] = rgb[1];
        buf[off + 2] = rgb[2];
        buf[off + 3] = 0xFF;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heatmap_config_dimensions() {
        let cfg = HeatmapConfig::new(800, 600);
        assert!(cfg.plot_width() > 0);
        assert!(cfg.plot_height() > 0);
        assert!(cfg.plot_width() < 800);
        assert!(cfg.plot_height() < 600);
    }

    #[test]
    fn test_render_heatmap_buffer_size() {
        let cfg = HeatmapConfig::new(640, 480);
        let mut buf = Vec::new();
        let data: Vec<f64> = (0..128 * 64).map(|i| (i as f64 * 0.01).sin()).collect();
        render_heatmap(
            &mut buf,
            &cfg,
            &data,
            128,
            64,
            (2.83, 1.0),
            (-1.0, 1.0),
            "theta",
            0,
            3,
            ColorMap::BlueWhiteRed,
        );
        assert_eq!(buf.len(), 640 * 480 * 4);
    }

    #[test]
    fn test_render_heatmap_draws_colored_pixels() {
        let cfg = HeatmapConfig::new(640, 480);
        let mut buf = Vec::new();
        // Create a gradient pattern
        let nx = 64;
        let nz = 32;
        let data: Vec<f64> = (0..nz * nx)
            .map(|i| {
                let j = i / nx;
                j as f64 / nz as f64
            })
            .collect();
        render_heatmap(
            &mut buf,
            &cfg,
            &data,
            nx,
            nz,
            (4.0, 2.0),
            (0.0, 1.0),
            "theta",
            0,
            1,
            ColorMap::TokyoNight,
        );
        // Check that the plot area has colored pixels (not just background)
        let x0 = cfg.margin_left;
        let y0 = cfg.margin_top;
        let pw = cfg.plot_width();
        let ph = cfg.plot_height();
        let mut bright_count = 0;
        for y in y0..y0 + ph {
            for x in x0..x0 + pw {
                let off = (y * cfg.frame_width + x) * 4;
                let lum = buf[off] as u32 + buf[off + 1] as u32 + buf[off + 2] as u32;
                if lum > 50 {
                    bright_count += 1;
                }
            }
        }
        assert!(
            bright_count > 100,
            "Heatmap should draw colored pixels, got {bright_count}"
        );
    }
}
