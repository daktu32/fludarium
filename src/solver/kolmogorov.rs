use crate::state::{idx_inner, N};

/// Inject sinusoidal body force for Kolmogorov flow.
/// F_x = amplitude * sin(2π * wavenumber * j / N)
/// Applied as: vx += dt * F_x for all grid points.
pub(super) fn inject_body_force(
    vx: &mut [f64],
    dt: f64,
    amplitude: f64,
    wavenumber: f64,
    nx: usize,
) {
    let two_pi_k_over_n = 2.0 * std::f64::consts::PI * wavenumber / N as f64;
    for j in 1..(N - 1) {
        let force = amplitude * (two_pi_k_over_n * j as f64).sin();
        let dv = dt * force;
        for i in 0..nx {
            vx[idx_inner(i, j, nx)] += dv;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inject_body_force_sinusoidal_profile() {
        let nx = N;
        let mut vx = vec![0.0; nx * N];
        let dt = 0.05;
        let amplitude = 0.1;
        let wavenumber = 4.0;

        inject_body_force(&mut vx, dt, amplitude, wavenumber, nx);

        // Check that force varies sinusoidally with y
        let mid_x = nx / 2;

        // At j where sin is positive, vx should be positive
        // wavenumber=4 means 4 full cycles in N rows
        // sin(2π * 4 * j / N) at j = N/16 (quarter of first cycle) should be sin(π/2) = 1
        let j_peak = N / 16;
        let expected = dt * amplitude * (2.0 * std::f64::consts::PI * wavenumber * j_peak as f64 / N as f64).sin();
        let actual = vx[idx_inner(mid_x, j_peak, nx)];
        assert!((actual - expected).abs() < 1e-10,
            "vx at y={} should be {}, got {}", j_peak, expected, actual);

        // At j = N/8 (half cycle), sin should be 0
        let j_zero = N / 8;
        let val = vx[idx_inner(mid_x, j_zero, nx)];
        assert!(val.abs() < 1e-10, "vx at y={} should be ~0, got {}", j_zero, val);
    }

    #[test]
    fn test_inject_body_force_accumulates() {
        let nx = N;
        let mut vx = vec![1.0; nx * N];
        inject_body_force(&mut vx, 0.05, 0.1, 4.0, nx);
        // Interior values should be 1.0 + force contribution
        let j = N / 16;
        let val = vx[idx_inner(nx / 2, j, nx)];
        assert!(val > 1.0, "Force should add to existing velocity, got {}", val);
    }
}
