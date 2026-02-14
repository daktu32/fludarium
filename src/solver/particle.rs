use crate::state::{idx, N};

/// Bilinearly interpolate velocity (vx, vy) at grid-relative position.
/// `i0`, `j0` are the lower-left cell indices; `sx`, `sy` are fractional offsets in [0, 1].
pub(super) fn interpolate_velocity(
    vx: &[f64], vy: &[f64],
    i0: i32, j0: i32,
    sx: f64, sy: f64,
    nx: usize,
) -> (f64, f64) {
    let i1 = i0 + 1;
    let j1 = j0 + 1;
    let vx_val = (1.0 - sx) * (1.0 - sy) * vx[idx(i0, j0, nx)]
        + sx * (1.0 - sy) * vx[idx(i1, j0, nx)]
        + (1.0 - sx) * sy * vx[idx(i0, j1, nx)]
        + sx * sy * vx[idx(i1, j1, nx)];
    let vy_val = (1.0 - sx) * (1.0 - sy) * vy[idx(i0, j0, nx)]
        + sx * (1.0 - sy) * vy[idx(i1, j0, nx)]
        + (1.0 - sx) * sy * vy[idx(i0, j1, nx)]
        + sx * sy * vy[idx(i1, j1, nx)];
    (vx_val, vy_val)
}

/// Advect particles through the velocity field using bilinear interpolation.
/// X wraps (periodic), Y reflects at walls.
pub(super) fn advect_particles(state: &mut crate::state::SimState, dt: f64) {
    let nx = state.nx;
    let dt0 = dt * (N - 2) as f64;
    let n_f = N as f64;
    let nx_f = nx as f64;

    for p in 0..state.particles_x.len() {
        let px = state.particles_x[p];
        let py = state.particles_y[p];

        // Bilinear interpolation of velocity at particle position
        let i0 = px.floor() as i32;
        let j0 = py.floor().max(0.0).min(n_f - 2.0) as i32;
        let sx = px - px.floor();
        let sy = py - j0 as f64;
        let (vx_interp, vy_interp) = interpolate_velocity(&state.vx, &state.vy, i0, j0, sx, sy, nx);

        // Move particle forward
        let new_x = px + dt0 * vx_interp;
        let mut new_y = py + dt0 * vy_interp;

        // Y: ping-pong reflect within interior [y_min, y_max].
        // Keep particles outside the 2-row Dirichlet boundary zone where
        // no-slip makes velocity = 0 and particles would get trapped.
        let y_min = 2.0;
        let y_max = n_f - 3.0;
        let y_range = y_max - y_min;
        if new_y < y_min || new_y > y_max {
            let mut t = (new_y - y_min) % (2.0 * y_range);
            if t < 0.0 {
                t += 2.0 * y_range;
            }
            new_y = if t <= y_range {
                y_min + t
            } else {
                y_max - (t - y_range)
            };
        }
        new_y = new_y.clamp(y_min, y_max);

        // X: wrap around (periodic)
        let new_x = ((new_x % nx_f) + nx_f) % nx_f;

        state.particles_x[p] = new_x;
        state.particles_y[p] = new_y;
    }
}

/// Advect particles for Karman flow.
/// Particles that exit right or enter cylinder are respawned at left.
pub(super) fn advect_particles_karman(state: &mut crate::state::SimState, dt: f64) {
    let nx = state.nx;
    let dt0 = dt * (N - 2) as f64;
    let n_f = N as f64;
    let nx_f = nx as f64;

    for p in 0..state.particles_x.len() {
        let px = state.particles_x[p];
        let py = state.particles_y[p];

        // Bilinear interpolation of velocity
        let i0 = px.floor().max(0.0).min(nx_f - 2.0) as i32;
        let j0 = py.floor().max(0.0).min(n_f - 2.0) as i32;
        let sx = px - i0 as f64;
        let sy = py - j0 as f64;
        let (vx_interp, vy_interp) = interpolate_velocity(&state.vx, &state.vy, i0, j0, sx, sy, nx);

        let new_x = px + dt0 * vx_interp;
        let new_y = (py + dt0 * vy_interp).clamp(2.0, n_f - 3.0);

        // Check if particle needs respawn (out of domain or near cylinder)
        let needs_respawn = new_x >= nx_f - 1.0 || new_x < 0.0 || {
            if let Some((cx, cy, r)) = state.cylinder {
                let dx = new_x - cx;
                let dy = new_y - cy;
                dx * dx + dy * dy < (r + 1.0) * (r + 1.0)
            } else {
                false
            }
        };

        if needs_respawn {
            // Respawn at left boundary with random y
            state.particles_x[p] = 2.0 + state.rng.next_f64().abs() * 2.0;
            state.particles_y[p] = 2.0 + (state.rng.next_f64() + 1.0) * 0.5 * (n_f - 5.0);
        } else {
            state.particles_x[p] = new_x;
            state.particles_y[p] = new_y;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::{idx, SimState, N};
    use crate::solver::{fluid_step, SolverParams};

    #[test]
    fn test_interpolate_velocity_uniform() {
        // Uniform velocity field -> interpolation should return exact value
        let n = 10;
        let len = n * n;
        let vx = vec![3.0; len];
        let vy = vec![-1.5; len];
        let (vx_val, vy_val) = interpolate_velocity(&vx, &vy, 2, 3, 0.5, 0.5, n);
        assert!((vx_val - 3.0).abs() < 1e-10);
        assert!((vy_val - (-1.5)).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate_velocity_corner_weights() {
        // Place known values at 4 corners, verify bilinear blend
        let n = 10;
        let len = n * n;
        let mut vx = vec![0.0; len];
        // Set corners of cell (2,3): (2,3)=1, (3,3)=2, (2,4)=3, (3,4)=4
        vx[idx(2, 3, n)] = 1.0;
        vx[idx(3, 3, n)] = 2.0;
        vx[idx(2, 4, n)] = 3.0;
        vx[idx(3, 4, n)] = 4.0;
        let vy = vec![0.0; len];
        // At center (sx=0.5, sy=0.5): (1+2+3+4)/4 = 2.5
        let (vx_val, _) = interpolate_velocity(&vx, &vy, 2, 3, 0.5, 0.5, n);
        assert!((vx_val - 2.5).abs() < 1e-10, "Expected 2.5, got {vx_val}");
        // At lower-left corner (sx=0, sy=0): should be 1.0
        let (vx_val, _) = interpolate_velocity(&vx, &vy, 2, 3, 0.0, 0.0, n);
        assert!((vx_val - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_advect_particles_zero_velocity() {
        let mut state = SimState::new(400, 0.15, N);
        // Zero out velocity
        state.vx.fill(0.0);
        state.vy.fill(0.0);
        let orig_x = state.particles_x.clone();
        let orig_y = state.particles_y.clone();

        advect_particles(&mut state, 0.02);

        // Particles should not move
        for i in 0..state.particles_x.len() {
            assert!(
                (state.particles_x[i] - orig_x[i]).abs() < 1e-10,
                "Particle {} x moved with zero velocity",
                i
            );
            assert!(
                (state.particles_y[i] - orig_y[i]).abs() < 1e-10,
                "Particle {} y moved with zero velocity",
                i
            );
        }
    }

    #[test]
    fn test_advect_particles_uniform_velocity() {
        let mut state = SimState::new(400, 0.15, N);
        // Uniform rightward velocity
        state.vx.fill(0.01);
        state.vy.fill(0.0);
        let orig_x = state.particles_x.clone();

        advect_particles(&mut state, 0.02);

        // All particles should have moved right
        let dt0 = 0.02 * (N - 2) as f64;
        let expected_dx = dt0 * 0.01;
        for i in 0..state.particles_x.len() {
            let dx = ((state.particles_x[i] - orig_x[i]) + N as f64) % N as f64;
            assert!(
                (dx - expected_dx).abs() < 1e-6 || (dx - expected_dx + N as f64).abs() < 1e-6,
                "Particle {} dx={} expected ~{}",
                i, dx, expected_dx
            );
        }
    }

    #[test]
    fn test_advect_particles_stay_in_domain() {
        let mut state = SimState::new(400, 0.15, N);
        let params = SolverParams::default();

        // Run 50 steps with active flow
        for _ in 0..50 {
            fluid_step(&mut state, &params);
        }

        // All particles should still be in domain
        let n_f = N as f64;
        for i in 0..state.particles_x.len() {
            let px = state.particles_x[i];
            let py = state.particles_y[i];
            assert!(px >= 0.0 && px < n_f, "Particle {} x out of range: {}", i, px);
            assert!(py >= 2.0 && py <= n_f - 3.0, "Particle {} y out of range: {}", i, py);
        }
    }
}
