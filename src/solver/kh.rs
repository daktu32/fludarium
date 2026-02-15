use crate::state::{idx_inner, N};

use super::SolverParams;

/// Weakly relax velocity toward the target shear profile near walls only.
/// Prevents numerical diffusion from erasing the driving shear, but leaves
/// the interior free so KH instability can grow at the interface.
pub fn reinject_shear(state: &mut crate::state::SimState, params: &SolverParams) {
    let nx = state.nx;
    let center_y = (N / 2) as f64;
    let delta = params.shear_thickness;
    let alpha = params.shear_relax;
    if alpha <= 0.0 {
        return;
    }
    // Only relax near walls (bottom 6 rows and top 6 rows).
    // Interior is left free for instability growth.
    let wall_rows = 6usize;
    for j in 0..N {
        // Distance from nearest wall (0 at wall, increases inward)
        let dist_from_wall = j.min(N - 1 - j);
        if dist_from_wall >= wall_rows {
            continue;
        }
        // Blend factor: strongest at wall, fading inward
        let fade = 1.0 - dist_from_wall as f64 / wall_rows as f64;
        let blend = alpha * params.dt * fade;
        let y_norm = (j as f64 - center_y) / delta;
        let target = params.shear_velocity * y_norm.tanh();
        for i in 0..nx {
            let ii = idx_inner(i, j, nx);
            state.vx[ii] += blend * (target - state.vx[ii]);
        }
    }
}

/// Inject dye tracer at left boundary columns matching the shear interface.
/// Not used in fluid_step_kh (periodic boundary needs no reinjection),
/// but kept for tests and potential future non-periodic variants.
#[cfg(test)]
pub fn inject_dye_kh(state: &mut crate::state::SimState) {
    let nx = state.nx;
    let center_y = (N / 2) as f64;
    let dye_width = 3.0;
    for j in 0..N {
        let y_norm = (j as f64 - center_y) / dye_width;
        let dye = 0.5 + 0.5 * y_norm.tanh();
        state.temperature[idx_inner(0, j, nx)] = dye;
        state.temperature[idx_inner(1, j, nx)] = dye;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::SolverParams;
    use crate::state::{idx, SimState, N};

    #[test]
    fn test_reinject_shear_restores_profile() {
        let params = SolverParams::default_kh();
        let mut state = SimState::new_kh(10, &params, N);
        // Zero out velocity
        for v in state.vx.iter_mut() {
            *v = 0.0;
        }
        // Reinject multiple times
        for _ in 0..100 {
            reinject_shear(&mut state, &params);
        }
        // Top should be positive, bottom should be negative
        let top_vx = state.vx[idx((N / 2) as i32, (N - 5) as i32, N)];
        let bot_vx = state.vx[idx((N / 2) as i32, 5, N)];
        assert!(top_vx > 0.0, "Top should have positive vx after reinjection, got {}", top_vx);
        assert!(bot_vx < 0.0, "Bottom should have negative vx after reinjection, got {}", bot_vx);
    }

    #[test]
    fn test_reinject_shear_zero_relax_noop() {
        let mut params = SolverParams::default_kh();
        params.shear_relax = 0.0;
        let mut state = SimState::new_kh(10, &params, N);
        let vx_before = state.vx.clone();
        reinject_shear(&mut state, &params);
        assert_eq!(state.vx, vx_before);
    }

    #[test]
    fn test_inject_dye_kh_interface() {
        let params = SolverParams::default_kh();
        let mut state = SimState::new_kh(10, &params, N);
        inject_dye_kh(&mut state);
        // Bottom-left should be near 0
        assert!(state.temperature[idx(0, 2, N)] < 0.2, "Bottom dye should be near 0");
        // Top-left should be near 1
        assert!(state.temperature[idx(0, (N - 3) as i32, N)] > 0.8, "Top dye should be near 1");
    }
}
