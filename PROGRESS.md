# fluvarium Progress

## Current Status: Dual-Model Fluid Simulator

103 tests passing. 2-thread pipeline (physics + render/display). Dual simulation models (Rayleigh-Benard convection + Kármán vortex street). N=80 grid with aspect-scaled NX for Kármán. Real-time parameter tuning via overlay panel. No external config files — all defaults in code.

## Completed

### Sub-issue 1: Project init
- `cargo init` with `icy_sixel`, `ctrlc`, `minifb` dependencies
- Module skeleton: `src/{solver,state,renderer,sixel}.rs`

### Sub-issue 2: SimState
- `SimState` struct with all grid fields + Xor128 PRNG
- `idx(x, y, nx)` with mod wrap-around, variable grid width
- `FrameSnapshot` for decoupling simulation state from rendering (includes velocity, trails, cylinder geometry)
- `FluidModel` enum: `RayleighBenard` / `KarmanVortex`
- Initial conditions: Gaussian hot spot at bottom center (RB), uniform inflow + cylinder obstacle (Kármán)
- Particle trail ring buffer (`TRAIL_LEN=8`) for trajectory visualization
- Fractional cylinder mask with smooth anti-aliased edges

### Sub-issue 3: CFD primitives
- `BoundaryConfig` enum dispatches per-model boundary conditions
- `set_bnd_rb`: periodic X, field_type 0/1 (Neumann), 2 (no-penetration wall), 3 (temperature Dirichlet with Gaussian profile)
- `set_bnd_karman`: Left Dirichlet inflow, right zero-gradient outflow, top/bottom no-slip walls
- `lin_solve`: Gauss-Seidel iteration
- `diffuse`: implicit diffusion

### Sub-issue 4: Advection & projection
- `advect`: Semi-Lagrangian reverse trace + bilinear interpolation, X clamping for Kármán
- `project`: divergence → pressure solve → gradient subtraction

### Sub-issue 5: fluid_step & buoyancy
- `SolverParams` with `default()` (RB) and `default_karman()` — all params in code
- Buoyancy: `vy += dt * buoyancy * (T - T_ambient)` where T_ambient=bottom_base
- Localized Gaussian heat source (`inject_heat_source`) with volumetric heating + Newtonian cooling
- Particle advection with bilinear velocity interpolation, periodic X wrap, ping-pong Y reflection
- `fluid_step_karman`: inflow injection, cylinder mask damping, dye tracer, vorticity confinement, wake perturbation
- Geometry-based particle respawn (distance check instead of mask threshold) to prevent surface sticking

### Sub-issue 6: Renderer
- Tokyo Night colormap (navy→blue→purple→pink→orange)
- Color bar with tick marks on the right side
- Adaptive contrast particles (bright on dark, dark on bright backgrounds)
- Dynamic `RenderConfig::fit()` for arbitrary pixel dimensions and configurable tile count
- Kármán: stretched to window, smooth anti-aliased cylinder rendering, optional vorticity visualization, particle trails
- **Status bar**: 5x7 bitmap font rendering current parameters
- **Font system**: nearest-neighbor resize for overlay text

### Sub-issue 7: Sixel output (test-only)
- `encode_sixel`: icy_sixel (test-only)
- `SixelEncoder`: custom encoder with fixed 64-color palette + 32KB RGB→palette LUT (test-only)

### Sub-issue 8: Main loop — minifb native window
- **Physics thread**: `fluid_step` / `fluid_step_karman` (configurable steps/frame) → `FrameSnapshot` via `sync_channel(1)`
- **Main thread**: `render()` → `rgba_to_argb()` → `window.update_with_buffer()` at 60fps target
- Ctrl+C handler with clean shutdown
- **Resizable window** with dynamic re-render on size change
- FPS counter in title bar
- **M key**: model switch with per-model parameter preservation (`rb_params` / `karman_params`)
- **V key**: toggle vorticity visualization (Kármán mode)

### Sub-issue 9: Interactive Overlay Parameter Panel
- **overlay.rs**: btop-style semi-transparent panel with darken background
  - `OverlayState` (visible, selected), `ParamDef` with get/set function pointers
  - RB: 7 params (visc, diff, dt, buoyancy, source_strength, cool_rate, bottom_base)
  - Kármán: 5 params (visc, diff, dt, inflow_vel, confinement)
- **Keyboard controls**: Space=toggle, Up/Down=navigate, Left/Right=adjust, Comma/Period=fine, R=reset, Escape=close/quit
- **Real-time tuning**: `mpsc::channel` sends updated `SolverParams` to physics thread

### Config simplification
- Removed YAML config system (`config.rs`, `serde`/`serde_yaml` dependencies)
- All defaults managed in code: `SolverParams::default()`, `default_karman()`, `PARAM_DEFS_*`
- Window size (1280×640), tiles, particles hardcoded in `main.rs`
- Startup model: Kármán vortex (visc=0.015, diff=0.003, dt=0.06, confinement=3.0, cylinder at vertical center)

## Performance Evolution

| Milestone | fps | Bottleneck |
|-----------|-----|------------|
| MVP single-thread (Sixel) | ~7 | icy_sixel encoding (130ms) |
| Custom Sixel encoder | ~7 | Physics N=256 (148ms) |
| N=256, dt=0.002 | ~14 | Physics (62ms) |
| N=128, half-size (Sixel) | ~55 | None (all stages <16ms) |
| N=128, full-width × half-height (Sixel) | ~34 | write (23ms) |
| **minifb native window** | **~60** | **None (vsync-limited)** |

## Key Debugging History

1. **Buoyancy sign error**: `vy -=` → `vy +=` (hot fluid was sinking)
2. **Flickering**: Added synchronized output (DEC 2026), single-buffer writes
3. **Temperature going all blue**: set_bnd Neumann BCs inside Gauss-Seidel iterations were erasing temperature boundaries → added field_type 3 (Dirichlet)
4. **Sharp red/blue interface**: Temperature advected with non-divergence-free velocity → moved projection before temperature advection
5. **No convection cells**: Uniform-in-x buoyancy killed by projection → switched from velocity noise to large-scale temperature perturbation at convection wavelengths
6. **Particle warping**: `drain_latest()` skipped intermediate frames → switched to `recv()` everywhere
7. **Ctrl+C deadlock**: 3-thread pipeline deadlocked on shutdown (blocked senders) → `drop(sixel_rx)` before joining threads
8. **Turbulent plume**: High buoyancy (60-200) caused Ra≈60,000+ with dominant numerical viscosity → physicist-tuned params (visc=0.008, diff=0.002, buoyancy=8.0) for clean convection
9. **Sixel bottleneck**: Sixel encode+write limited to ~34fps → switched to minifb native window for 60fps
10. **Particle sticking on cylinder**: mask threshold (>0.5) missed smooth edge zone → switched to geometry-based distance check with 1-cell margin
11. **Config/default mismatch**: YAML config values diverged from `default_karman()` → removed config system, single source of truth in code

## Dependency Cleanup
- Removed `libc` (terminal ioctl no longer needed with minifb)
- Removed `image` (unused)
- Removed `serde` + `serde_yaml` (YAML config removed)
- Sixel encoder gated behind `#[cfg(test)]`

## Test Summary
- **103 tests, all passing** (1 ignored: diagnostic)
- `cargo test` and `cargo build --release` both succeed with 0 warnings

## Next Steps
- Fine-tune convection dynamics (plume count, speed, visual appeal)
- Issue #10: Simulation presets
