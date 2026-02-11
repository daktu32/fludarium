# fluvarium Progress

## Current Status: MVP — Thermal Convection Working

34 tests passing. Rayleigh-Benard convection visually confirmed in WezTerm.

## Completed

### Sub-issue 1: Project init
- `cargo init` with `icy_sixel`, `image`, `ctrlc` dependencies
- Module skeleton: `src/{solver,state,renderer,sixel}.rs`

### Sub-issue 2: SimState
- `SimState` struct with all grid fields + Xor128 PRNG
- `idx(x, y)` with mod wrap-around
- Initial conditions: linear gradient + large-scale sinusoidal perturbation to seed convection
- 11 tests

### Sub-issue 3: CFD primitives
- `set_bnd`: field_type 0/1 (Neumann), 2 (no-penetration wall), 3 (temperature Dirichlet)
- `lin_solve`: Gauss-Seidel iteration
- `diffuse`: implicit diffusion

### Sub-issue 4: Advection & projection
- `advect`: Semi-Lagrangian reverse trace + bilinear interpolation
- `project`: divergence → pressure solve → gradient subtraction

### Sub-issue 5: fluid_step & buoyancy
- `SolverParams`: visc=0.0001, diff=0.0001, dt=0.02, buoyancy=0.1, decay=0.999
- Buoyancy: `vy += buoyancy * (T - T_ambient)` where T_ambient=0.5
- Large-scale thermal perturbation injection (modes 1-4) each step
- Step order: diffuse vel → project → advect vel → buoyancy → **project** → diffuse temp → advect temp → perturbation → decay
- Key insight: projection must happen BEFORE temperature advection (divergence-free velocity required)
- Key insight: horizontal temperature variation is essential — uniform-in-x buoyancy is killed by projection

### Sub-issue 6: Renderer
- `temperature_to_rgba`: Blue(0.0)→Cyan→Green→Yellow→Red(1.0) colormap
- Color bar with tick marks on the right side
- `render`: temperature field + color bar → RGBA buffer with y-axis flip
- 8 tests

### Sub-issue 7: Sixel output
- `encode_sixel`: icy_sixel with Atkinson dithering
- `output_frame`: synchronized output (DEC 2026 BSU/ESU) + single-buffer write
- 3 tests

### Sub-issue 8: Main loop
- Alternate screen + cursor hide + Sixel scrolling disable
- `ctrlc` AtomicBool flag
- 5 steps/frame, 50ms frame cap (~20 FPS)
- 1 integration test

## Key Debugging History

1. **Buoyancy sign error**: `vy -=` → `vy +=` (hot fluid was sinking)
2. **Flickering**: Added synchronized output (DEC 2026), single-buffer writes, Atkinson dithering
3. **Temperature going all blue**: set_bnd Neumann BCs inside Gauss-Seidel iterations were erasing temperature boundaries → added field_type 3 (Dirichlet)
4. **Sharp red/blue interface**: Temperature advected with non-divergence-free velocity → moved projection before temperature advection
5. **No convection cells**: Uniform-in-x buoyancy killed by projection → switched from velocity noise to large-scale temperature perturbation at convection wavelengths

## Test Summary
- **34 tests, all passing**
- `cargo test` and `cargo build --release` both succeed

## Next Steps
- Fine-tune convection dynamics (plume count, speed, visual appeal)
- Performance optimization if needed
