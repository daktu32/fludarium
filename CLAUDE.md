# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fluvarium is a fluid dynamics visualizer supporting multiple simulation models: Rayleigh-Benard convection (thermal plumes) and Karman vortex street (flow past a cylinder). Renders via minifb native window with real-time parameter tuning, or headless in terminal via iTerm2 Graphics Protocol (`--headless`). Written in Rust (edition 2024), single crate, no external C dependencies.

## Build & Test Commands

```bash
cargo build              # Debug build
cargo build --release    # Release build (recommended for running)
cargo test               # Run all tests (~108 passing, 1 ignored)
cargo test solver        # Run only solver module tests
cargo test state         # Run only state module tests
cargo test renderer      # Run only renderer module tests
cargo test overlay       # Run only overlay module tests
cargo test sixel         # Run only sixel module tests
cargo test test_name     # Run a single test by name
cargo run --release      # Run the simulation (GUI mode)
cargo run --release -- --headless  # Headless terminal mode (iTerm2/WezTerm)
```

## Architecture

2-thread pipeline with inter-thread communication:

```
Physics thread: SimState → fluid_step() → FrameSnapshot → sync_channel(1) ─┐
                     ↑                                                       │
              mpsc::channel (SolverParams updates from UI)                   │
                                                                             ↓
GUI mode:     keyboard input → FrameSnapshot → render() → overlay → rgba_to_argb() → minifb window
Headless mode:                  FrameSnapshot → render() → RGBA → PNG → base64 → iTerm2 escape → stdout
```

### Dual Simulation Model

The `FluidModel` enum (`state.rs`) selects between:
- **RayleighBenard**: Thermal convection with periodic X boundaries, hot bottom (Gaussian heat source), cold top. Grid is N×N, rendered with configurable tile count for horizontal repetition.
- **KarmanVortex**: Uniform inflow from left past a cylinder obstacle, outflow right. Grid is NX×N where NX scales with window aspect ratio. Supports dye and vorticity visualization, particle trails, and vorticity confinement.

Each model has its own `fluid_step` / `fluid_step_karman`, `set_bnd_rb` / `set_bnd_karman`, `SimState::new` / `SimState::new_karman`, and `SolverParams::default` / `SolverParams::default_karman`. The user can switch models at runtime with the M key; parameters are preserved per-model.

### Module Responsibilities

- **`state.rs`** — `FluidModel` enum, `SimState` struct (velocity fields `vx/vy/vx0/vy0`, `temperature`, scratch buffers, `Xor128` PRNG, particles, optional cylinder `mask`/`cylinder` geometry, particle trail ring buffer). Grid height is always `N`, width is `nx` (= N for RB, aspect-scaled for Karman). `idx(x, y, nx)` handles 2D→1D with mod wrapping. `FrameSnapshot` decouples simulation from rendering (includes temperature, velocity, particles, trails, cylinder info).
- **`solver.rs`** — Jos Stam "Stable Fluids" CFD. `BoundaryConfig` enum dispatches boundary conditions per model. Core functions: `diffuse`, `advect` (Semi-Lagrangian), `project` (pressure projection). `SolverParams` holds all tuning constants with `default()` (RB) and `default_karman()`. Karman-specific: `apply_cylinder_mask`, `vorticity_confinement`, `inject_dye_at_inflow`.
- **`overlay.rs`** — btop-style semi-transparent parameter panel. `OverlayState` (visible, selected), `ParamDef` with get/set function pointers. Model-specific param sets (`PARAM_DEFS_RB` with 7 params, `PARAM_DEFS_KARMAN` with 5 params). Teal gradient gauge bars, keyboard controls.
- **`renderer.rs`** — Maps fields to RGBA via 5-stop Tokyo Night colormap (navy→blue→purple→pink→orange). `RenderConfig::fit(w, h, tiles, sim_nx)` computes layout. RB: tiled horizontally with color bar. Karman: stretched to window, smooth cylinder rendering, optional vorticity field, particle trails. 5×7 bitmap font system with `draw_text` (1×) and `draw_text_sized` (scaled) for overlay/status bar.
- **`iterm2.rs`** — `Iterm2Encoder` for headless terminal rendering. RGBA → RGB → PNG (fast compression) → base64 → iTerm2 `\x1b]1337;File=...` escape sequence. Reusable buffers for zero per-frame allocation.
- **`sixel.rs`** — Custom `SixelEncoder` with fixed palette + LUT, and `encode_sixel` via icy_sixel. Both gated behind `#[cfg(test)]`.
- **`main.rs`** — Creates window with hardcoded defaults, spawns physics thread with `mpsc` channels for parameter/model-reset updates. Per-model params (`rb_params`/`karman_params`) preserved across M-key switches. Keyboard: Space=overlay, M=model switch, V=vorticity toggle, arrows/comma/period=param adjust, R=reset, Escape=close/quit.

### Coordinate System

- `y=0` is **bottom** (hot), `y=N-1` is **top** (cold) in simulation space
- Renderer flips Y for screen output (screen top = cold, screen bottom = hot)
- RB: X-axis wraps (periodic boundary), Y-axis has wall boundaries
- Karman: Left=inflow (Dirichlet), Right=outflow (zero-gradient), Top/Bottom=no-slip walls

### Key Constants

| Constant | Location | Value | Meaning |
|----------|----------|-------|---------|
| `N` | state.rs | 80 | Grid height (and width for RB) |
| `TRAIL_LEN` | state.rs | 8 | Particle trail ring buffer depth |

Window size (640×320), tiles (3), particles (400), and all physics params are hardcoded defaults in `main.rs` and `SolverParams`.

### Inter-thread Communication

| Channel | Direction | Type | Purpose |
|---------|-----------|------|---------|
| `sync_channel(1)` | physics → main | `FrameSnapshot` | Backpressure-limited frame delivery |
| `mpsc::channel` | main → physics | `SolverParams` | Real-time parameter updates from overlay |
| `mpsc::channel` | main → physics | `(FluidModel, SimState)` | Model switch with fresh state |

## Dependencies

- `icy_sixel` — Pure Rust Sixel encoder (test-only)
- `ctrlc` — Ctrl+C signal handler
- `minifb` — Lightweight native window for pixel buffer display
- `png` — Pure Rust PNG encoder for headless terminal output
- `base64` — Base64 encoding for iTerm2 protocol data payload


## Reference

CFD algorithm ported from [msakuta/cfd-wasm](https://github.com/msakuta/cfd-wasm). See `PRD.md` for full product requirements and `PROGRESS.md` for implementation status.
