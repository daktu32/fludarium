# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fluvarium is a terminal-based fluid dynamics visualizer. It runs a Rayleigh-Benard convection simulation (hot bottom, cold top) and renders it as Sixel graphics in WezTerm. Written in Rust (edition 2024), single crate, no external C dependencies.

## Build & Test Commands

```bash
cargo build              # Debug build
cargo build --release    # Release build (recommended for running)
cargo test               # Run all 32 tests
cargo test solver        # Run only solver module tests
cargo test state         # Run only state module tests
cargo test renderer      # Run only renderer module tests
cargo test sixel         # Run only sixel module tests
cargo test test_name     # Run a single test by name
cargo run --release      # Run the simulation (requires WezTerm with Sixel support)
```

## Architecture

Single-threaded pipeline running in a loop:

```
SimState → fluid_step() → render() → encode_sixel() → output_frame() → stdout
```

### Module Responsibilities

- **`state.rs`** — `SimState` struct (velocity fields `vx/vy/vx0/vy0`, `temperature`, `density`, scratch buffers `work/work2`, `Xor128` PRNG). Grid is `N=256`, flat `Vec<f64>` of `SIZE = N*N`. The `idx(x, y)` function handles 2D→1D with mod-N wrapping for horizontal periodicity.
- **`solver.rs`** — Jos Stam "Stable Fluids" CFD: `diffuse`, `advect` (Semi-Lagrangian), `project` (pressure projection), `fluid_step` (orchestrates one timestep). `SolverParams` holds all tuning constants. Boundary conditions via `set_bnd(field_type, ...)` where field_type 0=scalar, 1=vx, 2=vy(negate at walls), 3=temperature(Dirichlet).
- **`renderer.rs`** — Maps temperature field to RGBA via 5-stop colormap (Blue→Cyan→Green→Yellow→Red). Renders at 2x scale (`DISPLAY_SIZE=512`) with a color bar + tick marks. Output dimensions: `FRAME_WIDTH` x `FRAME_HEIGHT`.
- **`sixel.rs`** — Wraps `icy_sixel` for encoding. `output_frame` uses synchronized output (DEC 2026 mode) to prevent flicker.
- **`main.rs`** — Alternate screen setup, `ctrlc` handler, main loop with 5 sim steps per frame at ~20 FPS cap, terminal restore on exit.

### Coordinate System

- `y=0` is **bottom** (hot), `y=N-1` is **top** (cold) in simulation space
- Renderer flips Y for screen output (screen top = cold/blue, screen bottom = hot/red)
- X-axis wraps (periodic boundary), Y-axis has wall boundaries (top/bottom)

### Key Constants

| Constant | Location | Value | Meaning |
|----------|----------|-------|---------|
| `N` | state.rs | 256 | Grid dimension |
| `SIZE` | state.rs | 65536 | N*N total cells |
| `SCALE` | renderer.rs | 2 | Display upscale factor |
| `DISPLAY_SIZE` | renderer.rs | 512 | N*SCALE pixels |

## Dependencies

- `icy_sixel` — Pure Rust Sixel encoder
- `image` — Image buffer types (declared but minimally used; rendering is manual)
- `ctrlc` — Ctrl+C signal handler

## Reference

CFD algorithm ported from [msakuta/cfd-wasm](https://github.com/msakuta/cfd-wasm). See `PRD.md` for full product requirements and `PROGRESS.md` for implementation status.
