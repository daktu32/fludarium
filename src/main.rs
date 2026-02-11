mod renderer;
mod sixel;
mod solver;
mod state;

use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

fn main() {
    // Enter alternate screen and hide cursor
    let stdout = io::stdout();
    {
        let mut handle = stdout.lock();
        let _ = handle.write_all(b"\x1b[?1049h"); // alternate screen
        let _ = handle.write_all(b"\x1b[?25l"); // hide cursor
        let _ = handle.write_all(b"\x1b[?80h"); // disable Sixel scrolling
        let _ = handle.write_all(b"\x1b[2J"); // clear screen
        let _ = handle.flush();
    }

    // Set up Ctrl+C handler
    let running = Arc::new(AtomicBool::new(true));
    let r = running.clone();
    ctrlc::set_handler(move || {
        r.store(false, Ordering::SeqCst);
    })
    .expect("Error setting Ctrl+C handler");

    // Initialize simulation
    let mut sim = state::SimState::new();
    let params = solver::SolverParams::default();
    // Main loop with frame pacing
    let steps_per_frame = 5;
    let frame_duration = Duration::from_millis(50); // ~20 FPS cap
    while running.load(Ordering::SeqCst) {
        let frame_start = Instant::now();

        for _ in 0..steps_per_frame {
            solver::fluid_step(&mut sim, &params);
        }

        let rgba = renderer::render(&sim);

        match sixel::encode_sixel(&rgba, renderer::FRAME_WIDTH, renderer::FRAME_HEIGHT) {
            Ok(sixel_data) => {
                if sixel::output_frame(&sixel_data).is_err() {
                    break;
                }
            }
            Err(e) => {
                eprintln!("Sixel encoding error: {}", e);
                break;
            }
        }

        let elapsed = frame_start.elapsed();
        if elapsed < frame_duration {
            std::thread::sleep(frame_duration - elapsed);
        }
    }

    // Restore terminal
    {
        let mut handle = stdout.lock();
        let _ = handle.write_all(b"\x1b[?80l"); // re-enable Sixel scrolling
        let _ = handle.write_all(b"\x1b[?25h"); // show cursor
        let _ = handle.write_all(b"\x1b[?1049l"); // exit alternate screen
        let _ = handle.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pipeline_no_panic() {
        // Integration test: run 3 frames of the full pipeline
        let mut sim = state::SimState::new();
        let params = solver::SolverParams::default();

        for _ in 0..3 {
            solver::fluid_step(&mut sim, &params);
            let rgba = renderer::render(&sim);
            let result = sixel::encode_sixel(&rgba, renderer::FRAME_WIDTH, renderer::FRAME_HEIGHT);
            assert!(result.is_ok(), "Sixel encoding should succeed");
            let data = result.unwrap();
            assert!(!data.is_empty(), "Sixel output should not be empty");
        }
    }
}
