# Mandelbrot-Explorer

Mandelbrot deep-zoom renderer using perturbation, series approximation, reference
rebasing, BLA, exterior distance estimation, and a floatexp rescaled path for
correct rendering past double's ~1e320 underflow (validated to 1e1000+).

## Requirements

- MSVC (x64) with OpenMP and AVX2
- GMP via vcpkg: `vcpkg install gmp:x64-windows-static`

No FreeGLUT / GLEW / OpenGL / libpng / png++ are required: the interactive UI is
native Win32 (GDI) and PNG export uses the Windows Imaging Component. The build
scripts auto-detect Visual Studio, the Windows SDK, and vcpkg; set the
`VCPKG_ROOT` environment variable if your vcpkg lives somewhere unusual.

## Project layout

    src/engine/  core renderer (perturbation, floatexp, BLA, coloring, GMP)
    src/gui/     native Win32 explorer + navigator
    src/tools/   headless render / verify / diff utilities
    src/bench/   micro-benchmarks and demos
    scripts/     build_*.bat + msvcenv.bat (MSVC / SDK / vcpkg auto-detection)
    build/       compiler output

## Build

Build scripts live in `scripts\` and emit into `build\` (each `cd`s to the repo
root itself, so they can be run from anywhere).

- `scripts\build_gui.bat` — interactive native Win32 explorer -> `build\mandel_gui.exe`
  (single exe; `vcomp140.dll`, the OpenMP runtime, is copied alongside so the
  `build\` folder is copy-and-run portable).
- `scripts\build_verify.bat` — headless correctness/benchmark harness -> `build\verify.exe`
- `scripts\build_render.bat` — headless BMP renderer -> `build\render.exe`
- `scripts\build_pbench.bat` — clean-room perturbation/BLA bench -> `build\pbench.exe`
- `scripts\build_bench.bat` — BigFixed vs GMP micro-bench -> `build\bench_bigfixed.exe`

## Run

### Interactive explorer

    build\mandel_gui.exe

Controls:

- Drag: pan
- Mouse wheel / double-click: zoom in; right-click: zoom out
- `R`: reset, `Space`: re-render, `S`: save PNG, `C`: copy location, `+`/`-`: zoom
- Side panel: max iterations, color density, 3x supersampling, exterior distance
  estimation, palette preset, per-stop color (double-click a gradient marker;
  left-drag to move, right-click to delete, click the bar to add).
- `Export`: high-resolution PNG (monitor-aware resolution presets or custom size,
  live preview, progress, cancel; coordinate range matches the view).

### Headless render

    build\render.exe out.bmp W H cx cy scaleExp mxit [SS]
    build\render.exe                      (no args: renders a 1e1000 exterior demo)

### Verify / benchmark

    build\verify.exe [shallow|deep|ticktock|flake|exterior1000|parity1000|all] [W] [H]
