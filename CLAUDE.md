# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

A C++14 library and CLI for the traveling salesman problem. It does not implement its own TSP solver from scratch — it provides an `Instance`/`Solution` model plus C++ wrappers that shell out to external solvers: [Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) and [LKH](http://webhotel4.ruc.dk/~keld/research/LKH-3/). A native `local_search` algorithm is also being developed under `src/algorithms/`.

**Licensing**: this repo is MIT-licensed, but Concorde and LKH are not — they are separate executables with their own (non-MIT) licenses, invoked only via subprocess, and never vendored into this repository. Any change that would embed or redistribute their source must respect that.

## Git conventions

- When committing, use this repository's local git identity (`git config --local user.name` / `user.email`), not the machine's global/machine-wide identity. Do not read or rely on global git config for authorship in this repo.
- Never include anything Claude/AI-related in commit messages or PR descriptions: no `Co-Authored-By: Claude ...` line, no "Generated with Claude Code" line, no session IDs, no other AI attribution or signature.

## Build

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cmake --install build --config Release --prefix install
```

CMake options (top-level `CMakeLists.txt`):
- `TRAVELINGSALESMANSOLVER_BUILD_MAIN` (default `ON`) — builds the `travelingsalesmansolver` CLI executable.
- `TRAVELINGSALESMANSOLVER_BUILD_TEST` (default `ON`) — fetches Boost and googletest via `FetchContent` and adds `test/`. Note: `test/CMakeLists.txt` is currently empty — no gtest-based unit tests exist yet; correctness testing is done end-to-end through the CLI (see below), not through this CMake target.

Dependencies (`extern/CMakeLists.txt`) are fetched with `FetchContent`: Boost 1.84 (only needed for `TRAVELINGSALESMANSOLVER_BUILD_TEST`, and also linked into the CLI for `program_options`), googletest, and [fontanf/optimizationtools](https://github.com/fontanf/optimizationtools) (pinned to a specific commit; provides the `Parameters`/`Output`/logging/container base classes this project builds on).

Running the CLI requires the `concorde` and `LKH` executables to be present on `PATH` (not built by this project).

```shell
./install/bin/travelingsalesmansolver -v 1 -i ./data/tsplib/a280.tsp -a concorde
```

Only `-a concorde` and `-a lkh` are actually implemented in `src/main.cpp`'s `run()`; any other value (including the CLI's own default, `large-neighborhood-search`) throws `std::invalid_argument`.

Test data: `python3 scripts/download_data.py` (needs `pip install -r requirements.txt`).

## Testing

There's no unit test suite; correctness is checked by actually running the built CLI with both solvers against a fixed set of small TSPLIB instances and checking the process exit code (`scripts/run_tests.py`). CI (`.github/workflows/build.yml`) runs this on Linux, macOS, and Windows.

`scripts/install_solvers.py` downloads and builds Concorde and LKH for CI/local testing, since neither ships prebuilt for every platform and neither can be vendored into this repo (see licensing note above):
- LKH: built from source (`make`) on Linux/macOS; the official prebuilt `LKH-3.exe` is downloaded directly on Windows.
- Concorde: built from source on Linux/macOS, linked against a prebuilt QSopt LP library fetched for the current OS/architecture (`platform.system()`/`platform.machine()`). **Not installed on Windows** — Concorde's own documentation lists Windows as unsupported without a Cygwin toolchain, so Windows CI only exercises LKH (`scripts/run_tests.py --tests lkh`).

To run the same checks locally: `python3 scripts/install_solvers.py` (installs into `tools/bin`, prepend that to `PATH`), build+install as above, `python3 scripts/download_data.py`, then `python3 scripts/run_tests.py test_results`.

`Instance::read_tsplib` (`instance.cpp`) trims each line before matching TSPLIB field keywords (`optimizationtools::ltrim`/`rtrim`), so field lines with incidental leading/trailing whitespace (e.g. `data/tsplib/ulysses16.tsp`'s ` EOF` line) still parse correctly.

## Architecture

### Library layering (`src/CMakeLists.txt`, `src/algorithms/CMakeLists.txt`)

`TravelingSalesmanSolver_distances` → `TravelingSalesmanSolver_traveling_salesman` (`Instance`, `Solution`, `AlgorithmFormatter`) → per-solver libraries (`TravelingSalesmanSolver_concorde`, `TravelingSalesmanSolver_lkh`) → `TravelingSalesmanSolver_main` (the CLI).

### Distances abstraction and type dispatch

There isn't one `Distances` implementation — there are several (`distances_explicit`, `distances_explicit_triangle`, `distances_euc_2d`, `distances_euc_2d_unrounded`, `distances_ceil_2d`, `distances_geo`, `distances_att`, all under `include/travelingsalesmansolver/distances/`), chosen based on the input instance format. `Distances` (`distances.hpp`) holds an optional `unique_ptr` to each concrete type (only one is populated from parsing; the explicit ones can be lazily computed on demand via `compute_distances_explicit()`).

Algorithm code is written once as a template over the concrete distances type, then dispatched at runtime with the `FUNCTION_WITH_DISTANCES`/`FUNCTION_WITH_DISTANCES_R` macro family in `distances.hpp`, which expands to a ternary chain calling the template instantiation matching whichever concrete member is non-null. This is why every algorithm has both a templated header-only implementation and a small non-template `.cpp` entry point (see below) — the `.cpp` is what actually performs the dispatch via the macro.

### Algorithm module pattern

Each algorithm (`concorde`, `lkh`, `local_search`, under `include/travelingsalesmansolver/algorithms/` + `src/algorithms/`) follows the same shape:
- `XxxParameters : Parameters` and `XxxOutput : Output` structs (with `format()`/`to_json()` overrides).
- A template `xxx(const Distances&, const Instance&, const XxxParameters&)` defined inline in the header — this is the real implementation.
- A non-template `xxx(const Instance&, const XxxParameters&)` declared in the header and defined in the `.cpp`, which dispatches to the template via `FUNCTION_WITH_DISTANCES`.
- Each uses `AlgorithmFormatter` (wraps `optimizationtools::ComposeStream`) to print a consistent progress table and to call `update_solution`/`update_bound`, which in turn trigger `Parameters::new_solution_callback` (used by the CLI to write JSON/certificate output incrementally).

### Concorde / LKH subprocess wrappers

`concorde.hpp` and `lkh.hpp` both: write the instance to a temp file via `Instance::write()` (TSPLIB-style), build a solver-specific command line, run it with `std::system`, then parse the resulting tour file back into a `Solution`. LKH additionally writes/reads a parameters file and a candidate-edges file (`LkhParameters::candidate_file_content`, `read_candidates()`).

Temp files are created through the shared `make_temp_path()` helper in `include/travelingsalesmansolver/algorithms/temp_file.hpp` rather than `tmpnam` — it prefers a tmpfs (`/dev/shm`) on Linux to avoid real disk I/O, falls back to `$TMPDIR`/`/tmp` elsewhere (e.g. macOS), and uses the Win32 `GetTempPath`/`GetTempFileName` API on Windows. Any new external-solver wrapper following this pattern should reuse this helper instead of reintroducing `tmpnam`.
