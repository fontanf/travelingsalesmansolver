# Vendored EAX-TSP

The files in this directory (`cross.*`, `environment.*`, `evaluator.*`, `indi.*`,
`kopt.*`, `randomize.*`, `sort.*`) are adapted from Shujia Liu's C++
implementation of the Edge Assembly Crossover (EAX) genetic algorithm for the
TSP:

https://github.com/Sugia/GA-for-TSP

Licensed under the Apache License, Version 2.0; see the `LICENSE` file in
this directory for the full text.

Changes made to the original source:
- Renamed `.h` headers to `.hpp` and namespaced everything under
  `travelingsalesmansolver::eax_ga`.
- `TEvaluator`, `TCross`, `TKopt`, and `TEnvironment` are now templated on
  a `Distances` type and hold a `const Distances&` (via `TEvaluator`)
  instead of copying distances into an `int` matrix (`fEdgeDis`) or reading
  them from a TSPLIB file: `TEvaluator::distance(i, j)` looks them up
  directly on the original object, dispatched once per `eax()` call via
  the `Distances` template parameter, the same way every other algorithm
  in this library consumes distances. As part of this, every
  distance/gain/tour-length-valued field or local variable that used to be
  a plain `int` (`TIndi::fEvaluationValue`, `TEvaluator::fEdgeDis`,
  `TCross::fGainAB`/`fGainModi` and the `diff`/`gain`/`BestGain` locals
  that flow from them, `TKopt::fTourLength` and its `dis1`/`dis2` locals,
  `TEnvironment::fBestValue`) is now the shared `Distance` (`int64_t`)
  type, so distances scaled well beyond a 32-bit `int`'s range (e.g. the
  minmax mTSP paper's unrounded-`EUC_2D` convention, scaled by `1e6`) work
  correctly without the caller needing to know a scale exists at all.
  `TEvaluator`'s original file-based `setInstance(filename)` (reading a
  TSPLIB file and computing its own `EUC_2D`/`ATT`/`CEIL_2D` distances from
  coordinates) could not coexist with holding an external, generic
  `Distances` reference, and was already dead code for this integration
  (only the in-memory path was ever used), so it was removed along with
  the `x`/`y` coordinate storage it alone needed.
- `TEvaluator::getTour(TIndi&)` returns the 0-indexed tour without going
  through a temporary file.
- Added `fTimeLimit`, a wall-clock cutoff checked in
  `TEnvironment::terminationCondition()`.
- Removed debug `printf`/`std::cout` statements in `TEnvironment::doIt()`
  and `TEnvironment::getEdgeFreq()`.
- Since every class above is now a template, its member function bodies
  moved from the corresponding `.cpp` into the `.hpp` (matching how every
  other template in this library is structured) â€” `cross.cpp`,
  `environment.cpp`, `evaluator.cpp`, and `kopt.cpp` no longer exist.
  `indi.*`, `randomize.*`, and `sort.*` are untouched by this and remain
  separately-compiled, since none of them depend on `Distances`.
- `main.cpp`/`main.h` (the original's interactive driver) were not vendored;
  they are replaced by `travelingsalesmansolver::eax()` in
  `include/travelingsalesmansolver/algorithms/eax.hpp`.
