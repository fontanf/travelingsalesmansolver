# Vendored EAX-TSP

The files in this directory (`cross.*`, `environment.*`, `evaluator.*`, `indi.*`,
`kopt.*`, `randomize.*`, `sort.*`) are adapted from Shujia Liu's C++
implementation of the Edge Assembly Crossover (EAX) genetic algorithm for the
TSP:

https://github.com/Sugia/GA-for-TSP

Licensed under the Apache License, Version 2.0. A copy of the license is
available at:

http://www.apache.org/licenses/LICENSE-2.0

Changes made to the original source:
- Renamed `.h` headers to `.hpp` and namespaced everything under
  `travelingsalesmansolver::eax_ga`.
- Added `TEvaluator::setInstance(int, const vector<vector<int>>&)` and
  `TEvaluator::getTour(TIndi&)` to support building an instance from an
  in-memory distance matrix and extracting a tour without going through
  temporary files.
- Added `TEnvironment::defineFromEvaluator()` and a `fTimeLimit` wall-clock
  cutoff checked in `TEnvironment::terminationCondition()`.
- Removed debug `printf`/`std::cout` statements in `TEnvironment::doIt()`,
  `TEnvironment::getEdgeFreq()`, and `TEvaluator::setInstance()`.
- `main.cpp`/`main.h` (the original's interactive driver) were not vendored;
  they are replaced by `travelingsalesmansolver::eax()` in
  `include/travelingsalesmansolver/algorithms/eax.hpp`.
