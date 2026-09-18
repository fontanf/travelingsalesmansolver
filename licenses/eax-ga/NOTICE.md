# Vendored EAX-TSP

The implementation is adapted from Shujia Liu's C++ implementation of the
Edge Assembly Crossover (EAX) genetic algorithm for the TSP:

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

## Full rewrite into the project's coding style

The 11 files that used to make up this integration (`eax.hpp`/`eax.cpp` plus
the `algorithms/eax/{cross,environment,evaluator,indi,kopt,randomize,sort}.*`
subdirectories) were consolidated into just `include/travelingsalesmansolver/algorithms/eax.hpp`
(all `Distances`-templated code) and `src/algorithms/eax.cpp` (everything
else), and this `LICENSE`/`NOTICE.md` pair moved here, to
`licenses/eax-ga/`, since there is no longer an `algorithms/eax/`
subdirectory for them to live in.

At the same time, every type, method, and member field was renamed to this
project's conventions (PascalCase types with no `T`-prefix, e.g. `TEnvironment`
-> `Environment`, `TIndi` -> `Individual`, `TCross` -> `Cross`, `TKopt` ->
`KOpt`; snake_case methods, e.g. `doIt` -> `run`, `terminationCondition` ->
`termination_condition`; trailing-underscore private members with no
`f`-prefix, e.g. `fEvaluator` -> `evaluator_`), and every raw `new[]`/`delete[]`
array was replaced by `std::vector`/`std::array`, with the two-phase
`define()`-after-construction pattern replaced by ordinary constructors
(RAII). Ultra-local single/two-letter variables inside the crossover and
k-opt inner loops (`r1`/`r2`/`b1`/`b2`, `aa`/`bb`/`a1`/`b1`, `cem`/`ci`/`pr`/`st`,
etc. — these mirror the EAX paper's own red/blue-edge notation) were
deliberately left as-is, to avoid transcription risk in this delicate,
index-heavy code.

`TRandom`/`tRand` and `TSort`/`tSort` — two lazily-allocated global-pointer
singletons wrapping the C `rand()` function and a handful of stateless sort
routines — held no data members of their own, so they were replaced outright
by plain functions (`random_integer`, `random_permutation`, `sort_ascending`,
etc.), removing the process-wide global mutable state and the non-reentrancy
hazard around `InitSort()` that a comment in this file used to document
(`tSort` staying null until the GA reached a "Block2" eset stage, which only
small/quick-converging instances never hit). `InitURandom`/`InitSort`
themselves are gone; the equivalent seeding is `eax_ga::seed_random(seed)`,
called once from `eax()`.

The following genuinely dead code (unreachable from `eax()`, left over from
the original's un-vendored interactive `main.cpp` driver) was removed:
`Environment::printOn`/`writeBest`, `Evaluator::writeTo`/`checkValid`, and
`KOpt::checkDetail`/`checkValid`. A handful of write-only/never-read fields
(`Cross`'s `tBestTmp`, `exam`, `examFlag`, `chDis`, `stock`, a near-duplicate
`fNumOfABcycleInESet` alongside the actually-used `fNumOfABcycleInEset`,
`fJun`, `fOrd1`, `fOrd2`; `Environment`'s `fBestNumOfGen`,
`fBestAccumeratedNumCh`, `fTimeInit`, `fTimeEnd`) were also removed. The
`fInvNearList[i]` fixed-width-500 arrays in `KOpt` (silently truncated per
the original's own comment if a city's inverse near-list ever exceeded 500
entries) became unbounded `std::vector`s.

One pre-existing behavior was deliberately left unchanged despite looking
like a bug: `Cross::set_parents()` (formerly `TCross::setParents()`) computes
a local `int fDisAB = 0` that shadows the class member of the same name, so
the member is never actually updated by `set_parents()` — `Cross::run()`
(formerly `doIt()`)'s later read of that member (`2 * best_number_of_e_edges_
< distance_ab_`) sees a stale/previous-call value rather than the freshly
computed one. This looks like an upstream bug (a missing `this->`), but per
this rewrite's scope, behavior was preserved exactly rather than fixed:
`set_parents()`'s local is named `distance_ab_local` (distinct from the
member `distance_ab_`) precisely so the never-written-back member is
textually obvious rather than hidden behind an identically-named shadow; the
member itself is still never assigned by `set_parents()`, matching the
original's behavior byte-for-byte.

A first pass at this rewrite renamed `Cross`'s member fields and its class
declaration, but left its method *bodies* in the original vendored spacing
(`for( int j =0; j < Num; ++j ){`) with several parameters/locals still
using vendored-style names (`tPa1`, `tPa2`, `tKid`, `Num`, `cem`, `nearMax`).
A follow-up pass reformatted every method body to this project's spacing
convention and renamed those remaining identifiers (`tPa1`/`tPa2`/`tKid` ->
`parent1`/`parent2`/`child`, `Num` -> `number_of_candidates`, `cem` ->
`cycle_length`, `nearMax` -> `near_search_limit`, plus several member fields
that had only gained a trailing underscore without becoming descriptive,
e.g. `ci_` -> `current_city_`, `st_` -> `trace_start_`, `koritsu_many_` ->
`number_of_unbranched_`, `bunki_many_` -> `number_of_branching_`, `flag_st_`
-> `start_new_trace_`, `pr_type_` -> `traversal_type_`). Two more benign
local/member name collisions found during that pass (`form_ab_cycle()`'s and
`make_complete_sol()`'s own locals happening to shadow `trace_start_`/
`current_city_`/`random_pick_` after the first pass's renames, unlike the
`distance_ab_` case above, self-contained and never read back through the
member) were resolved by giving the locals distinct names entirely. Ultra-local
single/two-letter variables that mirror the EAX paper's own red/blue-edge
notation (`r1`/`r2`/`b1`/`b2`, `aa`/`bb`/`a1`/`b1`) were kept as-is throughout,
consistent with `KOpt`.
