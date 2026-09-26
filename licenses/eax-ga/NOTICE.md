# Vendored EAX-TSP

The implementation is adapted from Shujia Liu's C++ implementation of the
Edge Assembly Crossover (EAX) genetic algorithm for the TSP:

https://github.com/Sugia/GA-for-TSP

Licensed under the Apache License, Version 2.0; see the `LICENSE` file in
this directory for the full text.

Changes made to the original source:
- Renamed `.h` headers to `.hpp` and namespaced everything under an
  anonymous namespace nested inside `travelingsalesmansolver` (originally
  a named `travelingsalesmansolver::eax_ga`; see further down for why this
  changed).
- `TEvaluator`, `TCross`, `TKopt`, and `TEnvironment` are now templated on
  a `Distances` type and hold a `const Distances&` (via `TEvaluator`)
  instead of copying distances into an `int` matrix (`fEdgeDis`) or reading
  them from a TSPLIB file: `TEvaluator::distance(i, j)` looks them up
  directly on the original object, dispatched once per `local_search()` call via
  the `Distances` template parameter, the same way every other algorithm
  in this library consumes distances. As part of this, every
  distance/gain/tour-length-valued field or local variable that used to be
  a plain `int` (`TIndi::fEvaluationValue`, `TEvaluator::fEdgeDis`,
  `TCross::fGainAB`/`fGainModi` and the `diff`/`gain`/`BestGain` locals
  that flow from them, `TKopt::fTourLength` and its `dis1`/`dis2` locals,
  `TEnvironment::fBestValue`) is now the shared `Distance` (`int64_t`)
  type, so distances scaled well beyond a 32-bit `int`'s range (e.g. the
  minmax mTSP paper's unrounded Euclidean convention, `EUC_2D` scaled by `1e6`) work
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
  they are replaced by `travelingsalesmansolver::local_search()` (originally
  named `eax()`; see further down for why it was renamed) in
  `include/travelingsalesmansolver/algorithms/local_search.hpp`.

## Full rewrite into the project's coding style

The 11 files that used to make up this integration (`eax.hpp`/`eax.cpp` plus
the `algorithms/eax/{cross,environment,evaluator,indi,kopt,randomize,sort}.*`
subdirectories) were consolidated into just `include/travelingsalesmansolver/algorithms/local_search.hpp`
(all `Distances`-templated code) and `src/algorithms/local_search.cpp` (everything
else), and this `LICENSE`/`NOTICE.md` pair moved here, to
`licenses/eax-ga/`, since there is no longer an `algorithms/eax/`
subdirectory for them to live in.

At the same time, every type, method, and member field was renamed to this
project's conventions (PascalCase types with no `T`-prefix, e.g. `TEnvironment`
-> `Environment`, `TIndi` -> `Individual` (later renamed `LocalSearchSolution`), `TCross` -> `Cross`, `TKopt` ->
`KOpt`; snake_case methods, e.g. `doIt` -> `run`, `terminationCondition` ->
`termination_condition`; trailing-underscore private members with no
`f`-prefix, e.g. `fEvaluator` -> `evaluator_`), and every raw `new[]`/`delete[]`
array was replaced by `std::vector`/`std::array`, with the two-phase
`define()`-after-construction pattern replaced by ordinary constructors
(RAII).

`TRandom`/`tRand` and `TSort`/`tSort` — two lazily-allocated global-pointer
singletons wrapping the C `rand()` function and a handful of stateless sort
routines — held no data members of their own, so they were replaced outright
by the standard library, removing the process-wide global mutable state and the
non-reentrancy hazard around `InitSort()` that a comment in this file used to
document (`tSort` staying null until the GA reached a "Block2" eset stage,
which only small/quick-converging instances never hit). Random numbers are now
drawn from a `std::mt19937_64` owned by `LocalSearchData` and seeded with
`LocalSearchParameters::seed` (`std::uniform_int_distribution`,
`std::shuffle`), and sorting uses `std::sort`/`std::stable_sort`.
`InitURandom`/`InitSort` themselves are gone. As a consequence, the sequence
of random numbers, and hence the tours found for a given seed, differ from the
original's.

The following genuinely dead code (unreachable from `local_search()`, left over from
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

`Cross::set_parents()` (formerly `TCross::setParents()`) computed a local
`int fDisAB = 0` that shadowed the class member of the same name, so the
member was never updated and `Cross::run()`'s later read of it
(`2 * best_number_of_e_edges_ < distance_ab_`) used an uninitialized value.
The member, the local and that condition were removed; the check on
`child.length != parent2.length` remains.

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
`current_city_`/`random_pick_` after the first pass's renames, self-contained and never read back through the
member) were resolved by giving the locals distinct names entirely. That
pass deliberately kept the EAX paper's own red/blue-edge notation
(`r1`/`r2`/`b1`/`b2`, `aa`/`bb`/`a1`/`b1`) as-is.

A third pass superseded that last call: every vertex/city-id-valued
variable, parameter, and member across all five classes (`Individual`,
`Evaluator`, `KOpt`, `Cross`, `Environment`) and `eax.cpp` was retyped from
`int` to this project's `VertexId` (`int64_t`) — the type every
`Distances::distance(...)` overload already took, so these were previously
narrowed on every call — and every single/two-letter vertex-id variable was
renamed to `vertex_id`/`vertex_id_1`/`vertex_id_2`/... or, where the EAX
paper's red/blue-edge roles would otherwise be lost, to role-based names
(`red_vertex_id_1`/`red_vertex_id_2`/`blue_vertex_id_1`/`blue_vertex_id_2`
for `r1`/`r2`/`b1`/`b2`; `vertex_id_1..4` for `aa`/`bb`/`a1`/`b1`, matching
how they're stored in `modified_edge_[s][0..3]`). Position/index/count
variables that are *not* vertex ids (segment ids, AB-cycle indices, tour
positions such as `order_`'s indices, population/generation counters) were
deliberately left as plain `int`, including the mixed-semantics `ab_cycle_`
array (`ab_cycle_[c][0]` is a cycle-length count, `ab_cycle_[c][1..]` are
vertex ids) — see the source for the remaining per-variable judgment calls.
The nested `eax_ga` namespace also became anonymous in this pass, since
every template class in it is only ever instantiated from `local_search.cpp`'s
translation unit (via `FUNCTION_WITH_DISTANCES`) — `main.cpp` only calls
the non-template `local_search()` — so there is no ODR risk despite this
being a header.

## Renamed the public API from `eax` to `local_search`

`eax.hpp`/`eax.cpp` were renamed to `local_search.hpp`/`local_search.cpp`,
`EaxParameters` to `LocalSearchParameters`, the `eax()` functions to
`local_search()`, the CLI algorithm name from `-a eax` to `-a local_search`,
and the printed algorithm name from `EAX` to `Local Search`. This only
changes this project's own label for the entry point; the algorithm itself
is still literally the Edge Assembly Crossover genetic algorithm described
above, and this directory's name (`licenses/eax-ga/`) and the vendored
attribution in this file are unchanged since they describe the actual
upstream technique, not this project's chosen name for it.

## Flattened `Evaluator`/`KOpt`/`Cross`/`Environment` into `LocalSearchData` + free functions

Following the same pattern as `local_search_pfss_makespan.cpp` (a single
`struct LocalSearchData` holding all mutable working state, with every
operation a free function taking it in place of an implicit `this`), the
four classes were merged into one `template <typename Distances> struct
LocalSearchData` plus free functions in the same anonymous namespace.
`Individual` (now `LocalSearchSolution`) is unaffected (it was already a plain data holder with no
methods beyond `operator==`, which is unchanged). Unlike
`local_search_pfss_makespan`'s helpers, these free functions do not take a
`const Instance&`: they never needed one (distances/near-neighbor lists
were already cached in what was `Evaluator`), so it would have been an
unused parameter forced in purely to match the reference pattern's shape.

Flattening surfaced several same-named members from different classes that
would have silently aliased into one field (each class currently keeps its
own separate copy, so no behavior changed before this point, only after
flattening would collapsing them have become a real bug). These were
deliberately kept distinct (see `LocalSearchData`'s own comment in
`local_search.hpp` for the equivalent detail):
- `KOpt::number_of_segments_` (current number of segments in the 2-opt tree
  representation) vs `Cross::number_of_segments_` (path segments built by
  `make_complete_sol()`/`make_unit()`) -> `number_of_tree_segments` and
  `number_of_segments` respectively.
- `Cross::max_stagnation_` ("Block2" eset local search stagnation limit in
  `search_eset()`) vs `Environment::max_stagnation_` (generation-level
  stage-1 -> stage-2 stagnation limit) -> `eset_max_stagnation` and
  `max_stagnation` respectively. This is the one the user flagged up front
  as the most dangerous of the three: since "Block2" is only reached after
  enough stalled generations, a small/quick-converging instance (e.g.
  a280) would never exercise the colliding field and the corruption would
  only show up as quietly worse solutions on harder instances — the same
  general hazard class as the `tSort`/"Block2" non-reentrancy issue
  documented earlier in this file.
- `KOpt::next_city()`/`previous_city()` (`next_vertex()`/`previous_vertex()` after flattening)
  vs `Cross::current_city_`/`previous_city_` (AB-cycle trace state) — not
  found by inspection, only by mechanically stripping every member's
  trailing underscore and diffing all four classes' method and field names
  against each other pairwise. `Cross`'s fields became
  `trace_current_vertex`/`trace_previous_vertex` (with `trace_start_` also
  renamed to `trace_start` for consistency, though it had no collision).

Fields that were genuinely the same value in multiple classes (all set from
the same constructor argument, so unifying them is not a behavior change)
were merged into one field: `number_of_vertices` (was
`Evaluator::number_of_vertices`/`KOpt::number_of_vertices_`/
`Cross::number_of_vertices_`), `population_size` (was
`Cross::population_size_`/`Environment::population_size_`), and the
nearest-neighbor list size (was `Evaluator::max_near_cities_`/
`KOpt::max_near_cities_used_`, both hardcoded to `50`) into one file-scope
`constexpr int max_near_vertices`. `Evaluator` itself disappears as a
sub-object: its fields (`distances`, `near_vertices`, `number_of_vertices`)
are now direct fields of `LocalSearchData`, and call sites that used to say
`evaluator_.distance(a, b)` now say `distance(data, a, b)` via a small
inline free-function wrapper around `data.distances.distance(a, b)`.

Two now-redundant parameters were dropped as a direct, mechanical
consequence of the flattening, not a separate design change: `Cross`'s
`flags`/`edge_frequency` parameters (originally threaded in by reference
from `Environment`'s `fFlagC`/`fEdgeFreq` fields, since `Cross` had no
access to `Environment`'s data otherwise) are gone from `run_cross()`,
`set_parents()`, `set_ab_cycle()`, `increment_edge_freq()`,
`calc_adaptive_loss()`, and `calc_entropy_loss()`, which now read
`data.flags`/`data.edge_frequency` directly like every other field.

`LocalSearchData`'s constructor is kept minimal, matching the reference
pattern's spirit: it only binds the `const Distances&` reference member
(which cannot be default-constructed and assigned later) and sizes/computes
everything that depends solely on `number_of_vertices` (near-neighbor
lists, the 2-opt tree's working arrays, the inverse near-neighbor list —
this is genuine computation, not just sizing, so it stays in the
constructor rather than moving to `init()`). Everything that additionally
depends on `population_size`/`number_of_children` (formerly split across
`Cross`'s and `Environment`'s constructors) moved to a free function
`init(data, population_size, number_of_children)`, called once from
`local_search()`'s entry point — mirroring how `local_search_pfss_makespan()`
sizes `data.completion_times_0` etc. directly in its own entry point rather
than in a constructor.

Three methods were unambiguously renamed since, as free functions in one
shared namespace, they could no longer be disambiguated by class scope:
`Environment::run()` (the top-level generation loop, called from the
public `local_search()` entry point) stays plain `run()`; `KOpt::run()`
(local search on one solution) becomes `run_kopt()`; `Cross::run()`
(crossover for one mating pair) becomes `run_cross()`.
