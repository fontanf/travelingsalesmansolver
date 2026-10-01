# TravelingSalesmanSolver

A solver for the traveling salesman problem.

![travelingsalesman](img/travelingsalesman.png?raw=true)

[image source](https://commons.wikimedia.org/wiki/File:GLPK_solution_of_a_travelling_salesman_problem.svg)

## Implemented algorithms

- Construction algorithms:
  - Random permutation `--algorithm random-permutation`
  - Random walk on the candidate graph (as LKH's `WALK` initial tours) `--algorithm random-walk`
  - Greedy edge (Bentley's multiple fragment heuristic) `--algorithm greedy-edge`

- 2-opt local search, from a random permutation or an initial solution `--algorithm two-opt [--initial-solution solution.txt]`

- LKH-style Lin-Kernighan `--algorithm lin-kernighan`: sequential moves of up to 5 edges, non-sequential moves, trials merged by iterative partial transcription; written from the published descriptions (no LKH code). Its core is a problem-generic engine (`include/travelingsalesmansolver/lin_kernighan/engine.hpp`), meant to be reused by Lin-Kernighan algorithms for other problems: a problem class gives the costs of the edges guiding the search and, if needed, evaluates the moves exactly (see the documentation at the top of the file). Optionally, the search is guided by the costs penalized by the ascent of the alpha-nearness candidates, as LKH `--penalized-costs 1` (it didn't improve the results on TSPLIB instances).

- Genetic algorithm with edge assembly crossover (EAX), written from Nagata and Kobayashi, "A Powerful Genetic Algorithm Using Edge Assembly Crossover for the Traveling Salesman Problem" (INFORMS Journal on Computing, 2013): localized EAX (single AB-cycles) then EAX with the block2 strategy, entropy-preserving selection `--algorithm eax`. The initial tours of the population are random permutations or random walks `--initial-tours random-permutation|random-walk`, improved by a 2-opt (as in the paper), the first local search of the Lin-Kernighan, or nothing `--initial-local-search 2-opt|lin-kernighan|none`

- Wrappers calling external solvers through system calls (inputs and outputs through files; the `LKH` and `concorde` executables must be in the path):
  - [LKH](http://webhotel4.ruc.dk/~keld/research/LKH-3/) `--algorithm lkh`
  - [Concorde TSP Solver](https://www.math.uwaterloo.ca/tsp/concorde.html) `--algorithm concorde`

**Beware that even though this package is licensed under the terms of the MIT license, neither are the Concorde TSP Solver nor the LKH packages.**

Candidate edges (`include/travelingsalesmansolver/candidates/`), given to the algorithms as an optional argument `--candidates nearest-neighbor|alpha-nearness --number-of-candidates 5`:
* `nearest_neighbor_candidates`: the nearest neighbors of each vertex (default of `greedy_edge` and `two_opt`)
* `alpha_nearness_candidates`: alpha-nearness (Helsgaun, 2000), computed from 1-trees with vertex penalties optimized by a subgradient ascent, on a sparse graph (nearest neighbors, and nearest neighbors in each octant for coordinate-based distances) or on the complete graph (default of `random_walk` and `lin_kernighan`: 5 candidates; with the Lin-Kernighan, average gap to the optimum 0.07% instead of 0.31% with 10 nearest neighbors on rl5915, pla7397, rl11849 and usa13509, 60 s, 3 seeds; 0.01% instead of 0.22% on 7 instances of 783 to 4461 vertices, 30 s, 5 seeds)

## Usage (command line)

Compile:
```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cmake --install build --config Release --prefix install
```

Download data:
```shell
python3 scripts/download_data.py
```

Solve:
```shell
./install/bin/travelingsalesmansolver --verbosity-level 1 --input ./data/tsplib/a280.tsp --format tsplib --algorithm lin-kernighan --time-limit 1 --output a280_output.json --certificate a280_solution.txt
```
```
=====================================
       TravelingSalesmanSolver       
=====================================

Instance
--------
Number of vertices:  280

Algorithm
---------
Lin-Kernighan

Parameters
----------
Time limit:                          1
Messages
    Verbosity level:                 1
    Standard output:                 1
    File path:                       
    # streams:                       0
Logger
    Has logger:                      0
    Standard error:                  0
    File path:                       
Number of candidates:                5
Candidates:                          alpha-nearness
Penalized costs:                     0
Ascent initial period:               -1
Move type:                           5
Maximum depth:                       50
Perturbation:                        walks
Restricted search:                   1
Non-sequential moves:                1
Maximum number of trials:            -1
Kick segment length:                 50

    Time (s)       Value       Bound         Gap     Gap (%)                 Comment
    --------       -----       -----         ---     -------                 -------
       0.023         inf           0         inf         inf                        
       0.023        3093           0        3093         inf            initial tour
       0.041        2579           0        2579         inf           Lin-Kernighan

Final statistics
----------------
Value:                        2579
Bound:                        0
Absolute optimality gap:      2579
Relative optimality gap (%):  inf
Time (s):                     1.00489

Solution
--------
Number of vertices:  280 / 280 (100%)
Feasible:            1
Distance:            2579
```
