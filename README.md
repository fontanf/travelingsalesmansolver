# Traveling salesman solver

A solver for the traveling salesman problem.

![travelingsalesman](img/travelingsalesman.png?raw=true)

[image source](https://commons.wikimedia.org/wiki/File:GLPK_solution_of_a_travelling_salesman_problem.svg)

This package provides wrappers to call [Concord TSP Solver](https://www.math.uwaterloo.ca/tsp/concorde.html) and [LKH](http://webhotel4.ruc.dk/~keld/research/LKH-3/) from C++.

**Beware that even though this package is licensed under the terms of the MIT license, neither are the Concord TSP Solver nor the LKH packages.**

Both solvers are called through system calls. Inputs and outputs are handled through files. It is required to have the `concorde` and the `LKH` executables in the path.

## Usage (command line)

Compile:
```shell
bazel build -- //...
```

Examples:

```shell
./bazel-bin/travelingsalesmansolver/main -v 1 -i ./data/tsplib/a280.tsp -a concorde
```
```
=====================================
      Traveling salesman solver      
=====================================

Instance
--------
Number of vertices:  280

Algorithm
---------
Concorde

       T (s)          UB          LB         GAP     GAP (%)                 Comment
       -----          --          --         ---     -------                 -------
       0.000         inf           0         inf         inf                        
       0.614        2579           0        2579         inf          final solution
       0.614        2579        2579           0        0.00          final solution

Final statistics
----------------
Value:                         2579
Bound:                         2579
Absolute optimality gap:       0
Relative optimality gap (%):   0
Time (s):                      0.614178

Solution
--------
Number of vertices:  280 / 280 (100%)
Feasible:            1
Distance:            2579
```
