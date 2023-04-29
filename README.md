# Traveling salesman solver

A solver for the traveling salesman problem.

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
