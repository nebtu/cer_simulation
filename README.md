# Simulation of adaptive multi-arm trials with `adagraph`

This repo provides the code used for the simulation of some adaptive multi-arm
trials using a graph-based multiple-testing strategy with the
[adagraph](https://github.com/nebtu/adagraph) package.
The scenarios are taken from Mehta \[1\], and extended to also simulate mixed
continuous and binary endpoints.
An overview over the tested scenarios can be found in [results.html], as well as some results as produced by the current code.

## Usage

For simply reproducing the whole simulation setup, clone the repository, install all dependencies (as
listed in the [_targets.R] file) and run

```r
targets::tar_make_future(workers = 40)
```
while setting `workers` to the appropriate amount of parallel workers.
On 40 cores, the whole setup takes about 20 hours.

The exact scenarios being run are defined in [_targets.R](_targets.R) as the `power_analysis` and
`fwer_analysis` dataframes.
Inspect them after running `targets::tar_load_globals()`
To run a simulation with different or only a subset of the scenarios, change the
definition of those dataframes to include the same column names, with each row
defining one scenario.

## Structure

The simulation uses the [targets](https://books.ropensci.org/targets/) package.
It specifically makes use of [dynamic](https://books.ropensci.org/targets/dynamic.html) and
[static branching](https://books.ropensci.org/targets/static.html) for defining the different scenarios.

In the [_targets.R] file, the general definitions specifiying the workflow are
done. First, we define the `power_analyis` and `fwer_analysis` dataframes, which
specify the parameters used for the simulations.

Using the `tar_map()` function, we define now for each of those scenarios some
targets. These are used for generating the data generating function, the
adaption function, conducting the actual simulation and (split into two steps)
aggregating the results.
All defined targets can be viewed by calling `tar_visnetwork()` or as a
dataframe using `tar_manifest()`.

In the [sim_setup.R](R/sim_setup.R) file, we define the functions for running the
simulations. These are mostly calls to the appropriate `adagraph` functions.
The function `get_sim_adaption()` definies the exact adaption rules used.

In [process_sim.R](R/process_sim.R), functions for aggregating simulation results,
calculating the operation characteristics and displaying them are defined.

## References

\[1\] C. Mehta, A. Mukhopadhyay, and M. Posch, “Graph Based, Adaptive, Multi Arm, Multiple Endpoint, Two Stage Design,” Jan. 07, 2025, arXiv: arXiv:2501.03197. doi: 10.48550/arXiv.2501.03197.
