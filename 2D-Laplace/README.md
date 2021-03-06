# Prerequisite
1. Any MPI library (tested with OpenMPI 3.1.3 and MPICH 3.3)
2. Gnuplot

# Build
Run command:

``` make ```

This should make everything, including the following executibles:
Init, InitPar, Combine, Solve, SolveParRow, SolvePar, ReadGridParamsTest


# Clean up
Run command:

``` make clean ```

# Run
## Initialise boundary values
__Init__ and __InitPar__ executibles are responsible for initialising boundary values.
Both require a __parameter file__ in the same format as "parameters.in".

### __Init__ usage:

``` ./Init <ParamFile> <FunctionSelection> ```

Where \<ParamFile\> is the parameter file, and \<FunctionSelection\> is an integer from 0 to 3 (inclusive), which selects a
function we'd like to solve.

__Init__ outputs one single "initial.dat" file and also outputs a "solution.dat" file with the analytical solution

### __InitPar__ usage:

``` ./InitPar <NumPatchInX> <NumPatchInY> <ParamFile> <FunctionSelection> ```

Where \<NumPatchInX\> and \<NumPatchInY\> are integers (greater than 0), specifying how the global grid is partitioned into
local patches in both X and Y axes.
\<ParamFile\> is the parameter file, and \<FunctionSelection\> is an integer from 0 to 3 (inclusive), which selects a
function we'd like to solve.
  
__InitPar__ outputs \<NumPatchInX\> * \<NumPatchInY\> number of "initial.MPI_\<Rank\>.dat" and "solution.MPI_\<Rank\>.dat" files.

**Note that InitPar must be run with \<NumPatchInX\> * \<NumPatchInY\> number of processes.**

## Solve 2D laplace
__Solve__, __SolveParRow__ and __SolvePar__ solve the 2D laplace equations.
__Solve__ is a serial implementation; __SolveParRow__ parallelises across the rows and __SolvePar__ parallelises across both rows
and columns.

### __Solve__ usage:

``` ./Solve <ParamFile> ```

Where \<ParamFile\> is the parameter file.

__Solve__ then takes "initial.dat" as input file and outputs "laplace.dat".

### __SolveParRow__ usage:

``` ./SolveParRow <ParamFile> ```

Where \<ParamFile\> is the parameter file.
__SolveParRow__ then takes "initial.MPI_\<Rank\>.dat" as input files and outputs "laplace.MPI_\<Rank\>.dat".

**Note that SolveParRow must be run with the same number of processes as the number of local patches across the rows.
This also means that the number of processes must be the same as the number of "initial.MPI_\<Rank\>.dat" files.**

### __SolvePar__ usage:

``` ./SolvePar <NumPatchInX> <NumPatchInY> <ParamFile> ```

Where \<NumPatchInX\> and \<NumPatchInY\> are integers (greater than 0), specifying how the global grid is partitioned into
local patches in both X and Y axes. \<ParamFile\> is the parameter file.

__SolvePar__ then takes "initial.MPI_\<Rank\>.dat" as input files and outputs "laplace.MPI_\<Rank\>.dat".

**Note that SolveParRow must be run with the same number of processes as the number of local patches across the rows and
columns.  This also means that the number of processes must be the same as the number of "initial.MPI_\<Rank\>.dat"
file, and also be the same as \<NumPatchInX\> * \<NumPatchInY\>.**

## (Optional) Combine solutions together
__Combine__ combines multiple \*.dat files into a single \*.dat file.

### Combine usage:

``` ./Combine <NumPatchInX> <NumPatchInY> <DstDatFile> <SrcDatFile0>...<SrcDatFileN> ```

Where \<NumPatchInX\> and \<NumPatchInY\> are integers (greater than 0), specifying how the global grid is partitioned into
local patches in both X and Y axes.
\<DstDatFile\> is the name of the output combined .dat file
\<SrcDatFile0\>...\<SrcDatFileN\> is a list of space-separated input files.

## Plot results
TODO

# Source file descriptions
TODO
