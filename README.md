# `sun-hubbard`: Determinant Quantum Monte Carlo Package

A high-performance Determinant Quantum Monte Carlo (DQMC) package for simulating interacting fermionic lattice models, including $\mathrm{SU}(N)$ Hubbard model and spinless $t$-$V$ model. 
- Supports multiple lattice geometries: square, honeycomb, and cubic
- Features full MPI parallelization

## Project Structure

- `src/`: Source code files implementing core DQMC algorithms
- `lib/`: Library files containing random number generation, matrix operations, and FFT subroutines
- `utility/`: Utility scripts supporting DQMC simulations
- `examples/`: Example running scripts and data analysis tools
- `bin/` (if any): Compiled executables

## Compilation

### Prerequisites
Intel® Fortran Essentials or Intel® oneAPI HPC Toolkit ([Download from Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.k2rl4h))

### Build Instructions
```bash
# First compile the libraries
cd lib
cp make.sys.mpiifx make.sys
make

# Then compile the main program
cd ../src
cp make.sys.[lattice_type] make.sys  # Choose from: square, honeycomb, cubic
make ft # finite-temperature DQMC algorithm
make proj # projective DQMC algorithm
make # build the above two
```

## Running

### Local Execution
```bash
# Edit dqmc.in file to set simulation parameters
cp dqmc.in.sample.[lattice_type] dqmc.in

# Run the simulation
./dqmc_{ft,proj}
# Run the simulation with MPI
mpirun -np [number_of_cores] ./dqmc_{ft,proj}
```

Results will be written to `dqmc.out` and various data files for further analysis.

### HPC Cluster Execution
```bash
cd examples/run_honeycomb
# Check and modify parameters in cal_para.sh
vim cal_para.sh
./run_sjtusy.sh  # submit jobs to SJTU Siyuan cluster
```

## Preprocessor Options

The package behavior can be customized by modifying preprocessor options in the `make.sys` file under the `RUNMODE` variable.

### Update Algorithms

The following update algorithms are available for DQMC simulations:

- **Fast Update** (default): Efficient standard algorithm when no special update options are specified
- **Delay Update**: Add `-DDELAY` to enable the delay update algorithm (Ref: F. Sun and X. Y. Xu, Phys. Rev. B 109, 235140 (2024).)
- **Submatrix Update**: Add `-DSUBMATRIX` to enable the submatrix update algorithm (Ref: F. Sun and X. Xiao Yan, SciPost Physics 18, 055 (2025).)
  - For projective DQMC, the **Submatrix-LR** algorithm is used by default. Add `-DDELAYG` along with `-DSUBMATRIX` to use the deprecated **Submatrix-G** algorithm instead. 