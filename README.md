# fc2d_hps: ForestClaw2D Hierarchical Poincaré-Steklov Implementation

## Overview

This repo contains the ongoing implementation of the Hierarchical Poincaré-Steklov (HPS) Method into `ForestClaw`.

The HPS method is a direct solver for elliptic PDEs.

NOTE: This is not currently meant to be a stand-alone code; it is meant to be cloned into a `ForestClaw` instance under `${FORESTCLAW_ROOT}/src/solvers`.

## Features

Current implementation includes a serial, uniform implementation on a `p4est` grid.

Current patch solvers include `FISHPACK90` for solving Poisson's equation on a rectangular domain.

Future implementations include adaptive and parallel implementations of the HPS method, both on CPUs and GPUs using `PETSc` and `Kokkos`.

Future patch solvers include a variable coefficient Poisson solver and dispersive Shallow Water Equations (Serre-Green-Naghdi model) for 2D simulations.

## Installation

As mentioned above, this is meant to be included as part of a `ForestClaw` installation. There are uses of `ForestClaw` functions throughout the HPS functions. A stand-alone version is in the works as well.

To install, first clone `ForestClaw`: [ForestClaw on Github](https://github.com/ForestClaw/forestclaw). Then, switch to the `develop` branch; all `fc2d_hps` work is being done in tandem to the `develop` branch of `ForestClaw`:

```bash
git clone https://github.com/ForestClaw/forestclaw
git checkout develop
```

After cloning `ForestClaw` but before configuration and compiling, clone this repo into the solvers directory in `ForestClaw`:

```bash
cd ${FORESTCLAW_ROOT}/src/solvers
git clone https://github.com/ForestClaw/fc2d_hps
```

Once `fc2d_hps` is a directory in the `ForestClaw` solvers directory, navigate back to `${FORESTCLAW_ROOT}` and configure with `cmake`. To use `fc2d_hps` with `ForestClaw`, you need specific flags set when configuring `cmake`, namely:

```bash
-Dclawpatch=on
-Dmpi=on
-Dhps=on
```

See the `ForestClaw` documentation on other installation configuration. When I configure, I use the following script in a build directory parallel to the `ForestClaw` source directory:

```bash
FCLAW=../forestclaw-source

cmake configure \
	-DCMAKE_Fortran_COMPILER=gfortran \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_CXX_COMPILER=mpic++ \
	-Dclawpatch=on \
	-Dmpi=on \
  	-Dhps=on \
	${FCLAW}
```

Once configured properly, use `make` to compile:

```bash
make -j
```

You can run all unit tests via `make` as well:

```bash
make -j test
```

## Usage

Currently, there are two executables to work with, `hello` and `simple`. They are found under `${FORESTCLAW_ROOT}/src/solvers/fc2d_hps/examples/`

### `hello`

The `hello` program simply imports all necessary headers and prints out a welcoming message to check for all libraries.

### `simple`

The `simple` program solves a Poisson equation on a 2D rectangular domain. All options are found in the `fclaw_options.ini` file. To run:

```bash
./simple fclaw_options.ini
```