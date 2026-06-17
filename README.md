# Noise Inference from Single-Cell Expression (N.I.S.C.E.)

NISCE fit models of gene expression in single cells to estimate parameters from count data using maximum likelihood.
The input data is a set of counts (positive integers) for many genes (per rows) in many cells (per column). Counts can be raw read or UMI counts.
The underlying model supports heterogeneous cell coverage (normalization), and optionally dropout probabilities. Both Scaled-Beta (aka BetaPoisson) and Gamma (aka Negative Binomial) models are supported.
NISCE further allows partition models, with pre-defined cell clusters having distinct sets of parameters.

## Installation

### Dependencies

NISCE is written in C++20 and use builtin multithreading. It requires GCC (>= 12) or any compiler with support for `execution::par` (which is currently not the case of CLang).
NISCE also requires the libraries `boost::math_c99`, `TBB`, `Zlib` and Bio++ `bpp-core` (version $\geq 3.00$). 
Once the libraries are installed, NISCE can be compiled using the commands:

### Compilation

Clone the repository:

```bash
git clone https://github.com/jydu/nisce/
cd nisce
```

Prepare external build:

```bash
mkdir build
cd build
```

Run CMake:

```bash
CC=gcc-15 CXX=g++-15 cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
```

The CC and CXX variables are used to specify the compiler to use, in case GCC is not your default compiler. They can otherwise be omitted.
`CMAKE_INSTALL_PREFIX` tells where to install the compiled program. In this example, it will be installed in `$HOME/.local/bin/`.

Finally, compile and install NISCE:

```bash
make install
```

## Run the program

TODO
