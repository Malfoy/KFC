# KFC
kmer week experiment

[![Build Status](https://travis-ci.com/Malfoy/KFC.svg?branch=master)](https://travis-ci.com/Malfoy/KFC)


## Building

Prior to any compilation, please check your gcc version. It must be 8+.

To build a debug version (default in this early development stage):
``` bash
git clone https://github.com/Malfoy/KFC.git
cd KFC
make -j
```

For an optimized build add `DEBUG=0` to the make command:
```
make -j DEBUG=0
```

For building with optimizations while keeping asserts enabled, use `make -j DEBUG=0 ASSERTS=1`

## Development

If you have `clang-format` installed on your system, please setup the pre-commit hook:
```
.scripts/install_hooks.sh
```

## Benchmarks

Please use the snakefile to compare latest KMC version with KFC.
It's under construction.

