# KFC
kmer week experiment


## Building

To build a debug version (default in this early developpment stage):
``` bash
git clone https://github.com/Malfoy/KFC.git
cd KFC
make -j
```

For an optimized build add `DEBUG=0` to the make command:
```
make DEBUG=0 -j
```

For building with optimizations while keeping asserts enabled, use `DEBUG=0 ASSERTS=1`
