# Info
A part of AOBench, a benchmark by Syoyo Fujita to Futhark, a data-parallel research language.
Adapted from https://github.com/ispc/ispc/blob/main/examples/cpu/aobench/ao.ispc

The `aobench_novec.fut` variant does _not_ use the Futhark vector library, whereas the `aobench_withvec.fut` variant does. I got much better performance using the former.

# Build
1. Install futhark
2. `futhark pkg sync`

# Run
To make an image, use `futhark literate --backend=<some backend> aobench_novec.fut`.

To run as a benchmark, use `futhark bench --backend=<some backend> aobench_novec.fut`

# Pretty
![](pic.png)
