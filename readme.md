# BlackHole (Optimized)

An optimized implementation of the **BlackHole** community detection algorithm based on the original implementation provided by the authors of the paper:

> **BlackHole: Robust Community Detection Inspired by Graph Drawing**
> S. Lim, J. Kim, and J.-G. Lee, ICDE 2016.

The original implementation is available at:
https://github.com/kaist-dmlab/BlackHole

## Overview

This repository contains an optimized version of the original BlackHole implementation with a focus on improving runtime performance while preserving the original algorithmic behavior.

## How to run
The project contains a CMakeLists.txt in the source folder which needs to be built with CMake.

## How to use
In "main.cpp" the main algorithmic parameters are defined and can be tuned in addition to the input and output files. The input requires an edgelist starting with vertex indicies counting from 0. The output are the indices of the vertices in each community. The attraction and repulsion exponents are controlled in "projection_simulation.hpp" and may be tweaked. The default relative paths assume the build folder is placed in the main directory. Additionally the output folder must also be created manually or the path changed to an existing folder.

## Improvements

The following changes have been made:

* Performed numerous low-level performance optimizations throughout the codebase.
* Fixed a number of memory leaks.
* Improved documentation.
* Refactored and cleaned up the project structure to improve maintainability and readability.

## Performance

The optimizations provide significant speedups over the original implementation.

| Number of Vertices | Base: Projecting | Base: Clustering | Base: Total | Improved: Projecting | Improved: Clustering | Improved: Total |   Speedup |
| -----------------: | ---------------: | ---------------: | ----------: | -------------------: | -------------------: | --------------: | --------: |
|                10k |           60.0 s |            9.5 s |      69.5 s |               27.0 s |                1.4 s |          28.4 s | **2.45×** |
|                50k |          536.5 s |          213.0 s |     749.5 s |              189.2 s |               35.9 s |         225.1 s | **3.33×** |
|               100k |         1261.2 s |          862.2 s |    2123.4 s |              422.5 s |              160.2 s |         582.7 s | **3.64×** |

## Citation

If you use this implementation in your research, please cite the original BlackHole paper:

```bibtex
@inproceedings{blackhole,
  author    = {Sungsu Lim and Junghoon Kim and Jae-Gil Lee},
  title     = {BlackHole: Robust Community Detection Inspired by Graph Drawing},
  booktitle = {2016 IEEE 32nd International Conference on Data Engineering (ICDE)},
  pages     = {25--36},
  year      = {2016},
  doi       = {10.1109/ICDE.2016.7498226}
}
```

## License

This repository is based on the original BlackHole implementation. Please refer to the original project for its licensing terms and ensure compliance when using or redistributing this code.
