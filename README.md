# QuadriFlow: A Scalable and Robust Method for Quadrangulation

Source code for the paper:

Jingwei Huang, Yichao Zhou, Matthias Niessner, Jonathan Shewchuk and Leonidas Guibas. [**QuadriFlow: A Scalable and Robust Method for Quadrangulation**](http://stanford.edu/~jingweih/papers/quadriflow/quadriflow.pdf), The Eurographics Symposium on Geometry Processing (SGP) 2018.

## Processing Result

![QuadriFlow Results](https://github.com/hjwdzh/quadriflow/raw/master/img/result.jpg)

## Demo

The software supports cmake build for Linux/Mac/Windows systems. For linux and mac users, run **`sh demo.sh`** to build and try the QuadriFlow example, which converts `examples/Gargoyle_input.obj` to `examples/Gargoyle_quadriflow.obj`.

### Install

```
git clone --recursive -j8 git://github.com/hjwdzh/quadriflow
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=release
make -j
```

### QuadriFlow Software

We take a manifold triangle mesh `input.obj` and generate a manifold quad mesh `output.obj`. The face number increases linearly with the resolution controled by the user.

```
./quadriflow -i input.obj -o output.obj -f [resolution]
```

Here, the resolution is the desired number of faces in the quad mesh.

## Advanced Functions

### Min-cost Flow
By default, `quadriflow` uses the Boykov maximum flow solver from boost becuase it is faster.  To
enable the adaptive network simplex minimum-cost flow solver, you can enable the `-mcf` option:

```
./quadriflow -mcf -i input.obj -o output.obj -f [resolution]
```

### Sharp Preserving
By default, `quadriflow` does not explicitly detect and perserve the sharp edges in the model. To
enable this feature, uses

```
./quadriflow -sharp -i input.obj -o output.obj -f [resolution]
```

### SAT Flip Removal (Unix Only)
By default, `quadriflow` does not use the SAT solver to remove the flips in the interger offsets
map.  To remove the flips and guarantee a watertight result mesh, you can enable the SAT solver.
First, make sure that `minisat` and `timeout` is properly installed under your `${PATH}`.  The
former can be done by building `3rd/MapleCOMSPS_LRB/CMakeLists.txt` and copying `minisat` to `/usr/bin`.
In addition, `timeout` is included in coreutils. If you are using Mac, you can install it using
homebrew:
```
brew install coreutils
export PATH="/usr/local/opt/coreutils/libexec/gnubin:$PATH"
```

You can verify if those binaries are properly installed by executing
```
which minisat
which timeout
```

After that, you can enable SAT flip removal procedure by executing
```
./quadriflow -sat -i input.obj -o output.obj -f [resolution]
```

When using the SAT flip removal, we also suggest you enabling the verbose logging to understand
what is going on. You can build quadriflow with the following options:
```
cmake .. -DCMAKE_BUILD_TYPE=release -DBUILD_LOG=ON
```

### GUROBI Support (For Benchmark Purpose)

To use the Gurobi integer programming to solve the integer offset problem, you can build QuadriFlow with
```
cmake .. -DCMAKE_BUILD_TYPE=release -DBUILD_GUROBI=ON -DBUILD_LOG=ON
```
This override other solvers and should only be used for benchmark purpose.

## Dependencies
Require: Boost, GLM, Eigen
Optional: TBB, OpenMP

## Authors
- [Jingwei Huang](mailto:jingweih@stanford.edu)
- [Yichao Zhou](mailto:zyc@berkeley.edu)

&copy; Jingwei Huang, Stanford University

**IMPORTANT**: If you use this software please cite the following in any resulting publication:
```
@article {10.1111:cgf.13498,
    journal = {Computer Graphics Forum},
    title = {{QuadriFlow: A Scalable and Robust Method for Quadrangulation}},
    author = {Huang, Jingwei and Zhou, Yichao and Niessner, Matthias and Shewchuk, Jonathan Richard and Guibas, Leonidas J.},
    year = {2018},
    publisher = {The Eurographics Association and John Wiley & Sons Ltd.},
    ISSN = {1467-8659},
    DOI = {10.1111/cgf.13498}
}
```

## Triangle Manifold

In case you are dealing with a triangle mesh that is not a manifold, we implemented the software that converts any triangle mesh to watertight manifold. Please visit https://github.com/hjwdzh/Manifold for details.
