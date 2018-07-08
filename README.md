# QuadriFlow: A Scalable and Robust Method for Quadrangulation

Source code for the paper:

Jingwei Huang, Yichao Zhou, Matthias Niessner, Jonathan Shewchuk and Leonidas Guibas. [**QuadriFlow: A Scalable and Robust Method for Quadrangulation**](http://stanford.edu/~jingweih/papers/quadriflow/quadriflow.pdf), The Eurographics Symposium on Geometry Processing (SGP) 2018.

## Demo

The software supports cmake build for Linux/Mac/Windows systems. For linux and mac users, run **sh demo.sh** to build and try the QuadriFlow example, which converts examples/Gargoyle_input.obj to examples/Gargoyle_quadriflow.obj.

### Install

git clone --recursive -j8 git://github.com/hjwdzh/quadriflow

mkdir build

cd build

cmake ..

make

### QuadriFlow Software

We take a manifold triangle mesh "input.obj" and generate a manifold quad mesh "output.obj". The face number increases linearly with the resolution controled by the user.

./quadriflow -i input.obj -o output.obj -f [resolution]

## Dependencies
Require: Boost, GLM, Eigen
Optional: TBB, OpenMP

## Authors
- [Jingwei Huang](mailto:jingweih@stanford.edu)

&copy; Jingwei Huang, Stanford University


**IMPORTANT**: If you use this software please cite the following in any resulting publication:
```
@inproceedings{huang2018quadriflow,
  title={QuadriFlow: A Scalable and Robust Method for Quadrangulation},
  author={Huang, Jingwei and Zhou, Yichao and Nie{\ss}ner, Matthias and Shewchuk, Jonathan and Guibas, Leonidas},
  booktitle={Proceedings of the Eleventh Eurographics/ACMSIGGRAPH Symposium on Geometry Processing},
  pages={1--14},
  year={2018},
  organization={Eurographics Association}
}
```


