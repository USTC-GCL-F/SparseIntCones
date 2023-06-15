# Practical Integer-Constrained Cone Construction for Conformal Parameterizations

This is a C++ implement of the following paper:

The code is written by Zheng Zhang using Microsoft Visual Studio 2019, and has passed the test on Windows 10.

## External Libraries
All the required libraries are the latest versions.

* [Mosek](https://www.mosek.com/)
* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/)
* [CGAL](https://www.cgal.org/)

## Build and Run

### build
```
cd IntCones
```
**Edit SparseCones/CmakeLists.txt** to set the solver and configuration path.(line 39, 46 and 47)

```
mkdir build && cd build
cmake ..
```

### run
```
exeSparseCones.exe example.obj output --normBound=sigma
```

### example
```
exeSparseCones.exe dog.obj dog --normBound=0.2
```

## Output Files
* ouput.txt: The running time of our algorithm and distortion are recorded.