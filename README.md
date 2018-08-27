# GraPly
Gravity fields and gravity gradient tensors (GGT) caused by polyhedrons with polynomial density constrasts (up to cubic order). Singularity-free analytical solutions are used in the calculation. For gravity, observation sites can be located inside, on, or outside the 3D mass bodies. For gravity gradient tensor, observation sites can be located everywhere except on the egdes or at the corners of polyhedrons (though on a facet the GGT values are evaluated as average of left-handed and righ-handed limits).


## Building
'cd' to the 'GraPly' directory and type 'make' to compile.  If everything goes well, you will get three executable program: **GraPly**, **GraTet** and **Sites**

- **GraPly** is a program to calculate the gravity and GGT of (convex) general polyhedrons using singularity-free closed-form solutions. 
- **GraTet** is a program for forward modelling of gravity fields and GGTs using tetrahedral unstructured grids.
- **Sites** generates regular grid of measuring points in a text file.

To remove the program binaries and object files, just type 'make clean'.

## How to use
For **GraPly**, the command is 
```
GraPly model_file observation_file output_file field_flag
```
The 'field_flag' can be 'g' or 'ggt'.

The 'model_file' is a file containing descriptions about polyhedrons.
```
First line: <Number of polyhedrons n>
Following lines: <Description of polyhedron 0>
                 <Description of polyhedron 1>
                 ...
                 <Description of polyhedron n-1>
```
The format for <descriptoin of polyhedron i> is

```
One line: <index of this polyhedron> <number of nodes> <number of facets>
Following lines list nodes:
  <index of this node> <x> <y> <z>
  ...
Following lines list facets:
  <index of this facet> <number of corners> <corner 1> <corner 2> ...
Folowing lines list polynomial coefficients of polynomial density contrast:
  a000
  a100  a010  a001
  a200  a020  a002  a101  a011  a110
  a003  a012  a021  a030  a102  a111  a120  a201  a210  a300
```
The polynomial density contrast is



## Citation
- Zhengyong Ren, Chaojian Chen, Kejia Pan, Thomas Kalscheuer, Hansruedi Maurer, and Jingtian Tang. Gravity Anomalies of Arbitrary 3D Polyhedral Bodies with Horizontal and Vertical Mass Contrasts. Surveys in Geophysics, 38(2):479–502, 2017.
- Zhengyong Ren, Yiyuan Zhong, Chaojian Chen, Jingtian Tang, Thomas Kalscheuer, Hansruedi Maurer, and Yang Li. Gravity Gradient     Tensor of Arbitrary 3D Polyhedral Bodies with up to Third-Order Polynomial Horizontal and Vertical Mass Contrasts. Surveys in Geophysics, 39(5):901–935, 2018. ([Open access](https://link.springer.com/article/10.1007/s10712-018-9467-1))
- Zhengyong Ren, Yiyuan Zhong, Chaojian Chen, Jingtian Tang, and Kejia Pan. Gravity anomalies of arbitrary 3D polyhedral bodies with horizontal and vertical mass contrasts up to cubic order. Geophysics, 83(1):G1–G13, 2018.


## Notes
It may not work correctly when dealing with non-convex polyhedrons, because the automated calculation of facet normal vectors in the codes (Polyhedral.cpp, Integral.cpp) is invalid for non-convex polyhedron.
