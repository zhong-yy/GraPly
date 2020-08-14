# GraPly
It's a C++ program for computing gravity fields and gravity gradient tensors (GGT) caused by polyhedrons with polynomial density contrasts (up to cubic order). Singularity-free analytical solutions are used in the calculation. For gravity, observation sites can be located inside, on, or outside the 3D mass bodies. For gravity gradient tensor, observation sites can be located everywhere except on the egdes or at the corners of polyhedrons (though on a facet the GGT values are evaluated as average of left-handed and righ-handed limits).


## 1 Building
'cd' to the 'GraPly' directory and type 'make' to compile.  If everything goes well, you will get three executable programs: **GraPly**, **GraTet** and **Sites**

- **GraPly** is a program to calculate the gravity and GGT of (convex) general polyhedrons using singularity-free closed-form solutions. 
- **GraTet** is a program for forward modelling of gravity fields and GGTs using tetrahedral unstructured grids.
- **Sites** generates regular grid of measuring points in a ASCII file.

To remove the program binaries and object files, just type 'make clean'.

## 2 How to use
### 2.1 GraPly
The command to use  '**GraPly**' is 
```
GraPly model_file observation_file output_file field_flag
```
The 'field_flag' can be 'g' or 'ggt'. If 'g' is used, the gravity field will be calculated; if 'ggt' is used, the gravity gradient tensor field will be calculated.

**Model file**
The 'model_file' is a file containing descriptions about polyhedrons. The format for the model file is 
```
*********************************************************************************
First line: <Number of polyhedrons n> 
Following lines: <Description of polyhedron 0>
                 <Description of polyhedron 1>
                 ...
                 <Description of polyhedron n-1>
*********************************************************************************
```
where the format for <descriptoin of polyhedron i\> is
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

Given the polynomial coefficients in the model file, the polynomial density contrast is
```
a000
+ a100*x   + a010*y     + a001*z
+ a200*x^2 + a020*y^2   + a002*z^2   + a101*x*z + a011*y*z + a110*x*y
+ a003*z^3 + a012*y*z^2 + a021*y^2*z + a030*y^3 + a102*x*z^2+ a111*x*y*z+ a120*x*y^2+ a201*x^2*z+ a210*x^2*y+ a300x^3
```
i.e.

![](http://latex.codecogs.com/gif.latex?\\lambda=a_{000}+a_{100}x+a_{010}y+a_{001}z+a_{200}x^2+a_{020}y^2+a_{002}z^2+a_{101}xz+a_{011}yz+a_{110}xy+a_{003}z^3+a_{012}yz^2+a_{021}y^2z+a_{030}y^3+a_{102}xz^2+a_{111}xyz+a_{120}xy^2+a_{201}x^2z+a_{210}x^2y+a_{300}x^3)

![](https://latex.codecogs.com/gif.latex?\dpi{400}\lambda=a_{000}+a_{100}x+a_{010}y+a_{001}z+a_{200}x^2+a_{020}y^2+a_{002}z^2+a_{101}xz+a_{011}yz+a_{110}xy+a_{003}z^3+a_{012}yz^2+a_{021}y^2z+a_{030}y^3+a_{102}xz^2+a_{111}xyz+a_{120}xy^2+a_{201}x^2z+a_{210}x^2y+a_{300}x^3)

In the model file, all polyhedrons, facets of a single polyhedron,and corners of a single facet are numbered from zero. The units for x, y, z are kilometer, and the unit for density contrast is kilogram per cubic meter.  

**Observation point file**
The format for files containing observation points is given as
```
*********************************************************************************
First line: <Number of points>
Following lines list coordinates of points:
            <index of this point> <x> <y> <z>
            ...
*********************************************************************************
```
The obervation points can be numbered from zero or one.  Comments in the files are prefixed by charactor '#'.

An example is given in the folder '/Examples/Octahedron/' .

### 2.2 GraTet
The command for '**GraTet**' is
```
GraTet configuration_file
```
The configuration file contains information about mesh files and density contrasts of different regions. The format is
```
*********************************************************************************
<input file name>
<observation point file name>
<order of polynomial density contrast> #0, 1, 2 or 3
<number of regions> # agree with the number of regions specified in *.poly file
<marker> #region marker
<polynomial coefficients of density contrast>
...      # next region
*********************************************************************************
```
The polynomial coefficients given in the file as
```
a000
a100  a010  a001
a200  a020  a002  a101  a011  a110
a003  a012  a021  a030  a102  a111  a120  a201  a210  a300
```
will define the following density contrast:
```
a000
+ a100*x   + a010*y     + a001*z
+ a200*x^2 + a020*y^2   + a002*z^2   + a101*x*z + a011*y*z + a110*x*y
+ a003*z^3 + a012*y*z^2 + a021*y^2*z + a030*y^3 + a102*x*z^2+ a111*x*y*z+ a120*x*y^2+ a201*x^2*z+ a210*x^2*y+ a300x^3
```

An example is given in the folder '/Examples/Tetrahedral_grid/' .

## 2.3 Sites
**Sites** is a tool to generate regular grid of measuring points. 

The command is
```
Sites x_start:x_interval:x_end y_start:y_interval:y_end z output_file_name
```

To generate a measuring profile with x=[0,10], y=0, z=-0.001 and 0.5 point interval, and write results to file 'out', type
```
./Sites 0:0.5:10 0 -0.001 out
```

To generate a grid with z=-0.001, x=[0,10], y=[0,10], point interval 0.5, and write results to file 'out', type
```
./Sites 0:0.5:10 0:0.5:10 -0.001 out
```

To generate a grid with z=[-5,1], x=[-5,5], y=[0,10], point interval 0.5, and write results to file 'out', type
```
./Sites -5:0.5:5 0:0.5:10 -5:0.5:1 out
```

## 3 Citation
Use the following citations to reference our codes:

- Zhengyong Ren, Chaojian Chen, Kejia Pan, Thomas Kalscheuer, Hansruedi Maurer, and Jingtian Tang. Gravity Anomalies of Arbitrary 3D Polyhedral Bodies with Horizontal and Vertical Mass Contrasts. Surveys in Geophysics, 38(2):479–502, 2017.  DOI: https://doi.org/10.1007/s10712-016-9395-x

- Zhengyong Ren, Yiyuan Zhong, Chaojian Chen, Jingtian Tang, and Kejia Pan. Gravity anomalies of arbitrary 3D polyhedral bodies with horizontal and vertical mass contrasts up to cubic order. Geophysics, 83(1):G1–G13, 2018. DOI: https://doi.org/10.1190/geo2017-0219.1

- Zhengyong Ren, Yiyuan Zhong, Chaojian Chen, Jingtian Tang, Thomas Kalscheuer, Hansruedi Maurer, and Yang Li. Gravity Gradient     Tensor of Arbitrary 3D Polyhedral Bodies with up to Third-Order Polynomial Horizontal and Vertical Mass Contrasts. Surveys in Geophysics, 39(5):901–935, 2018. ([Open access](https://link.springer.com/article/10.1007/s10712-018-9467-1))  DOI: https://doi.org/10.1007/s10712-018-9467-1 

If you are a BibTex or BibLaTex user, see [citation.bib](https://github.com/zhong-yy/GraPly/blob/master/References/citation.bib).

## Notices / todos
It may not work correctly when dealing with non-convex polyhedrons, because the automated calculation of facet normal vectors in the codes (Polyhedral.cpp, Integral.cpp) is invalid for non-convex polyhedron, although the closed-form formulae in the above cited papers are fit for both convex and non-convex polyhedrons. An alternative is to decompose the non-convex polyhedron into a couple of convex polyhedrons.
