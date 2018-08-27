1. Copy 'GraPly' and 'Sites' to the current directory.
2. Download and install 3D mesh generator tetgen
http://www.wias-berlin.de/software/index.jsp?id=TetGen&lang=1
3. Copy 'tetgen' to the current directory (or set environment variable)
4. Run script 'mesh_generating.sh' to generate tetrahedral mesh and measuring points
```
sh mesh_generating.sh
```
5. Run script 'calculate.sh' to calculate the gravity field and gravity gradient tensor
```
sh calculate.sh
```
