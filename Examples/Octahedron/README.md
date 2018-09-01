https://github.com/zhong-yy/GraPly/blob/master/Examples/Octahedron/README.md

The file 'octahedron_model' describes the geometry and density contrast of a octahedral mass body.

The file 'vertex_edge_interior' contains 3 observation points: a corner point (0, 0, 0), an edge point (1, 0, 1) and an interior point (0,0,1).

The file 'profile1' contains observation points on a profile outside the octahedron, in a range of x=[-2,2] km, y=0.5 km, and z=-0.001 km. The file 'profile1' contains a measuring profile passing the octahedron, with x=[-2,2] km, y=0.5 km, and z=-0.001 km. The point interval is 0.1 km.

Assuming you have built the program 'GraPly' using makefile in the '/GraPly' directory, first copy 'GraPly' to the current directory.

1. To calculate the gravity field at a corner point, a edge point and an interior point, type
```
./GraPly octahedron_model vertex_edge_interior g_out_vertex_edge_interior g
```
The results are written to 'g_out_vertex_edge_interior'.

2. To calculate the gravity and gravity gradient tensor on profile 1 type
```
./GraPly octahedron_model profile1 g_out_profile1 g
./GraPly octahedron_model profile1 ggt_out_profile1 ggt
```
The gravity fields are output to 'g_out_profile1', and the gravity gradient tensors are output to 'ggt_out_profile1'.

3. To calculate the gravity and gravity gradient tensor on profile 2, type
```
./GraPly octahedron_model profile2 g_out_profile2 g
./GraPly octahedron_model profile2 ggt_out_profile2 ggt
```
The gravity fields are output to 'g_out_profile2', and the gravity gradient tensors are output to 'ggt_out_profile2'.

4. use python to plot results on profile1 and profile2 
```
python plot.py
```

Altanatively, you can use the script octahedron.sh to do all the above things.
