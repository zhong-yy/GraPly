# calculate the gravity fields of the octahedron at a vertex point (0,0,0), an edge point (1,0,1),
# and a point inside the octahedron(0,0,1). Output results to file 'g_out_vertex_edge_interior'
./GraPly octahedron_model vertex_edge_interior g_out_vertex_edge_interior g


# calculate gravity and gravity gradient on profile 1, which is outside the octahedron
# output gravity field to file 'g_out_profile1', gravity gradient tensor field to file 'ggt_out_profile1'
./GraPly octahedron_model profile1 g_out_profile1 g
./GraPly octahedron_model profile1 ggt_out_profile1 ggt

# calculate gravity and gravity gradient on profile 2, which is passing through the octahedron
# output gravity field to file 'g_out_profile2', gravity gradient tensor field to file 'ggt_out_profile2'
./GraPly octahedron_model profile2 g_out_profile2 g
./GraPly octahedron_model profile2 ggt_out_profile2 ggt

#plot results on profile1 and profile2
python plot.py
