#include "tet_mesh.h"



tet_mesh::tet_mesh()
{
}


tet_mesh::~tet_mesh()
{
}

void tet_mesh::clear()
{
	assert(tets.size() == n_tets);
	assert(mesh_nodes.size() == n_nodes);
	for (unsigned int i = 0; i < n_tets; i++) {
		if (tets[i] != NULL){
			delete tets[i];
			tets[i] = NULL;
		}
	}
	this->tets.clear();

	for (unsigned int j = 0; j < n_nodes; j++) {
		if (mesh_nodes[j] != NULL) {
			delete mesh_nodes[j];
			mesh_nodes[j] = NULL;
		}
	}
	this->mesh_nodes.clear();
}

void tet_mesh::read_mesh(string grid_name)
{
	string name_node, name_ele;
	string dummy = grid_name;
	name_node = dummy + string(".node");
	name_ele = dummy + string(".ele");

	ifstream node_stream(name_node.c_str());
	ifstream ele_stream(name_ele.c_str());
	assert(node_stream.good());
	assert(ele_stream.good());
	cout << "Reading mesh with tetrahedal elements ....\n";
	nodes_in(node_stream);
	tets_in(ele_stream);

	return;

}

void tet_mesh::nodes_in(std::istream & node_stream)
{
	assert(node_stream.good());// Check input buffer
	unsigned int dimension = 0, nAttributes = 0, BoundaryMarkers = 0;
	node_stream >> n_nodes          // Read the number of nodes from the stream
		>> dimension        // Read the dimension from the stream
		>> nAttributes      // Read the number of attributes from stream
		>> BoundaryMarkers; // Read if or not boundary markers are included in *.node (0 or 1)
	cout << "Number of nodes: " << n_nodes<<endl;
	assert(dimension == 3);
	assert(n_nodes>0);
	assert(nAttributes == 0);
	assert(BoundaryMarkers == 0);//boundary marker is currently not necessary
	this->mesh_nodes.clear();
	this->mesh_nodes.resize(n_nodes);
	unsigned int node_lab = 0;
	double x, y, z;
	int marker;
	for (unsigned int i = 0; i<n_nodes; i++)
	{
		// Check input buffer
		assert(node_stream.good());
		node_stream >> node_lab  // node number
			>> x         // x-coordinate value
			>> y         // y-coordinate value
			>> z;		 // z-coordinate value
		Node* newnode = new Node(x, y, z);
		this->mesh_nodes[i] = newnode;
	}
	return;

}

void tet_mesh::tets_in(std::istream & ele_stream)
{
	// Check input buffer
	assert(ele_stream.good());
	unsigned int tet_lab = 0, n_nodes_tet = 0, nAttri = 0;
	//n_nodes_tet: number of nodes associated with each element
	this->n_tets = 0;
	int tet_marker = 0;
	ele_stream >> n_tets     // Read the number of trirahedrons from the stream.
		>> n_nodes_tet // Linear or second Tet(4 or 10, defaults to 4).
		>> nAttri;     // Read the number of  region attributes from stream.
	cout << "Number of elements: " << n_tets<<endl;
	assert(n_tets>0);
	assert(n_nodes_tet == 4);
	assert(nAttri == 1);

	this->tets.clear();
	this->tets.resize(n_tets);

	for (unsigned int i = 0; i<n_tets; i++)
	{
		assert(ele_stream.good());
		this->tets[i] = new Tetrahedron;
		// Read the label
		ele_stream >> tet_lab;
		this->tets[i]->set_id(tet_lab);
		// Read node labels
		for (unsigned int j = 0; j<4; j++)
		{
			unsigned long int node_label;
			ele_stream >> node_label;
			// Assign node to trient
			this->tets[i]->set_node(j, this->mesh_nodes[node_label]);
		}
		this->tets[i]->init();//sort the nodes, and compute normals, and set face nodes pointers

		// Read attributes from the stream.
		ele_stream >> tet_marker;
		this->tets[i]->set_marker(tet_marker);// Only 1 atrri.
	}
	return;



}