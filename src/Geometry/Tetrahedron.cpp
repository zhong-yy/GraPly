#include "Tetrahedron.h"



Tetrahedron::Tetrahedron(unsigned int id)
{
	this->id = id;
	this->node_number = 4;
	this->face_number = 4;
	this->node.resize(4);
	for (int i = 0; i < node_number; i++) {
		node[i] = NULL;
	}

	this->face_node_number.resize(4);
	this->face_node.resize(4);
	this->face_normal_vector.resize(4);
	this->faces.resize(4);

	for (int i = 0; i < face_number; i++) {
		face_node_number[i] = 3;
		face_node[i].resize(3);
	}
	//face 0
	face_node[0][0] = 0;
	face_node[0][1] = 1;
	face_node[0][2] = 2;

	//face 1
	face_node[1][0] = 1;
	face_node[1][1] = 2;
	face_node[1][2] = 3;

	//face 2
	face_node[2][0] = 0;
	face_node[2][1] = 1;
	face_node[2][2] = 3;

	//face 3
	face_node[3][0] = 0;
	face_node[3][1] = 2;
	face_node[3][2] = 3;
	ini_state = false;
}
Tetrahedron::~Tetrahedron()
{
}

void Tetrahedron::set_four_nodes(Node * n0, Node * n1, Node * n2, Node * n3)
{
	node[0] = n0;
	node[1] = n1;
	node[2] = n2;
	node[3] = n3;
	this->sort_and_compute_normal_vector();
	this->set_faces_from_indices();
	ini_state = true;
}