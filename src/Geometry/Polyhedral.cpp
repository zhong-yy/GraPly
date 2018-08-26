#include "Polyhedral.h"

Polyhedral::Polyhedral() { ini_state = false; }

Polyhedral::Polyhedral(int n) : face_number(n)
{
	face_normal_vector.resize(n);
	ini_state = false;
}

Polyhedral::~Polyhedral()
{
}
void Polyhedral::init()
{
	this->sort_and_compute_normal_vector();
	this->set_faces_from_indices();
	ini_state = true;
}
//compute unit normal vectors and sort
void Polyhedral::sort_and_compute_normal_vector()
{
	Point center(0, 0, 0);
	for (int i = 0; i < node_number; i++)
	{
		center += *node[i];
	}
	center = center / node_number;

	Point a, b; //temporary variables

	for (int i = 0; i < face_number; i++)
	{
		//https://stackoverflow.com/questions/32274127/how-to-efficiently-determine-the-normal-to-a-polygon-in-3d-space
		//https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
		Point *current, *next;
		Point sum_normal(0.0, 0.0, 0.0);
		int nodeN = face_node_number[i];
		/*Newell method*/
		for (int k = 0; k < nodeN; k++)
		{
			current = node[face_node[i][k]];
			next = node[face_node[i][(k + 1) % nodeN]];
			sum_normal.x += (current->y - next->y) * (current->z + next->z);
			sum_normal.y += (current->z - next->z) * (current->x + next->x);
			sum_normal.z += (current->x - next->x) * (current->y + next->y);
		}
		
		sum_normal.setUnit();
		face_normal_vector[i] = sum_normal;

		// a = *node[face_node[i][1]] - *node[face_node[i][0]];
		// b = *node[face_node[i][2]] - *node[face_node[i][0]];
		// face_normal_vector[i] = unitCross(a, b);
		double flag = face_normal_vector[i] * (*node[face_node[i][0]] - center);
		if (flag < 0.0)
		{
			std::reverse(face_node[i].begin(), face_node[i].end());
			face_normal_vector[i].reverse();
		}
	}
	return;
}

void Polyhedral::set_faces_from_indices()
{
	faces.resize(face_number);
	for (int i = 0; i < face_number; i++)
	{
		faces[i].resize(face_node_number[i]);

		for (int k = 0; k < face_node_number[i]; k++)
		{
			faces[i][k] = node[face_node[i][k]];
		}
	}
}

void Polyhedral::set_node(unsigned int i, Node *node)
{
	assert(i < node_number);
	assert(node != NULL);
	this->node[i] = node;
}
