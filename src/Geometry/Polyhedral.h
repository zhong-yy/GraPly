#pragma once
#include<vector>
#include<algorithm>
#include<cassert>
#include"Point.h"
#include"Node.h"
class Polyhedral
{
public:
	Polyhedral();
	Polyhedral(int n);
	//Polyhedral(unsigned int node_num, unsigned int face_num);


	~Polyhedral();

	void init();

	void sort_and_compute_normal_vector();

	//this need to be called after calling "sort_and_compute_normal_vector()";
	void set_faces_from_indices();


	void set_node(unsigned int i, Node* node);
	

public:
	/*********geometrical parameters************/
	unsigned int node_number;
	unsigned int face_number;

	vector<Node*> node;

	vector<int> face_node_number;// number of nodes on each face
	vector<vector<unsigned> > face_node;//indices for nodes on each face
	vector<vector<Node*> > faces;
	vector<Point> face_normal_vector;
	//Point center;
	


	/*********physical parameters************/

	double const_density;
	double lin[3];
	double qua[6];
	double cub[10];
	/* rho=const_density
	 *     +lin[0]*x+lin[1]*y+lin[2]*z
	 *     +qua[0]*x^2+qua[1]*y^2+qua[2]*z^2+qua[3]*xz+qua[4]*yz+qua[5]*xy
	 *     +(cub[0]*z^3+cub[1]*yz^2+cub[2]*y^2*z+cub[3]*y^3)
	 *     +(cub[4]x*z^2+cub[5]*xyz+cub[6]xy^2)+(cub[7]*x^2*z+cub[8]*x^2*y)+cub[9]*x^3
	*/
protected:
	bool ini_state;

};

