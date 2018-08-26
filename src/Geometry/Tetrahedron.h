#pragma once
#include"Node.h"
#include"Polyhedral.h"
#include<cassert>
class Tetrahedron:public Polyhedral
{
public:
	Tetrahedron(unsigned int id= INVALID_UNIT);
	void set_four_nodes(Node* n0, Node* n1, Node* n2, Node* n3);
	void set_id(unsigned int id) { this->id = id; }
	void set_marker(int marker) { this->marker = marker; }
	int get_marker() { return marker; }
	~Tetrahedron();
	
private:
	unsigned int id;
	int marker;
	bool ini_state;
};

