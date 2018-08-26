#include "Node.h"



Node::Node()
{
}

Node::Node(double x, double y, double z, unsigned int id, int marker):Point(x,y,z)
{
	this->id = id;
	this->marker = marker;
}



Node::~Node()
{
}
