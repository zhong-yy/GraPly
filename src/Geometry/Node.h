#pragma once
#include"Point.h"
#include"parameter.h"
class Node: public Point
{
public:
	Node();
	Node(double x, double y, double z, unsigned int id=-1, int marker=-99);
	~Node();
private:
	unsigned int id;
	int marker;
	//double rho;
};

