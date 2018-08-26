#pragma once
#include"Point.h"
#include<cassert>
class Dyadic
{
public:
	Dyadic();
	Dyadic(Point a, Point b);
	Dyadic(double xx, double yx, double zx, double xy, double yy, double zy, double xz, double yz, double zz);
	
	~Dyadic();
	//operations of dyadics
	Dyadic operator+(const Dyadic& d);
	
	//subtract
	Dyadic operator - (const Dyadic& p)const;
	
	//plus
	Dyadic operator += (const Dyadic& d);
	
	//multiply
	Dyadic operator * (const double a)const;
	
	//multiply
	friend Dyadic operator*(double a, const Dyadic& d);
	
	//set to zero
	void set();
	
	//show the values
	void display();
	
	Point operator[](int i)const;
	//entries
	double xx, yx, zx;
	double xy, yy, zy;
	double xz, yz, zz;
};

