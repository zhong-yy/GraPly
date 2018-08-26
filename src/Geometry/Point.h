#pragma once
/***************************************************
*Point
*  by CSU,ZHONG YIYUAN
*     2017,1,18
**************************************************/
#include<cmath>
#include<iostream>
#include"parameter.h"
using namespace std;
class Point;
Point operator*(double a, const Point & p);
Point cross(const Point& v1, const Point& v2);
Point unitCross(const Point& v1, const Point& v2);
class Point
{
public:
	Point();
	~Point();
	//constructer
	Point(double x1, double y1, double z1) :x(x1), y(y1), z(z1) {}


	//reset the point/vector value
	void setPoint(double x0, double y0, double z0);

	//distance from p1
	double distance(const Point& p1)const;
	
	/***********operator*****************/
	//addition and subtraction
	Point operator + (const Point& v)const;
	Point operator - (const Point& v)const;

	//dot product and scalar multiplication
	Point operator * (const double a)const;
	friend Point operator*(double a,const Point& p);
	double operator * (const Point& v)const;

	Point operator -()const;
	void reverse();

	//cross product
	friend Point cross(const Point& v1, const Point& v2);
	friend Point unitCross(const Point& v1, const Point& v2);

	//compound assignment
	bool operator ==(const Point& p);
	Point operator += (const Point& v);
	Point operator *= (double a);

	//divided by a real
	Point operator/(double a);
	Point operator/=(double a);

	// triple product
	friend double tripleProduct(const Point& v1, const Point& v2,const Point& v3);


	/*******other********/
	//get the length of vector (x,y,z)
	double size()const;

	Point getUnit()const;//get the unit vector
	void setUnit();//to set itself to be a unit vector

	void zero();//reset

	//get the projection onto the plane defined by p1,p2 and p3
	Point projection(const Point& p1, const Point& p2, const Point& p3) const;
	
	//get the projection point on to the straight line passing through  points p1 and p2
	Point projection_line(const Point& p1, const Point& p2)const;
	
	//displaying the coordinate of the vector or point, which is used to debug
	void display() const{
		cout << "(" << x << "," << y << "," << z << ")\n";
	}
public:
	//components
	double x, y, z;
};