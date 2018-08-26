#include "Point.h"
Point::Point()
{
}

Point::~Point()
{
}


void Point::setPoint(double x0, double y0, double z0)
{
	x = x0;y = y0;z = z0;
}

double Point::distance(const Point & p1)const
{
	//double x0 = p1.x - x;
	//double y0 = p1.y - y;
	//double z0 = p1.z - z;
	return sqrt((p1.x - x )*(p1.x - x) +(p1.y - y)*(p1.y - y)+(p1.z - z)*(p1.z - z));
}

Point Point::operator+(const Point & v)const
{
	return Point(x+v.x,y+v.y,z+v.z);
}

Point Point::operator-(const Point & v)const
{
	return Point(x-v.x,y-v.y,z-v.z);
}

Point Point::operator*(const double a)const
{
	return Point(x*a,y*a,z*a);
}

double Point::operator*(const Point & v)const
{
	return (x*v.x+y*v.y+z*v.z);
}

Point Point::operator-() const
{
	return Point(-x,-y,-z);
}

void Point::reverse()
{
	x = -x;
	y = -y;
	z = -z;

}

bool Point::operator==(const Point & p)
{
	Point r = (*this) - p;
	return r.size()<TOLERANCE;
}

Point Point::operator+=(const Point & v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return Point(x,y,z);
}

Point Point::operator*=(double a)
{
	x *= a;
	y *= a;
	z *= a;
	return Point(x,y,z);
}

Point Point::operator/(double a)
{
	return Point(x/a,y/a,z/a);
}

Point Point::operator/=(double a)
{
	x = x / a;
	y = y / a;
	z = z / a;
	return Point(x,y,z);
}

double tripleProduct(const Point & v1, const Point & v2, const Point & v3)
{
	return (v1.x*v2.y*v3.z+v1.y*v2.z*v3.x+v1.z*v2.x*v3.y
		-v1.z*v2.y*v3.x-v1.x*v2.z*v3.y-v1.y*v2.x*v3.z);
}

double Point::size() const
{
	return sqrt(x*x+y*y+z*z);
}

Point Point::getUnit() const
{
	return Point(x / sqrt(x*x + y*y + z*z),y / sqrt(x*x + y*y + z*z), z / sqrt(x*x + y*y + z*z));
}

void Point::setUnit()
{
	double length=sqrt(x*x + y*y + z*z);
	x = x / length;
	y = y / length;
	z = z / length;
	return;
}
void Point::zero() { x = 0.; y = 0.; z = 0.; }
//projection onto plane pasing through p1,p2 and p3
Point Point::projection(const Point& p1, const Point& p2, const Point& p3) const
{
	//normal vector of the plane
	Point n = cross(p3-p1, p2-p1);
	//plane equation:Ax+By+Cz+D=0;
	double A=n.x,B=n.y,C=n.z;
	double D=-(A*p1.x+B*p1.y+C*p1.z);
	double t = -(D+A*x+B*y+C*z) / (A*A+B*B+C*C);
	return Point(A*t+x,B*t+y,C*t+z);
}

Point Point::projection_line(const Point& p1, const Point& p2)const
{
	Point t = (p2 - p1).getUnit();
	Point prj = p1 + (((*this) - p1)*t)*t;
	return prj;
}

//v1 cross v2
Point cross(const Point & v1, const Point & v2)
{
	double i, j, k;
	double x1 = v1.x, x2 = v2.x;
	double y1 = v1.y, y2 = v2.y;
	double z1 = v1.z, z2 = v2.z;
	i = y1*z2 - y2*z1;
	j = x2*z1 - x1*z2;
	k = x1*y2 - x2*y1;
	return Point(i, j, k);
}
//return unit vector of cross product
Point unitCross(const Point & v1, const Point & v2)
{
	double i, j, k;
	double x1 = v1.x, x2 = v2.x;
	double y1 = v1.y, y2 = v2.y;
	double z1 = v1.z, z2 = v2.z;
	i = y1*z2 - y2*z1;
	j = x2*z1 - x1*z2;
	k = x1*y2 - x2*y1;
	Point cro(i, j, k);
	return cro.getUnit();
}
Point operator*(double a, const Point & p)
{
	return Point(a*p.x,a*p.y,a*p.z);
}
