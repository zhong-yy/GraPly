#include "Dyadic.h"



Dyadic::Dyadic()
{
	xx = 0; yx = 0; zx = 0;
	xy = 0; yy = 0; zy = 0;
	xz = 0; yz = 0; zz = 0;
}

Dyadic::Dyadic(Point a, Point b)
{
	xx = a.x*b.x; yx = a.y*b.x; zx = a.z*b.x;
	xy = a.x*b.y; yy = a.y*b.y; zy = a.z*b.y;
	xz = a.x*b.z; yz = a.y*b.z; zz = a.z*b.z;
}

Dyadic::Dyadic(double xx, double yx, double zx, double xy, double yy, double zy, double xz, double yz, double zz)
{
	this->xx = xx; this->yx = yx; this->zx = zx;
	this->xy = xy; this->yy = yy; this->zy = zy;
	this->xz = xz; this->yz = yz; this->zz = zz;
}


Dyadic::~Dyadic()
{
}

Dyadic Dyadic::operator+(const Dyadic & d)
{

	return Dyadic(this->xx+d.xx,this->yx+d.yx,this->zx+d.zx,
		this->xy+d.xy,this->yy+d.yy,this->zy+d.zy,
		this->xz+d.xz,this->yz+d.yz,this->zz+d.zz);
}

Dyadic Dyadic::operator-(const Dyadic & p) const
{
	return Dyadic(xx-p.xx,yx-p.yx,zx-p.zx,xy-p.xy,yy-p.yy,zy-p.zy,xz-p.xz,yz-p.yz,zz-p.zz);
}

Dyadic Dyadic::operator+=(const Dyadic & d)
{
	xx += d.xx;
	yx += d.yx;
	zx += d.zx;
	xy += d.xy;
	yy += d.yy;
	zy += d.zy;
	xz += d.xz;
	yz += d.yz;
	zz += d.zz;
	return Dyadic(xx,yx,zx,xy,yy,zy,xz,yz,zz);
}

Dyadic Dyadic::operator*(const double a) const
{
	return Dyadic(xx*a,yx*a,zx*a,xy*a,yy*a,zy*a,xz*a,yz*a,zz*a);
}

// set to zero
void Dyadic::set()
{
	xx = 0; yx = 0; zx = 0;
	xy = 0; yy = 0; zy = 0;
	xz = 0; yz = 0; zz = 0;
}

//show
void Dyadic::display()
{
	cout << xx << "," << yx << "," << zx << '\n'
		<< xy << "," << yy <<","<< zy << '\n'
		<< xz << "," << yz << "," << zz << '\n';
}

//get a line of elements
Point Dyadic::operator[](int i)const
{
	Point a;
	assert(i >= 0 && i <= 2);
	switch (i) {
	case 0:
		a.setPoint(xx, yx, zx); break;
	case 1:
		a.setPoint(xy, yy, zy); break;
	case 2:
		a.setPoint(xz, yz, zz); break;
	}
	return a;
}



Dyadic operator*(double a, const Dyadic & d)
{
	return Dyadic(a*d.xx,a*d.yx,a*d.zx,a*d.xy,a*d.yy,a*d.zy,a*d.xz,a*d.yz,a*d.zz);
}
