#include "Integral.h"

Integral::Integral()
{
}


Integral::~Integral()
{
}

double Integral::Line_R1(const Point& observation, const Point& v0, const Point& v1)
{
	//projection onto the line v0-v1 from the observation point 
	Point prj = observation.projection_line(v0, v1);
	Point t = (v1 - v0).getUnit();//unit tangent vector
	double D0 = observation.distance(prj); //distance from D0 to the line
	double s1, s0, R1, R0;
	s1 = t*(v1 - observation);//parameter
	s0 = t*(v0 - observation);
	R1 = observation.distance(v1);
	R0 = observation.distance(v0);
	double itg=0.0;
	//Table of integrals, series and products,equation(2.262.1),(Gradshteyn,1994)
	if (abs(s0 + R0)< TOLERANCE) {
		itg = 0.5*(s1*R1 - s0*R0);
	}
	else {
		itg=0.5*(D0*D0*log((s1+R1)/(s0+R0))+s1*R1 - s0*R0);
	}
	return itg;
	
}

double Integral::Line_R3(const Point& observation, const Point& v0, const Point& v1)
{
	//projection onto the line v0-v1 from the observation point 
	Point prj = observation.projection_line(v0, v1);
	Point t = (v1 - v0).getUnit();//unit tangent vector
	double D0 = observation.distance(prj);//distance from D0 to the line
	double s1, s0, R1, R0;
	s1 = t*(v1 - observation);//parameter s
	s0 = t*(v0 - observation);
	R1 = observation.distance(v1);
	R0 = observation.distance(v0);
	double itg = 0.0;
	//Table of integrals, series and products,equation(2.260.2),(Gradshteyn,1994)
	if (abs(s0 + R0) < TOLERANCE) {
		itg = 0.25*(s1*R1*R1*R1-s0*R0*R0*R0);
	}
	else {
		double itg_R1 = 0.5*(D0*D0*log((s1 + R1) / (s0 + R0)) + s1*R1 - s0*R0);
		itg = 0.25*(s1*R1*R1*R1 - s0*R0*R0*R0)
			+0.75*(D0*D0*itg_R1);
	}
	return itg;
}

double Integral::Line_R_1(const Point& observation, const Point& v0, const Point& v1)
{
	Point t = (v1 - v0).getUnit();
	Point prj = observation.projection_line(v0, v1);
	double L0 = observation.distance(prj);
	double s0, s1, R0, R1,result;
	s0 = (v0 - observation)*t;
	s1 = (v1 - observation)*t;
	R0 = v0.distance(observation);
	R1 = v1.distance(observation);
	
	if (abs(s0 + R0) < TOLERANCE) {
		result = abs(log(s1 / s0));
		assert(s0*s1 > 0);
	}
	else {
		result = log((R1 + s1) / (R0 + s0));
	}


	return result;
}

double Integral::Line_Bj1(const Point& observation, const Point& v0, const Point& v1, const Point& unit_normal)
{
	Point t = (v1 - v0).getUnit();
	Point mj = unitCross(t,unit_normal);
	double mj_value = (v0 - observation)*mj;
	if (abs(mj_value) < TOLERANCE) {
		return 0.0;
	}
	double h = (observation - v0)*unit_normal;
	double s0, s1, R0, R1;
	s1 = t*(v1 - observation);//parameter
	s0 = t*(v0 - observation);
	R1 = observation.distance(v1);
	R0 = observation.distance(v0);
	double bj1 = abs(h)*(atan((abs(h)*s1) / (mj_value*R1)) - atan((abs(h)*s0) / (mj_value*R0)))
		+ mj_value*log((s1 + R1) / (s0 + R0));
	//here,the case of mj_value=0 has been handled before
	return bj1;
}

double Integral::Line_Bj3(const Point& observation, const Point& v0, const Point& v1, const Point& unit_normal)
{
	Point t = (v1 - v0).getUnit();
	Point mj = cross(t, unit_normal);
	double mj_value = (v0 - observation)*mj;
	
	//judge if the mj_value is equal to zero
	if (abs(mj_value) < TOLERANCE) {
		return 0.0;
	}
	double s0, s1, R0, R1,h;
	
	//local coordinate s
	s1 = t*(v1 - observation);
	s0 = t*(v0 - observation);
	R1 = observation.distance(v1);
	//distance R
	R0 = observation.distance(v0);
	
	//height from observation point to the plane containing the face
	h = abs((observation - v0)*unit_normal);
	double bj3 = h*h*h*(atan((h*s1) / (mj_value*R1)) - atan((h*s0) / (mj_value*R0)))
		+ 0.5*mj_value*(3.0 * h*h + mj_value*mj_value)*log((s1 + R1) / (s0 + R0))
		+ 0.5*mj_value*(s1*R1 - s0*R0);
	//here,the case of mj_value=0 has been handled before
	return bj3;
}

Point Integral::Line_rR1(const Point& ob, const Point& v0, const Point& v1)
{
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);//projection on edge from point ob 
	Point L = prj - ob;
	double D0 = L.size();
	double s1, s0, R1, R0;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = ob.distance(v0);
	R1 = ob.distance(v1);
	double temp1, temp2;
	//cout << "here" << endl;
	if (abs(s0 + R0) < TOLERANCE)
		temp1 =0.5* (s1*R1 - s0*R0);
	else {
		temp1 = 0.5*(D0*D0 * log((s1 + R1) / (s0 + R0)) + s1*R1 - s0*R0);
	}
	temp2 = (R1*R1*R1 - R0*R0*R0) / 3.0;
	return (L*temp1+t*temp2);
}

Point Integral::Line_rR_1(const Point& ob, const Point& v0, const Point& v1)
{
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	Point result;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	if (abs(s0 + R0) < TOLERANCE) {
		result = L*abs(log(s1 / s0)) + t*(R1 - R0);
	}
	else {
		result = L*log((s1 + R1) / (s0 + R0)) + t*(R1 - R0);
	}
	return result;
}

double Integral::Line_x2R(const Point& ob, const Point& v0, const Point& v1)
{
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	double result = 0.0;
	if (abs(s0 + R0) < TOLERANCE) {
		result = (t.x)*(t.x)*0.25 * (s1*R1*R1*R1 - s0*R0*R0*R0)
			+ 2.0*(L.x)*(t.x)*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.x)*(L.x)*(s1*R1 - s0*R0);
	}
	else {
		result = (t.x)*(t.x)*0.125*(2*(s1*R1*R1*R1 - s0*R0*R0*R0) -L0*L0*(s1*R1-s0*R0)-L0*L0*L0*L0*log((s1+R1)/(s0+R0)))
			+2.0*(L.x)*(t.x)*(R1*R1*R1-R0*R0*R0)/3.0
			+0.5*(L.x)*(L.x)*(s1*R1-s0*R0+L0*L0*log((s1+R1)/(s0+R0)));
	}
	return result;
}

double Integral::Line_y2R(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	if (abs(s0 + R0) < TOLERANCE) {
		result = (t.y)*(t.y)*0.25 * (s1*R1*R1*R1 - s0*R0*R0*R0)
			+ 2.0*(L.y)*(t.y)*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.y)*(L.y)*(s1*R1 - s0*R0);
	}
	else {
		result = (t.y)*(t.y)*0.125*(2 * (s1*R1*R1*R1 - s0*R0*R0*R0) - L0*L0*(s1*R1 - s0*R0) - L0*L0*L0*L0*log((s1 + R1) / (s0 + R0)))
			+ 2.0*(L.y)*(t.y)*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.y)*(L.y)*(s1*R1 - s0*R0 + L0*L0*log((s1 + R1) / (s0 + R0)));
	}
	return result;
}

double Integral::Line_z2R(const Point& ob, const Point& v0,const Point& v1)
{
	double result = 0.0;
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	if (abs(s0 + R0) < TOLERANCE) {
		result = (t.z)*(t.z)*0.25 * (s1*R1*R1*R1 - s0*R0*R0*R0)
			+ 2.0*(L.z)*(t.z)*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.z)*(L.z)*(s1*R1 - s0*R0);
	}
	else {
		result = (t.z)*(t.z)*0.125*(2 * (s1*R1*R1*R1 - s0*R0*R0*R0) - L0*L0*(s1*R1 - s0*R0) - L0*L0*L0*L0*log((s1 + R1) / (s0 + R0)))
			+ 2.0*(L.z)*(t.z)*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.z)*(L.z)*(s1*R1 - s0*R0 + L0*L0*log((s1 + R1) / (s0 + R0)));
	}
	return result;
}

double Integral::Line_xzR(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	if (abs(s0 + R0) < TOLERANCE) {
		result = (t.x)*(t.z)*0.25 * (s1*R1*R1*R1 - s0*R0*R0*R0)
			+ ((L.x)*(t.z)+(L.z)*(t.x))*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.x)*(L.z)*(s1*R1 - s0*R0);
	}
	else {
		result = (t.x)*(t.z)*0.125*(2 * (s1*R1*R1*R1 - s0*R0*R0*R0) - L0*L0*(s1*R1 - s0*R0) - L0*L0*L0*L0*log((s1 + R1) / (s0 + R0)))
			+ ((L.x)*(t.z)+ (L.z)*(t.x))*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.x)*(L.z)*(s1*R1 - s0*R0 + L0*L0*log((s1 + R1) / (s0 + R0)));
	}
	return result;
}

double Integral::Line_yzR(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	if (abs(s0 + R0) < TOLERANCE) {
		result = (t.y)*(t.z)*0.25 * (s1*R1*R1*R1 - s0*R0*R0*R0)
			+ ((L.y)*(t.z) + (L.z)*(t.y))*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.y)*(L.z)*(s1*R1 - s0*R0);
	}
	else {
		result = (t.y)*(t.z)*0.125*(2 * (s1*R1*R1*R1 - s0*R0*R0*R0) - L0*L0*(s1*R1 - s0*R0) - L0*L0*L0*L0*log((s1 + R1) / (s0 + R0)))
			+ ((L.y)*(t.z) + (L.z)*(t.y))*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.y)*(L.z)*(s1*R1 - s0*R0 + L0*L0*log((s1 + R1) / (s0 + R0)));
	}
	return result;
}

double Integral::Line_xyR(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point t = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*t;
	s1 = (v1 - ob)*t;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	if (abs(s0 + R0) < TOLERANCE) {
		result = (t.x)*(t.y)*0.25 * (s1*R1*R1*R1 - s0*R0*R0*R0)
			+ ((L.x)*(t.y) + (L.y)*(t.x))*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.x)*(L.y)*(s1*R1 - s0*R0);
	}
	else {
		result = (t.x)*(t.y)*0.125*(2 * (s1*R1*R1*R1 - s0*R0*R0*R0) - L0*L0*(s1*R1 - s0*R0) - L0*L0*L0*L0*log((s1 + R1) / (s0 + R0)))
			+ ((L.x)*(t.y) + (L.y)*(t.x))*(R1*R1*R1 - R0*R0*R0) / 3.0
			+ 0.5*(L.x)*(L.y)*(s1*R1 - s0*R0 + L0*L0*log((s1 + R1) / (s0 + R0)));
	}
	return result;
}

double Integral::Line_x2R_1(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point e = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*e;
	s1 = (v1 - ob)*e;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	double Lx = L.x;
	double ex = e.x;
	if (abs(s0 + R0) < TOLERANCE) {
		assert(s0*s1>0);
		result = ex*ex*0.5*(s1*R1 - s0*R0) + 2 * Lx*ex*(R1 - R0) + Lx*Lx*abs(log(s1 / s0));
	}
	else {
		result = ex*ex*0.5*((s1*R1 - s0*R0) - L0*L0*log((s1 + R1) / (s0 + R0))) + 2 * Lx*ex*(R1 - R0) + Lx*Lx*log((s1 + R1) / (s0 + R0));
	}
	return result;
}

double Integral::Line_y2R_1(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point e = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*e;
	s1 = (v1 - ob)*e;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	double Ly = L.y;
	double ey = e.y;
	if (abs(s0 + R0) < TOLERANCE) {
		assert(s0*s1>0);
		result = ey*ey*0.5*(s1*R1 - s0*R0) + 2 * Ly*ey*(R1 - R0) + Ly*Ly*abs(log(s1 / s0));
	}
	else {
		result = ey*ey*0.5*((s1*R1 - s0*R0) - L0*L0*log((s1 + R1) / (s0 + R0))) + 2 * Ly*ey*(R1 - R0) + Ly*Ly*log((s1 + R1) / (s0 + R0));
	}
	return result;
}

double Integral::Line_z2R_1(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point e = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*e;
	s1 = (v1 - ob)*e;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	double Lz = L.z;
	double ez = e.z;
	if (abs(s0 + R0) < TOLERANCE) {
		assert(s0*s1>0);
		result = ez*ez*0.5*(s1*R1 - s0*R0) + 2 * Lz*ez*(R1 - R0) + Lz*Lz*abs(log(s1 / s0));
	}
	else {
		result = ez*ez*0.5*((s1*R1 - s0*R0) - L0*L0*log((s1 + R1) / (s0 + R0))) + 2 * Lz*ez*(R1 - R0) + Lz*Lz*log((s1 + R1) / (s0 + R0));
	}
	return result;
}

double Integral::Line_xzR_1(const Point& ob, const Point& v0,const Point& v1)
{
	double result = 0.0;
	Point e = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*e;
	s1 = (v1 - ob)*e;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	double Lx = L.x,Lz=L.z;
	double ex = e.x,ez=e.z;
	if (abs(s0 + R0) < TOLERANCE) {
		assert(s0*s1>0);
		result = ex*ez*0.5*(s1*R1 - s0*R0) + (Lx*ez + Lz*ex)*(R1 - R0) + Lx*Lz*abs(log(s1 / s0));
	}
	else {
		result = ex*ez*0.5*((s1*R1 - s0*R0) - L0*L0*log((s1 + R1) / (s0 + R0))) + (Lx*ez + Lz*ex)*(R1 - R0) + Lx*Lz*log((s1 + R1) / (s0 + R0));
	}
	return result;
}
double Integral::Line_yzR_1(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point e = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*e;
	s1 = (v1 - ob)*e;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	double Ly = L.y, Lz = L.z;
	double ey = e.y, ez = e.z;
	if (abs(s0 + R0) < TOLERANCE) {
		assert(s0*s1>0);
		result = ey*ez*0.5*(s1*R1 - s0*R0) + (Ly*ez + Lz*ey)*(R1 - R0) + Ly*Lz*abs(log(s1 / s0));
	}
	else {
		result = ey*ez*0.5*((s1*R1 - s0*R0) - L0*L0*log((s1 + R1) / (s0 + R0))) + (Ly*ez + Lz*ey)*(R1 - R0) + Ly*Lz*log((s1 + R1) / (s0 + R0));
	}
	return result;
}
double Integral::Line_xyR_1(const Point& ob, const Point& v0, const Point& v1)
{
	double result = 0.0;
	Point e = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;
	s0 = (v0 - ob)*e;
	s1 = (v1 - ob)*e;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);
	double Ly = L.y, Lx = L.x;
	double ey = e.y, ex = e.x;
	if (abs(s0 + R0) < TOLERANCE) {
		assert(s0*s1>0);
		result = ey*ex*0.5*(s1*R1 - s0*R0) + (Ly*ex + Lx*ey)*(R1 - R0) + Ly*Lx*abs(log(s1 / s0));
	}
	else {
		result = ey*ex*0.5*((s1*R1 - s0*R0) - L0*L0*log((s1 + R1) / (s0 + R0))) + (Ly*ex + Lx*ey)*(R1 - R0) + Ly*Lx*log((s1 + R1) / (s0 + R0));
	}
	return result;
}
double Integral::Line_f3R_1(const Point& ob, const Point& v0, const Point& v1,string s)
{
	Point e = (v1 - v0).getUnit();
	Point prj = ob.projection_line(v0, v1);
	Point L = prj - ob;
	double L0 = L.size();
	double s0, s1, R0, R1;

	s0 = (v0 - ob)*e;
	s1 = (v1 - ob)*e;
	R0 = v0.distance(ob);
	R1 = v1.distance(ob);


	//double Lx = L.x, Ly = L.y, Lz = L.z;
	//double ex = e.x, ey = e.y, ez = e.z;
	double Li[3] = { L.x,L.y,L.z };
	double ei[3] = { e.x,e.y,e.z };
	double itg0 = 0.0, itg1 = 0.0, itg2 = 0.0,itg3=0.0;
	double result = 0.0;

	if (abs(s0 + R0) < TOLERANCE) {
		assert(s0*s1>0);
		itg0 = abs(log(s1 / s0));
		itg2 = 0.5*(s1*R1 - s0*R0);
	}
	else {
		itg0 = log((s1 + R1) / (s0 + R0));
		itg2 = 0.5*((s1*R1 - s0*R0) - L0*L0*log((s1 + R1) / (s0 + R0)));
	}
	itg1 = R1 - R0;
	itg3 = ((s1*s1 - 2.0*L0*L0)*R1 - (s0*s0 - 2.0*L0*L0)*R0) / 3.0;
	int p = 0, q = 0, t = 0;
	if (s == "x3") {}
	else if (s == "y3") {
		p = 1; q = 1; t = 1;
	}
	else if (s == "z3") {
		p = 2; q = 2; t = 2;
	}
	else if ((s == "xy2")||(s=="y2x")) {
		p = 0; q = 1; t = 1;
	}
	else if ((s == "x2y") || (s == "yx2")) {
		p = 0; q = 0; t = 1;
	}
	else if ((s == "xz2") || (s == "z2x")) {
		p = 0; q = 2; t = 2;
	}
	else if ((s == "x2z") ||(s == "zx2")) {
		p = 0; q = 0; t = 2;
	}
	else if ((s == "yz2") || (s == "z2y")) {
		p = 1; q = 2; t = 2;
	}
	else if ((s == "y2z" || (s == "zy2"))) {
		p = 1; q = 1; t = 2;
	}
	else if ((s == "xyz")) {
		p = 0; q = 1; t = 2;
	}
	else {
		cout << "error in argument s in Integral::Integral::Line_f3R_1(Point ob, Point v0, Point v1,string s)\n";
		std::abort();
	}
	result = Li[p] * Li[q] * Li[t] * itg0
		+ (Li[p] * Li[q] * ei[t] + Li[p] * Li[t] * ei[q] + Li[q] * Li[t] * ei[p])*itg1
		+ (Li[t] * ei[p] * ei[q] + Li[p] * ei[q] * ei[t] + Li[q] * ei[p] * ei[t])*itg2
		+ ei[p] * ei[q] * ei[t] * itg3;
	return result;
}
double Integral::getBeta(const Point& observation, const vector<Node*>& face)
{
	int n = face.size();//node number of the face
	assert(n >= 3);
	Point mj,t;
	Point unit_normal = unitCross(*face[1] - *face[0], *face[2] - *face[1]);
	//unit_normal.display();
	double beta = 0.0;
	double s1, s0, mj_value;
	//face.push_back(face[0]);
	int i_next = 0;
	for (int i = 0; i < n; i++) {
		i_next = (i + 1) % n;
		t=(*face[i_next] - *face[i]).getUnit();
		s0 = (*face[i] - observation)*t;
		s1 = (*face[i_next] - observation)*t;
		mj = cross(t, unit_normal);
		mj_value = mj * (*face[i] - observation);
		double mp = mj *(*face[i_next] - observation);
		double tt = 0.;
		if (abs(mp)< TOLERANCE) { tt = 0.; }
		else if (mp>0) { tt = 1.; }
		else { tt = -1.; }
		if (abs(mj_value) < TOLERANCE) { beta += 0.0; }
		else {
		//	cout <<(i+1)<<'\t'<< mj_value << endl;
			beta =beta+ tt*(atan(s1 / abs(mj_value)) - atan(s0 / abs(mj_value)));
		}
	}
	return beta;
}



double Integral::Surface_R1(const Point& observation, const vector<Node*>& face_node, int n)
{
	//assuming that face nodes have been arranged correctly
	int nn = face_node.size();//node number of the face
	assert(n ==nn);	
	double beta = getBeta(observation, face_node);//solid angle
	double result = 0.0;
	//unit normal vector
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (observation - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	for (int i = 0; i < n; i++) {
		result= result+1./3.*Line_Bj3(observation, *face_node[i], *face_node[(i + 1)%n], unit_normal);
	}

	result = result - beta/3.0*abs(h)*abs(h)*abs(h);
	return result;
}

double Integral::Surface_R_1(const Point& observation, const vector<Node*>& face_node,int n)
{
	int nn = face_node.size();//node number of the face
	assert(n == nn);
	
	double beta = getBeta(observation, face_node);//solid angle
	//unit normal vector
	Point unit_normal= unitCross(*face_node[1] - *face_node[0], *face_node[n-1] - *face_node[0]);
	double h = (observation - *face_node[0])*unit_normal;
	//cout<<h<<endl;
	//face_node.push_back(face_node[0]);
	double result = 0.0;
	unsigned int next = 0;
	for (int j = 0; j < n; j++) {
		next = (j + 1) % n;
		Point t = (*face_node[next] - *face_node[j]).getUnit();
		Point mj = unitCross(t, unit_normal);
		double mj_value = (*face_node[j] - observation)*mj;
		double bj1=0.0;
		if (abs(mj_value) < TOLERANCE) {
			bj1 = 0.0;
		}
		else {
			double s0 = t*(*face_node[j] - observation);
			double s1 = t*(*face_node[next] - observation);
			double R0 = observation.distance(*face_node[j]);
			double R1 = observation.distance(*face_node[next]);
			bj1 = abs(h)*(atan((abs(h)*s1) / (mj_value*R1)) - atan((abs(h)*s0) / (mj_value*R0)))
				+ mj_value*log((s1 + R1) / (s0 + R0));
		}
		
		result += bj1;
		//result += Line_Bj1(observation, face_node[i], face_node[i + 1], unit_normal);
	}
	result = result - beta*abs(h);
	//cout << result << endl;
	return result;
}

Point Integral::Surface_rR1(const Point& observation, const vector<Node*>& face_node, int n)
{
	Point i_rR1(0.0,0.0,0.0);
	Point prj = observation.projection(*face_node[0], *face_node[1], *face_node[2]);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	Point mj;
	double i_R1 = Surface_R1(observation, face_node, n);
	i_rR1 += (prj - observation)*i_R1;

	//face_node.push_back(face_node[0]);
	unsigned int next = 0;
	for (int i = 0; i < n ; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next]-*face_node[i], unit_normal);
		i_rR1 +=1./3.* mj*Line_R3(observation, *face_node[i], *face_node[next]);
	}
	return i_rR1;
}

Point Integral::Surface_rR_1(const Point& observation, const vector<Node*>& face_node, int n)
{
	Point rR_1(0.0, 0.0, 0.0);
	Point prj = observation.projection(*face_node[0], *face_node[1], *face_node[2]);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[0]);
	Point mj;
	double R_1 = Surface_R_1(observation, face_node, n);
	rR_1 = rR_1 + (prj - observation)*R_1;
	//face_node.push_back(face_node[0]);
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next] - *face_node[i], unit_normal);
		rR_1 = rR_1 + mj*Line_R1(observation, *face_node[i], *face_node[next]);
	}	
	return rR_1;
	
}

Point Integral::Surface_r2R_1(const Point& observation, const vector<Node*>& face_node, int n)
{
	double Rds = Surface_R1(observation, face_node, n);
	Point rR_1ds = Surface_rR_1(observation, face_node, n);

	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[0]);
	
	Point rRdl(0,0,0),mj;
	Point r2R_1(0,0,0);
	for (int i = 0; i < n-1; i++) {
		rRdl = Line_rR1(observation, *face_node[i], *face_node[i+1]);
		//rRdl.display();
		mj = unitCross(*face_node[i + 1] - *face_node[i],unit_normal);
		r2R_1.x = r2R_1.x + mj.x*rRdl.x;
		r2R_1.y = r2R_1.y + mj.y*rRdl.y;
		r2R_1.z = r2R_1.z + mj.z*rRdl.z;
	}
	//the last edge
	rRdl = Line_rR1(observation, *face_node[n-1], *face_node[0]);
	mj = unitCross(*face_node[0] - *face_node[n-1], unit_normal);
	r2R_1.x = r2R_1.x + mj.x*rRdl.x;
	r2R_1.y = r2R_1.y + mj.y*rRdl.y;
	r2R_1.z = r2R_1.z + mj.z*rRdl.z;
	double h = (observation - *face_node[0])*unit_normal;
	r2R_1.x = r2R_1.x - unit_normal.x*h*(rR_1ds.x) - (1 - (unit_normal.x)*(unit_normal.x))*Rds;
	r2R_1.y = r2R_1.y - unit_normal.y*h*(rR_1ds.y) - (1 - (unit_normal.y)*(unit_normal.y))*Rds;
	r2R_1.z = r2R_1.z - unit_normal.z*h*(rR_1ds.z) - (1 - (unit_normal.z)*(unit_normal.z))*Rds;
	
	return r2R_1;
}

double Integral::Surface_xzR_1(const Point& observation, const vector<Node*>& face_node, int n)
{
	double Rds = Surface_R1(observation, face_node, n);
	Point rR_1ds = Surface_rR_1(observation, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (observation - *face_node[0])*unit_normal;
	Point rRdl, mj;
	double xzR_1=0.0;
	for (int i = 0; i < n - 1; i++) {
		rRdl = Line_rR1(observation, *face_node[i], *face_node[i + 1]);
		mj = unitCross(*face_node[i+1]-*face_node[i],unit_normal);
		xzR_1 += mj.z*rRdl.x;
	}
	rRdl = Line_rR1(observation, *face_node[n-1], *face_node[0]);
	mj = unitCross(*face_node[0] - *face_node[n-1], unit_normal);
	xzR_1 += mj.z*rRdl.x;
	xzR_1 = xzR_1 - unit_normal.z*h*rR_1ds.x + unit_normal.z*unit_normal.x*Rds;
	return xzR_1;
}

double Integral::Surface_yzR_1(const Point& observation, const vector<Node*>& face_node, int n)
{
	double Rds = Surface_R1(observation, face_node, n);
	Point rR_1ds = Surface_rR_1(observation, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (observation - *face_node[0])*unit_normal;
	Point rRdl, mj;
	double yzR_1=0.0;
	for (int i = 0; i < n - 1; i++) {
		rRdl = Line_rR1(observation, *face_node[i], *face_node[i + 1]);
		mj = unitCross(*face_node[i + 1] - *face_node[i], unit_normal);
		yzR_1 += mj.z*rRdl.y;
	}
	rRdl = Line_rR1(observation, *face_node[n-1], *face_node[0]);
	mj = unitCross(*face_node[0] - *face_node[n-1], unit_normal);
	yzR_1 += mj.z*rRdl.y;
	yzR_1 = yzR_1 - unit_normal.z*h*rR_1ds.y + unit_normal.z*unit_normal.y*Rds;
	return yzR_1;
}

double Integral::Surface_xyR_1(const Point& observation, const vector<Node*> &face_node, int n)
{
	double Rds = Surface_R1(observation, face_node, n);
	Point rR_1ds = Surface_rR_1(observation, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (observation - *face_node[0])*unit_normal;
	Point rRdl, mj;
	double xyR_1 = 0.0;
	for (int i = 0; i < n - 1; i++) {
		rRdl = Line_rR1(observation, *face_node[i], *face_node[i + 1]);
		mj = unitCross(*face_node[i + 1] - *face_node[i], unit_normal);
		xyR_1 += mj.y*rRdl.x;
	}
	rRdl = Line_rR1(observation, *face_node[n-1], *face_node[0]);
	mj = unitCross(*face_node[0] - *face_node[n-1], unit_normal);
	xyR_1 += mj.y*rRdl.x;
	xyR_1 = xyR_1 - unit_normal.y*h*rR_1ds.x + unit_normal.x*unit_normal.y*Rds;
	return xyR_1;
}

double Integral::Surface_R_3(const Point& ob, const vector<Node*>& face_node, int n)
{
	//face_node.push_back(face_node[0]);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	assert(!(abs(h) < TOLERANCE));
	double R1 = 0.0, R0 = 0.0, s1 = 0.0, s0 = 0.0,m0=0.0;
	Point e,mj;
	double result = 0.0;
	//Point prj = ob.projection_line(v0, v1);
	//Point L = prj - ob;
	double L0 = 0.0;
	Point prj;
	unsigned next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		e = (*face_node[next] - *face_node[i]).getUnit();
		mj= unitCross(*face_node[next] - *face_node[i], unit_normal);
		prj = ob.projection_line(*face_node[next], *face_node[i]);
		L0 = (prj - ob).size();
		s1 = (*face_node[next] - ob)*e;
		s0 = (*face_node[i] - ob)*e;
		R1 = ob.distance(*face_node[next]);
		R0 = ob.distance(*face_node[i]);
		m0 = (*face_node[i] - ob)*mj;
		result += atan((m0*s1)/(L0*L0+abs(h)*R1)) - atan((m0*s0)/(L0*L0+abs(h)*R0));
	}
	result = result / abs(h);
	return result;
}

Point Integral::Surface_rR_3(const Point& ob, const vector<Node*>& face_node, int n)
{
		
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	Point result(0, 0, 0);
	Point mj;
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj=unitCross(*face_node[next] - *face_node[i], unit_normal);
		result = result - mj*Line_R_1(ob, *face_node[i], *face_node[next]);
	}
	if (abs(h) < TOLERANCE) {}
	else{
		double sR_3 = Surface_R_3(ob, face_node, n);
		result = result - unit_normal*h*sR_3;
	}
	//if (abs(h) < TOLERANCE) {
	//	cout << endl;
	//	result.display();
	//	for (int i = 0; i < face_node.size(); i++) {
	//		(*face_node[i]).display();
	//	}
	//	cout << endl;
	//}
	return result;
}

double Integral::Surface_x2R_3(const Point& ob, const vector<Node*>& face_node, int n)
{
	Point S_rR_3 = Surface_rR_3(ob, face_node, n);
	double S_R_1 = Surface_R_1(ob, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	Point mj(0,0,0);
	double result = 0.0;
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next] - *face_node[i], unit_normal);
		result = result + (mj.x)*(Line_rR_1(ob, *face_node[i], *face_node[next]).x);
	}
	result = result + (unit_normal.x)*h*(S_rR_3.x) - (1 - (unit_normal.x)*(unit_normal.x))*S_R_1;
	result = -result;
	return result;
}

double Integral::Surface_y2R_3(const Point& ob, const vector<Node*>& face_node, int n)
{
	Point S_rR_3 = Surface_rR_3(ob, face_node, n);
	double S_R_1 = Surface_R_1(ob, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	Point mj(0, 0, 0);
	double result = 0.0;
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next] - *face_node[i], unit_normal);
		result = result + (mj.y)*(Line_rR_1(ob, *face_node[i], *face_node[next]).y);
	}
	result = result + (unit_normal.y)*h*(S_rR_3.y) - (1 - (unit_normal.y)*(unit_normal.y))*S_R_1;
	result = -result;
	return result;
}

double Integral::Surface_z2R_3(const Point& ob, const vector<Node*>& face_node, int n)
{
	Point S_rR_3 = Surface_rR_3(ob, face_node, n);
	double S_R_1 = Surface_R_1(ob, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	Point mj(0, 0, 0);
	double result = 0.0;
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next] - *face_node[i], unit_normal);
		result = result + (mj.z)*(Line_rR_1(ob, *face_node[i], *face_node[next]).z);
	}
	result = result + (unit_normal.z)*h*(S_rR_3.z) - (1 - (unit_normal.z)*(unit_normal.z))*S_R_1;
	result = -result;
	return result;
}

double Integral::Surface_xzR_3(const Point& ob, const vector<Node*>& face_node, int n)
{
	Point S_rR_3 = Surface_rR_3(ob, face_node, n);
	double S_R_1 = Surface_R_1(ob, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	Point mj(0, 0, 0);
	double result = 0.0;
	unsigned next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next] - *face_node[i], unit_normal);
		result = result + (mj.z)*(Line_rR_1(ob, *face_node[i], *face_node[next]).x);
	}
	result = result + (unit_normal.z)*h*(S_rR_3.x) + (unit_normal.x)*(unit_normal.z)*S_R_1;
	result = -result;
	return result;
}

double Integral::Surface_yzR_3(const Point& ob, const vector<Node*>& face_node, int n)
{
	Point S_rR_3 = Surface_rR_3(ob, face_node, n);
	double S_R_1 = Surface_R_1(ob, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	Point mj(0, 0, 0);
	double result = 0.0;
	unsigned next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next] - *face_node[i], unit_normal);
		result = result + (mj.z)*(Line_rR_1(ob, *face_node[i], *face_node[next]).y);
	}
	result = result + (unit_normal.z)*h*(S_rR_3.y) + (unit_normal.y)*(unit_normal.z)*S_R_1;
	result = -result;
	return result;
}

double Integral::Surface_xyR_3(const Point& ob, const vector<Node*>& face_node, int n)
{
	Point S_rR_3 = Surface_rR_3(ob, face_node, n);
	double S_R_1 = Surface_R_1(ob, face_node, n);
	Point unit_normal = unitCross(*face_node[1] - *face_node[0], *face_node[2] - *face_node[1]);
	double h = (ob - *face_node[0])*unit_normal;
	//face_node.push_back(face_node[0]);
	Point mj(0, 0, 0);
	double result = 0.0;
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mj = unitCross(*face_node[next] - *face_node[i], unit_normal);
		result = result + (mj.x)*(Line_rR_1(ob, *face_node[i], *face_node[next]).y);
	}
	result = result + (unit_normal.x)*h*(S_rR_3.y) + (unit_normal.y)*(unit_normal.x)*S_R_1;
	result = -result;
	return result;
}

void Integral::Surface_r3R_1(const Point& ob, const vector<Node*>& face, int n,
	double & x3R_1, double & y3R_1, double & z3R_1, 
	double & x2yR_1, double & x2zR_1, double & y2xR_1, double & y2zR_1, double & z2xR_1, double & z2yR_1,
	double & xyzR_1)
{
	int nn = face.size();
	assert(n == nn);
	Point r2R_1 = Surface_r2R_1(ob, face, n);
	double xyR_1= Surface_xyR_1(ob, face, n);
	//double xzR_1 = Surface_xzR_1(ob, face, n);
	//double yzR_1 = Surface_yzR_1(ob, face, n);
	Point rR = Surface_rR1(ob, face, n);
	vector<double> x2Rdl(n);
	vector<double> y2Rdl(n);
	vector<double> z2Rdl(n);
	vector<double> xyRdl(n);
	vector<Point> mij;
	Point unit_normal = unitCross(*face[1] - *face[0], *face[2] - *face[0]);
	//cout << "BB" << face.size() << endl;
	//face.push_back(face[0]);
	//cout << "BB"<<face.size() << endl;
	mij.resize(n);
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mij[i] = unitCross(*face[next] - *face[i],unit_normal);
		x2Rdl[i] = Line_x2R(ob, *face[i], *face[next]);
		y2Rdl[i] = Line_y2R(ob, *face[i], *face[next]);
		z2Rdl[i] = Line_z2R(ob, *face[i], *face[next]);
		xyRdl[i] = Line_xyR(ob, *face[i], *face[next]);
	}
	double h = (ob - *face[0])*unit_normal;
	
	x3R_1 = -unit_normal.x*h*r2R_1.x - 2 * (1 - (unit_normal.x)*(unit_normal.x))*rR.x;
	y3R_1 = -unit_normal.y*h*r2R_1.y - 2 * (1 - (unit_normal.y)*(unit_normal.y))*rR.y;
	z3R_1 = -unit_normal.z*h*r2R_1.z - 2 * (1 - (unit_normal.z)*(unit_normal.z))*rR.z;	
	x2yR_1 = -unit_normal.y*h*r2R_1.x + 2 * (unit_normal.y)*(unit_normal.x)*rR.x;
	x2zR_1 = -unit_normal.z*h*r2R_1.x + 2 * (unit_normal.z)*(unit_normal.x)*rR.x;
	y2xR_1 = -unit_normal.x*h*r2R_1.y + 2 * (unit_normal.x)*(unit_normal.y)*rR.y;
	y2zR_1 = -unit_normal.z*h*r2R_1.y + 2 * (unit_normal.z)*(unit_normal.y)*rR.y;
	z2xR_1 = -unit_normal.x*h*r2R_1.z + 2 * (unit_normal.x)*(unit_normal.z)*rR.z;
	z2yR_1 = -unit_normal.y*h*r2R_1.z + 2 * (unit_normal.y)*(unit_normal.z)*rR.z;
	xyzR_1 = -unit_normal.z*h*xyR_1
		+unit_normal.z*unit_normal.y*rR.x+unit_normal.z*unit_normal.x*rR.y;
		
	for (int i = 0; i < n; i++) {
		x3R_1 = x3R_1 + mij[i].x*x2Rdl[i];
		y3R_1 = y3R_1 + mij[i].y*y2Rdl[i];
		z3R_1 = z3R_1 + mij[i].z*z2Rdl[i];
		x2yR_1 = x2yR_1 + mij[i].y*x2Rdl[i];
		x2zR_1 = x2zR_1 + mij[i].z*x2Rdl[i];
		y2xR_1 = y2xR_1 + mij[i].x*y2Rdl[i];
		y2zR_1 = y2zR_1 + mij[i].z*y2Rdl[i];
		z2xR_1 = z2xR_1 + mij[i].x*z2Rdl[i];
		z2yR_1 = z2yR_1 + mij[i].y*z2Rdl[i];
		xyzR_1 = xyzR_1 + mij[i].z*xyRdl[i];
	}
	return ;
}

void Integral::Surface_f3R_3(const Point& ob, const vector<Node*>& face, int n,
	double & x3R_3, double & y3R_3, double & z3R_3,
	double & x2yR_3, double & x2zR_3,
	double & y2xR_3, double & y2zR_3,
	double & z2xR_3, double & z2yR_3,
	double & xyzR_3)
{
	int nn = face.size();
	assert(n == nn);

	double x2R_3 = Surface_x2R_3(ob, face, n);
	double y2R_3 = Surface_y2R_3(ob, face, n);
	double z2R_3 = Surface_z2R_3(ob, face, n);
	double xyR_3 = Surface_xyR_3(ob, face, n);

	Point rR_1 = Surface_rR_1(ob, face, n);

	vector<double> x2R_1(n);
	vector<double> y2R_1(n);
	vector<double> z2R_1(n);
	vector<double> xyR_1(n);
	vector<Point> mij(n);
	//face.push_back(face[0]);
	Point unit_normal = unitCross(*face[1] - *face[0], *face[2] - *face[0]);
	unsigned int next = 0;
	for (int i = 0; i < n; i++) {
		next = (i + 1) % n;
		mij[i] = unitCross(*face[next] - *face[i], unit_normal);
		x2R_1[i] = Line_x2R_1(ob, *face[i], *face[next]);
		y2R_1[i] = Line_y2R_1(ob, *face[i], *face[next]);
		z2R_1[i] = Line_z2R_1(ob, *face[i], *face[next]);
		xyR_1[i] = Line_xyR_1(ob, *face[i], *face[next]);
	}
	double h = (ob - *face[0])*unit_normal;

	x3R_3 = unit_normal.x*h*x2R_3 - 2 * (1 - (unit_normal.x)*(unit_normal.x))*rR_1.x;
	y3R_3 = unit_normal.y*h*y2R_3 - 2 * (1 - (unit_normal.y)*(unit_normal.y))*rR_1.y;
	z3R_3 = unit_normal.z*h*z2R_3 - 2 * (1 - (unit_normal.z)*(unit_normal.z))*rR_1.z;
	x2yR_3 = unit_normal.y*h*x2R_3 + 2 * (unit_normal.y)*(unit_normal.x)*rR_1.x;
	x2zR_3 = unit_normal.z*h*x2R_3 + 2 * (unit_normal.z)*(unit_normal.x)*rR_1.x;
	y2xR_3 = unit_normal.x*h*y2R_3 +	2 * (unit_normal.x)*(unit_normal.y)*rR_1.y;
	y2zR_3 = unit_normal.z*h*y2R_3 + 2 * (unit_normal.z)*(unit_normal.y)*rR_1.y;
	z2xR_3 = unit_normal.x*h*z2R_3 + 2 * (unit_normal.x)*(unit_normal.z)*rR_1.z;
	z2yR_3 = unit_normal.y*h*z2R_3 + 2 * (unit_normal.y)*(unit_normal.z)*rR_1.z;
	xyzR_3 = unit_normal.z*h*xyR_3
		+unit_normal.z*unit_normal.y*rR_1.x+unit_normal.z*unit_normal.x*rR_1.y;

	for (int i = 0; i < n; i++) {
		x3R_3 = x3R_3 + mij[i].x*x2R_1[i];
		y3R_3 = y3R_3 + mij[i].y*y2R_1[i];
		z3R_3 = z3R_3 + mij[i].z*z2R_1[i];
		x2yR_3 = x2yR_3 + mij[i].y*x2R_1[i];
		x2zR_3 = x2zR_3 + mij[i].z*x2R_1[i];
		y2xR_3 = y2xR_3 + mij[i].x*y2R_1[i];
		y2zR_3 = y2zR_3 + mij[i].z*y2R_1[i];
		z2xR_3 = z2xR_3 + mij[i].x*z2R_1[i];
		z2yR_3 = z2yR_3 + mij[i].y*z2R_1[i];
		xyzR_3 = xyzR_3 + mij[i].z*xyR_1[i];
	}
	
	x3R_3 = -x3R_3;
	y3R_3 = -y3R_3;
	z3R_3 = -z3R_3;
	x2yR_3 = -x2yR_3;
	x2zR_3 = -x2zR_3;
	y2xR_3 = -y2xR_3;
	y2zR_3 = -y2zR_3;
	z2xR_3 = -z2xR_3;
	z2yR_3 = -z2yR_3;
	xyzR_3 = -xyzR_3;
	return;
}



