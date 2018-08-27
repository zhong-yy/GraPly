#pragma once
/*****************************************************************
*Gravity class is to calculate gravity and gravity gradient of a 
*single polyhedral mass body with (cubic) polynomial density contrast
*
*Copyright 2017
*Zhong, Yiyuan 
*Central South University
******************************************************************/
#include"Integral.h"
#include"Polyhedral.h"
#include"Dyadic.h"
class Gravity
{
public:
	Gravity();
	~Gravity();
	
	/******Computation of gravity field due to a polyhedral body with polynomial density contrast*******/
	Point g_const(const Polyhedral& phl, const Point& ob, const double& rho);
	//constant density, x-component of gravity field, unit:mGal
	double gx_const(const Polyhedral& phl, const Point& ob, const double& rho);
	//constant density, y-component of gravity field, unit:mGal
	double gy_const(const Polyhedral& phl, const Point& ob, const double& rho);
	//constant density, x-component of gravity field, unit:mGal
	double gz_const(const Polyhedral& phl, const Point& ob, const double& rho);
	
	//linearly varying density in x,y,z direction, unit:mGal
	Point g_1st(const Polyhedral& phl, const Point& ob, 
		const double& a, const double& b, const double& c);

	double gx_1st(const Polyhedral& phl, const Point& ob, 
		const double& a, const double& b, const double& c);

	double gy_1st(const Polyhedral& phl, const Point& ob, 
		const double& a, const double& b, const double& c);
	double gz_1st(const Polyhedral& phl, const Point& ob, 
		const double& a, const double& b, const double& c);
	
	

	//quadraticly varying density in x y z direction, unit:mGal
	//rho:a*x^2+by^2+c*z^2+d*xz+e*yz+f*xy
	Point g_2nd(const Polyhedral& poh, const Point& ob, 
		const double& a, const double& b, const double& c, 
		const double& d, const double& e, const double& f);
	double gx_2nd(const Polyhedral& poh, const Point& ob, 
		const double& a, const double& b, const double& c, const double& d, const double& e, const double& f);
	double gy_2nd(const Polyhedral& poh, const Point& ob, 
		const double& a, const double& b, const double& c, const double& d, const double& e, const double& f);
	double gz_2nd(const Polyhedral& poh, const Point& ob, 
		double a, double b, double c, double d, double e, double f);

	Point g_3rd(const Polyhedral& po, const Point& ob, double a[10]);


	/************computation of gravity gradient tensor**************************************************/

	//constant density contrast
	Dyadic tensor_const(const Polyhedral& phl, const Point& ob, const double& rho);

	//linear density contrast
	Dyadic tensor_1st(const Polyhedral& phl, const Point& ob, const double& a, const double &b, const double& c);

	//quadratic density contrast
	Dyadic tensor_2nd(const Polyhedral& phl, const Point& ob, double a, double b, double c, double d, double e, double f);

	//cubic density contrast
	Dyadic tensor_3rd(const Polyhedral& po, const Point& ob, double a[10]);



private:
	//local coordinate system where ob is the origin
	//but the parameter "ob" is described by global coordinates
	Point g_1st0(const Polyhedral& phl, const Point& ob, double a, double b, double c);
	double gx_1st0(const Polyhedral& phl, const Point& ob, double a, double b, double c);
	double gy_1st0(const Polyhedral& phl, const Point& ob, double a, double b, double c);
	double gz_1st0(const Polyhedral& phl, const Point& ob,double a,double b,double c);

	Point g_2nd0(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f);
	double gx_2nd0(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f);
	double gy_2nd0(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f);
	double gz_2nd0(const Polyhedral& poh, const Point& ob, double a,double b,double c,double d,double e,double f);

	Point g_3rd0(const Polyhedral& po, const Point& ob,double a[10]);

	//in local coordinate system
	Dyadic tensor_1st0(const Polyhedral& phl, const Point& ob, double a, double b, double c);
	Dyadic tensor_2nd0(const Polyhedral& phl, const Point& ob, double a, double b, double c, double d, double e, double f);
	Dyadic tensor_3rd0(const Polyhedral& po, const Point& ob, double a[10]);


};
