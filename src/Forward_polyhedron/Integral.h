#pragma once

/*****************************************************************
 * Integral class is to calculte some useful line integrals and surface
 * integrals in terms of closed-form formulae.
 *
 * Copyright 2017
 * Zhong, Yiyuan 
 * Central South University
 ******************************************************************/
#include<vector>
#include<cassert>
#include <cstdlib>
#include"Point.h"
#include"Node.h"
#include<string>

using namespace std;
class Integral
{
public:
	Integral();
	~Integral();


	//caculate the solid angle of projection onto face from observation site
	static double getBeta(const Point& observation, const vector<Node*>& face);



	/****************************************Line integrals*****************************************/
	//line integral of R, v0 and v1 are the endpoints of the edge
	static double Line_R1(const Point& observation, const Point& v0, const Point& v1);

	//line integral of R^3
	static double Line_R3(const Point& observation, const Point& v0, const Point& v1);

	/*calculate the line integral of 1/R,observation site cannot locate on edge v0-v1,
	but it can locate on the extension of the edge*/
	static double Line_R_1(const Point& observation, const Point& v0, const Point& v1);

	//Bj:mj*integral(R/(mj+s)^2)dl,see equation(51)-(52) Ren(2017) in "Survey in Geophysics"
	//normal is the unit vector normal to the plane, which is used to determine the plane
	static double Line_Bj1(const Point& observation,const Point& v0, const Point& v1, const Point& normal);
	static double Line_Bj3(const Point& observation,const Point& v0, const Point& v1, const Point& unit_normal);
	
	//calculate the line integral of rR,where r=(x,y,z) and r'=(0,0,0), 
	//ob means observation point, v0 and v1 are the endpoints
	static Point Line_rR1(const Point& ob, const Point& v0, const Point& v1);
	
	//calculate the line integral of r/R,where r=(x,y,z) and r'=(0,0,0).
	//observation site cannot locate on edge v0-v1,though it can locate on the extension of the edge
	static Point Line_rR_1(const Point& ob, const Point& v0, const Point& v1);

	//calculate the line integral of  x^2*R, y^2*R, z^2*R, xz*R,yz*R,xy*R
	static double Line_x2R(const Point& ob,const Point& v0,const Point& v1);
	static double Line_y2R(const Point& ob, const Point& v0, const Point& v1);
	static double Line_z2R(const Point& ob, const Point& v0, const Point& v1);
	static double Line_xzR(const Point& ob, const Point& v0, const Point& v1);
	static double Line_yzR(const Point& ob, const Point& v0, const Point& v1);
	static double Line_xyR(const Point& ob,const Point& v0, const Point& v1);

	//calculate the line integral of  x^2/R, y^2/R, z^2/R, xz/R,yz/R,xy/R
	static double Line_x2R_1(const Point& ob, const Point& v0,const Point& v1);
	static double Line_y2R_1(const Point& ob, const Point& v0,const Point& v1);
	static double Line_z2R_1(const Point& ob, const Point& v0,const Point& v1);
	static double Line_xzR_1(const Point& ob, const Point& v0,const Point& v1);
	static double Line_yzR_1(const Point& ob, const Point& v0,const Point& v1);
	static double Line_xyR_1(const Point& ob, const Point& v0,const Point& v1);

	/*To calculate the line integral of xp*xq*xt/R^3, where xp, xq, xt can be x y or z.
	 *The string s is used to specify the integrand.
	 *Possible values of s are "x3","y3","z3","x2y","x2z","y2x","y2z","z2x","z2y","xyz"
	 *For example, Line_f3R_1(ob, v0, v1,"xyz") mean the integrand is xyz/R
	 **/
	static double Line_f3R_1(const Point& ob, const Point& v0, const Point& v1,string s);





	/***************************************Surface integrals**************************************/
	//calculate the surface integral of R*ds,face_node is the nodes on the face,n is the number of nodes
	static double Surface_R1(const Point& observation, const vector<Node*>& face_node, int n);

	//calculate the surface integral of (1/R)*ds
	static double Surface_R_1(const Point& observation, const vector<Node*>& face_node, int n);

	/*To calculate the surface integral of (rR)*ds.
	 *Return value is a vector of (integral{xR}ds,integral{yR}ds,integral{zR}ds)
	 */
	static Point Surface_rR1(const Point& observation, const vector<Node*>& face_node, int n);

	/*To compute the surface integra of (r/R)*ds
	 *Return value is a vector of (integral{x/R}ds,integral{y/R}ds,integral{z/R}ds)
	 */
	static Point Surface_rR_1(const Point& observation, const vector<Node*>& face_node, int n);

	/*To compute the surface integral of (f^2/R)*ds, where f=x, y or z.
	*Return value is a vector of (integral{x^2/R}ds,integral{y^2/R}ds,integral{z^2/R}ds)
	*/
	static Point Surface_r2R_1(const Point& observation, const vector<Node*>& face_node, int n);

	/*To compute the surface integral of (xz/R)ds */
	static double Surface_xzR_1(const Point& observation, const vector<Node*>& face_node, int n);
	/*To compute the surface integral of (yz/R)ds */
	static double Surface_yzR_1(const Point& observation, const vector<Node*>& face_node, int n);
	/*To compute the surface integral of (xy/R)ds */
	static double Surface_xyR_1(const Point& observation, const vector<Node*>& face_node, int n);
	
	/*To compute the surface integral of (1/R^3)ds */
	static double Surface_R_3(const Point& ob, const vector<Node*>& face_node, int n);

	/*To compute the surface integral of (r/R^3)ds ,r=(x,y,z)
	return a vector of (integral{x/R^3}ds,integral{y/R^3}ds,integral{z/R^3}ds)
	*/
	static Point Surface_rR_3(const Point& ob, const vector<Node*>& face_node, int n);

	//To compute the surface integral of x2/R^3
	static double Surface_x2R_3(const Point& ob, const vector<Node*>& face_node, int n);

	//To compute the surface integral of y2/R^3
	static double Surface_y2R_3(const Point& ob, const vector<Node*>& face_node, int n);

	//To compute the surface integral of z2/R^3
	static double Surface_z2R_3(const Point& ob, const vector<Node*>& face_node, int n);

	//To compute the surface integral of xz/R^3
	static double Surface_xzR_3(const Point& ob, const vector<Node*>& face_node, int n);

	//To compute the surface integral of yz/R^3
	static double Surface_yzR_3(const Point& ob, const vector<Node*>& face_node, int n);

	//To compute the surface integral of xy/R^3
	static double Surface_xyR_3(const Point& ob, const vector<Node*>& face_node, int n);

	//To compute the surface integral of (xp*xq*xt/R)ds, where xp, xq, xt can be x or y or z.
	static void Surface_r3R_1(const Point& ob, const vector<Node*>& face, int n,
		double& x3R_1, double& y3R_1, double &z3R_1,
		double &x2yR_1, double& x2zR_1,
		double &y2xR_1, double& y2zR_1,
		double &z2xR_1, double& z2yR_1,
		double &xyzR_1);

	//To compute the surfae integral of (xp*xq*xt/R^3)ds
	static void Surface_f3R_3(const Point& ob, const vector<Node*>& face, int n,
		double& x3R_3, double& y3R_3, double &z3R_3,
		double &x2yR_3, double& x2zR_3,
		double &y2xR_3, double& y2zR_3,
		double &z2xR_3, double& z2yR_3,
		double &xyzR_3);
};

