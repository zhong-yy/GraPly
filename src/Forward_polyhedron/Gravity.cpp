#include "Gravity.h"



Gravity::Gravity()
{
}


Gravity::~Gravity()
{
}

Point Gravity::g_const(const Polyhedral& phl, const Point& ob, const double& rho)
{
	Point g(0, 0, 0);
	//vector<Point> face;
	double iR_1 = 0;
	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++) {
		//	//cout << i<<","<<phl.face_node[i][k];
		//	//(*phl.node[phl.face_node[i][k]]).display();
		//	face[k] = *(phl.node[phl.face_node[i][k]]);
		//}
		iR_1 = Integral::Surface_R_1(ob, phl.faces[i], phl.face_node_number[i]);
		g = g - 1.*G*1e8*rho*(phl.face_normal_vector[i])*iR_1;
	}
	return g;
}

double Gravity::gx_const(const Polyhedral& phl, const Point& ob, const double& rho)
{
	double gx = 0.0;
	//vector<Point> face;
	double iR_1 = 0;
	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++)
		//	face[k] = *(phl.node[phl.face_node[i][k]]);
		iR_1 = Integral::Surface_R_1(ob, phl.faces[i], phl.face_node_number[i]);
		gx = gx - 1.*G*1e8*rho*(phl.face_normal_vector[i].x*1.)*iR_1;
	}
	return gx;
}

double Gravity::gy_const(const Polyhedral& phl, const Point& ob, const double& rho)
{
	double gy = 0.0;
	//vector<Point> face;
	double iR_1 = 0;

	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++)
		//	face[k] = *(phl.node[phl.face_node[i][k]]);
		iR_1 = Integral::Surface_R_1(ob, phl.faces[i], phl.face_node_number[i]);
		gy = gy - 1.*G*1e8*rho*(phl.face_normal_vector[i].y*1.)*iR_1;

	}
	return gy;
}

double Gravity::gz_const(const Polyhedral& phl, const Point& ob, const double& rho)
{
	double gz=0.0;
	//vector<Point> face;
	double iR_1=0;
	//double temp;
	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++)
		//	face[k] = *(phl.node[phl.face_node[i][k]]);
		iR_1 = Integral::Surface_R_1(ob, phl.faces[i], phl.face_node_number[i]);
		gz =gz - 1.*G*1e8*rho*(phl.face_normal_vector[i].z*1.)*iR_1;

	}
	return gz;
}

Point Gravity::g_1st(const Polyhedral& pol, const Point& ob, 
					const double& a, const double& b, const double& c)
{
	Point g(0.0, 0.0, 0.0);
	g = g_1st0(pol, ob, a, b, c);
	g = g + g_const(pol, ob, a*ob.x + b*ob.y + c*ob.z);
	return g;
}

double Gravity::gx_1st(const Polyhedral& pol, const Point& ob,
					   const double& a, const double& b, const double& c)
{
	double gx = 0.0;
	gx = gx_1st0(pol, ob, a, b, c);
	gx = gx + gx_const(pol, ob, a*ob.x + b*ob.y + c*ob.z);
	return gx;
}

double Gravity::gy_1st(const Polyhedral& pol, const Point& ob,
					const double& a,const double& b, const double& c)
{
	double gy = 0.0;
	gy = gy_1st0(pol, ob, a, b, c);
	gy = gy + gy_const(pol, ob, a*ob.x + b*ob.y + c*ob.z);
	return gy;
}

double Gravity::gz_1st(const Polyhedral& pol, const Point& ob,
					   const double& a, const double& b, const double& c)
{
	double g111=0.0;
	g111 = gz_1st0(pol, ob, a, b, c);
	g111 = g111 + gz_const(pol, ob, a*ob.x + b*ob.y + c*ob.z);
	return g111;
}



Point Gravity::g_2nd(const Polyhedral& poh, const Point& ob,
	const double& a, const double& b, const double& c, const double& d, const double& e, const double& f)
{
	Point g(0., 0., 0.);
	double x0 = ob.x, y0 = ob.y, z0 = ob.z;
	g += g_2nd0(poh, ob, a, b, c, d, e, f);
	g += g_1st0(poh, ob, (2.0 * a*x0 + d*z0 + f*y0),
		(2.0 * b*y0 + e*z0 + f*x0),
		(2.0 * c*z0 + d*x0 + e*y0));
	g += g_const(poh, ob, a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*z0 + e*y0*z0 + f*x0*y0);
	return g;
}

double Gravity::gx_2nd(const Polyhedral& poh, const Point& ob,
	const double& a, const double& b, const double& c, const double& d, const double& e, const double& f)
{
	double gx = 0.0;
	double x0 = ob.x, y0 = ob.y, z0 = ob.z;
	gx += gx_2nd0(poh, ob, a, b, c, d, e, f);
	gx += gx_1st0(poh, ob, (2.0 * a*x0 + d*z0 + f*y0),
		(2.0 * b*y0 + e*z0 + f*x0),
		(2.0 * c*z0 + d*x0 + e*y0));
	gx += gx_const(poh, ob, a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*z0 + e*y0*z0 + f*x0*y0);
	return gx;
}

double Gravity::gy_2nd(const Polyhedral& poh, const Point& ob, 
	const double& a, const double& b, const double& c, const double& d, const double& e, const double& f)
{
	double gy = 0.0;
	double x0 = ob.x, y0 = ob.y, z0 = ob.z;
	gy += gy_2nd0(poh, ob, a, b, c, d, e, f);
	gy += gy_1st0(poh, ob, (2.0 * a*x0 + d*z0 + f*y0),
		(2.0 * b*y0 + e*z0 + f*x0),
		(2.0 * c*z0 + d*x0 + e*y0));
	gy += gy_const(poh, ob, a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*z0 + e*y0*z0 + f*x0*y0);
	return gy;
}

double Gravity::gz_2nd(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f)
{
	double gz = 0.0;
	double x0 = ob.x, y0 = ob.y, z0 = ob.z;
	gz += gz_2nd0(poh, ob, a, b, c, d, e, f);
	gz += gz_1st0(poh, ob, (2.0 * a*x0 + d*z0 + f*y0), 
		(2.0 * b*y0 + e*z0 + f*x0), 
		(2.0 * c*z0 + d*x0 + e*y0));
	gz += gz_const(poh, ob, a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*z0 + e*y0*z0 + f*x0*y0);
	return gz;
}

Point Gravity::g_3rd(const Polyhedral& po, const Point& ob, double a[10])
{
	Point g(0., 0., 0.);
	double x0 = ob.x, y0 = ob.y, z0 = ob.z;
	g = g + g_3rd0(po, ob, a);
	g = g + g_2nd0(po, ob, a[7] * z0 + a[8] * y0 + 3 * a[9] * x0,
		a[2] * z0 + 3 * a[3] * y0 + a[6] * x0,
		3 * a[0] * z0 + a[1] * y0 + a[4] * x0,
		2 * a[4] * z0 + a[5] * y0 + 2 * a[7] * x0,
		2 * a[1] * z0 + 2 * a[2] * y0 + a[5] * x0,
		a[5] * z0 + 2 * a[6] * y0 + 2 * a[8] * x0);
	g = g + g_1st0(po, ob, a[4] * z0*z0 + a[5] * y0*z0 + a[6] * y0*y0 + 2 * a[7] * x0*z0 + 2 * a[8] * x0*y0 + 3 * a[9] * x0*x0,
		a[1] * z0*z0 + 2 * a[2] * y0*z0 + 3 * a[3] * y0*y0 + a[5] * x0*z0 + 2 * a[6] * x0*y0 + a[8] * x0*x0,
		3 * a[0] * z0*z0 + a[2] * y0*y0 + 2 * a[1] * y0*z0 + 2 * a[4] * x0*z0 + a[5] * x0*y0 + a[7] * x0*x0);
	g = g + g_const(po, ob, a[0] * z0*z0*z0 + a[1] * y0*z0*z0 + a[2] * y0*y0*z0 + a[3] * y0*y0*y0 
		+ a[4] * x0*z0*z0 + a[5] * x0*y0*z0 + a[6] * x0*y0*y0 
		+ a[7] * x0*x0*z0 + a[8] * x0*x0*y0 + a[9] * x0*x0*x0);
	return g;
}

Point Gravity::g_1st0(const Polyhedral& pol, const Point& ob, double a, double b, double c)
{
	Point A(a, b, c);
	Point g( 0.0,0.0,0.0);
	double h = 0.0;
	//vector<Point> face;
	Point rR_1;

	for (int i = 0; i < pol.face_number; i++) {
		//face.resize(pol.face_node_number[i]);
		//for (int k = 0; k < pol.face_node_number[i]; k++) {
		//	face[k] = *pol.node[pol.face_node[i][k]];
		//}
		const vector<Node*>& face = pol.faces[i];
		h = pol.face_normal_vector[i] * (ob - *face[0]);
		rR_1 = Integral::Surface_rR_1(ob, face, pol.face_node_number[i]);
		g = g + pol.face_normal_vector[i]*(A*rR_1)
			+ 0.5*h*Integral::Surface_R_1(ob, face, pol.face_node_number[i])*A;
	}
	g = -G*1e8*g;//unit:mGal
	return g;
}

double Gravity::gx_1st0(const Polyhedral& pol, const Point& ob, double a, double b, double c)
{
	Point A(a, b, c);
	double gx = 0.0;
	double h = 0.0;
	//vector<Point> face;
	Point rR_1;

	for (int i = 0; i < pol.face_number; i++) {
		//face.resize(pol.face_node_number[i]);
		//for (int k = 0; k < pol.face_node_number[i]; k++) {
		//	face[k] = *pol.node[pol.face_node[i][k]];
		//}
		const vector<Node*>& face = pol.faces[i];
		h = pol.face_normal_vector[i] * (ob - *face[0]);
		rR_1 = Integral::Surface_rR_1(ob, face, pol.face_node_number[i]);
		gx = gx + pol.face_normal_vector[i].x*(A*rR_1)
			+ 0.5*a*h*Integral::Surface_R_1(ob, face, pol.face_node_number[i]);
	}
	gx = -G*1e8*gx;//unit:mGal
	return gx;
}

double Gravity::gy_1st0(const Polyhedral& pol, const Point& ob, double a, double b, double c)
{
	Point A(a, b, c);
	double gy = 0.0;
	double h = 0.0;
	//vector<Point> face;
	Point rR_1;

	for (int i = 0; i < pol.face_number; i++) {
		//face.resize(pol.face_node_number[i]);
		//for (int k = 0; k < pol.face_node_number[i]; k++) {
		//	face[k] = *pol.node[pol.face_node[i][k]];
		//}
		const vector<Node*>& face = pol.faces[i];
		h = pol.face_normal_vector[i] * (ob - *face[0]);
		rR_1 = Integral::Surface_rR_1(ob, face, pol.face_node_number[i]);
		gy = gy + pol.face_normal_vector[i].y*(A*rR_1)
			+ 0.5*b*h*Integral::Surface_R_1(ob, face, pol.face_node_number[i]);
	}
	gy = -G*1e8*gy;//unit:mGal
	return gy;
}

double Gravity::gz_1st0(const Polyhedral& pol, const Point& ob,double a,double b,double c)
{
	Point A(a, b, c);
	double gz=0.0;
	double h=0.0;
	//vector<Point> face;
	Point rR_1;
	for (int i = 0; i < pol.face_number; i++) {
		//face.resize(pol.face_node_number[i]);
		//for (int k = 0; k < pol.face_node_number[i]; k++) {
		//	face[k] = *pol.node[pol.face_node[i][k]];
		//}		

		const vector<Node*>& face = pol.faces[i];
		h = pol.face_normal_vector[i] * (ob - *face[0]);
		rR_1 = Integral::Surface_rR_1(ob, face, pol .face_node_number[i]);
		gz = gz + pol.face_normal_vector[i].z*(A*rR_1) 
			+ 0.5*c*h*Integral::Surface_R_1(ob, face, pol.face_node_number[i]);
	}
	gz = -G*1e8*gz;//unit:mGal
	return gz;
}

Point Gravity::g_2nd0(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f)
{
	Point g(0.0, 0.0, 0.0);
	//vector<Point> face;
	Point r2R_1;
	Point A(a, b, c);
	Point n(0, 0, 0);
	Point x(1., 0., 0.);
	Point y(0., 1., 0.);
	Point z(0., 0., 1.);
	double xzR_1, yzR_1, xyR_1, Rds;
	for (int i = 0; i < poh.face_number; i++) {
		//face.resize(poh.face_node_number[i]);
		//for (int k = 0; k < poh.face_node_number[i]; k++)
		//	face[k] = *poh.node[poh.face_node[i][k]];
		const vector<Node*>& face = poh.faces[i];
		n = poh.face_normal_vector[i];
		r2R_1 = Integral::Surface_r2R_1(ob, face, poh.face_node_number[i]);
		xzR_1 = Integral::Surface_xzR_1(ob, face, poh.face_node_number[i]);
		yzR_1 = Integral::Surface_yzR_1(ob, face, poh.face_node_number[i]);
		xyR_1 = Integral::Surface_xyR_1(ob, face, poh.face_node_number[i]);
		Rds = Integral::Surface_R1(ob, face, poh.face_node_number[i]);
		g = g + n*(A*r2R_1 + d*xzR_1 + e*yzR_1 + f*xyR_1)
			- (x*(2.0*a*n.x + d*n.z + f*n.y)+ y*(2.0*b*n.y + e*n.z + f*n.x)
				+ z*(2.0*c*n.z + d*n.x + e*n.y))*Rds;
	}
	g = -G*1e8*g;//unit:mGal
	return g;
}

double Gravity::gx_2nd0(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f)
{
	double gx = 0.0;
	//vector<Point> face;
	Point r2R_1;
	Point A(a, b, c);
	Point n(0, 0, 0);
	double xzR_1, yzR_1, xyR_1, Rds;
	for (int i = 0; i < poh.face_number; i++) {
		//face.resize(poh.face_node_number[i]);
		//for (int k = 0; k < poh.face_node_number[i]; k++)
		//	face[k] = *poh.node[poh.face_node[i][k]];
		const vector<Node*>& face = poh.faces[i];
		n = poh.face_normal_vector[i];
		r2R_1 = Integral::Surface_r2R_1(ob, face, poh.face_node_number[i]);
		xzR_1 = Integral::Surface_xzR_1(ob, face, poh.face_node_number[i]);
		yzR_1 = Integral::Surface_yzR_1(ob, face, poh.face_node_number[i]);
		xyR_1 = Integral::Surface_xyR_1(ob, face, poh.face_node_number[i]);
		Rds = Integral::Surface_R1(ob, face, poh.face_node_number[i]);
		gx = gx + n.x*(A*r2R_1 + d*xzR_1 + e*yzR_1 + f*xyR_1)
			- (2.0*a*n.x + d*n.z + f*n.y)*Rds;
	}
	gx = -G*1e8*gx;//unit:mGal
	return gx;
}

double Gravity::gy_2nd0(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f)
{
	double gy = 0.0;
	//vector<Point> face;
	Point r2R_1;
	Point A(a, b, c);
	Point n(0, 0, 0);
	double xzR_1, yzR_1, xyR_1, Rds;
	for (int i = 0; i < poh.face_number; i++) {
		//face.resize(poh.face_node_number[i]);
		//for (int k = 0; k < poh.face_node_number[i]; k++)
		//	face[k] = *poh.node[poh.face_node[i][k]];
		const vector<Node*>& face = poh.faces[i];
		n = poh.face_normal_vector[i];
		r2R_1 = Integral::Surface_r2R_1(ob, face, poh.face_node_number[i]);
		xzR_1 = Integral::Surface_xzR_1(ob, face, poh.face_node_number[i]);
		yzR_1 = Integral::Surface_yzR_1(ob, face, poh.face_node_number[i]);
		xyR_1 = Integral::Surface_xyR_1(ob, face, poh.face_node_number[i]);
		Rds = Integral::Surface_R1(ob, face, poh.face_node_number[i]);
		gy = gy + n.y*(A*r2R_1 + d*xzR_1 + e*yzR_1 + f*xyR_1)
			- (2.0*b*n.y + e*n.z + f*n.x)*Rds;
	}
	gy = -G*1e8*gy;//unit:mGal
	return gy;
}

double Gravity::gz_2nd0(const Polyhedral& poh, const Point& ob, double a, double b, double c, double d, double e, double f)
{
	double gz = 0.0;
	//vector<Point> face;
	Point r2R_1;
	Point A(a, b, c);
	Point n(0, 0, 0);
	double xzR_1, yzR_1, xyR_1,Rds;
	for (int i = 0; i < poh.face_number; i++) {
		//face.resize(poh.face_node_number[i]);
		//for (int k = 0; k < poh.face_node_number[i]; k++)
		//	face[k] = *poh.node[poh.face_node[i][k]];
		const vector<Node*>& face = poh.faces[i];
		n = poh.face_normal_vector[i];
		r2R_1 = Integral::Surface_r2R_1(ob, face, poh.face_node_number[i]);
		xzR_1=Integral::Surface_xzR_1(ob, face, poh.face_node_number[i]);
		yzR_1=Integral::Surface_yzR_1(ob, face, poh.face_node_number[i]);
		xyR_1=Integral::Surface_xyR_1(ob, face, poh.face_node_number[i]);
		Rds=Integral::Surface_R1(ob, face, poh.face_node_number[i]);
		gz = gz + n.z*(A*r2R_1+d*xzR_1+e*yzR_1+f*xyR_1) 
			- (2.0*c*n.z+d*n.x+e*n.y)*Rds;
	}
	gz = -G*1e8*gz;//unit:mGal
	return gz;
}

Point Gravity::g_3rd0(const Polyhedral& po, const Point& ob, double a[10])
{
	Point g(0.0, 0.0, 0.0);
	Point n(0, 0, 0);
	//vector<Point> face;
	double Rds = 0.0;
	Point rRds(0.0,0.0,0.0);
	double r3ds[10];
	double lamd_R = 0.0,h=0.0;
	for (int i = 0; i < po.face_number; i++) {
		//face.resize(po.face_node_number[i]);
		//for (int k = 0; k < po.face_node_number[i]; k++)
		//	face[k] = *po.node[po.face_node[i][k]];
		const vector<Node*>& face = po.faces[i];
		n = po.face_normal_vector[i];
		Integral::Surface_r3R_1(ob, face, po.face_node_number[i],
			r3ds[9],r3ds[3],r3ds[0],
			r3ds[8],r3ds[7],
			r3ds[6],r3ds[2],
			r3ds[4],r3ds[1],
			r3ds[5]);
		lamd_R = 0.0;
		for (int j = 0; j < 10; j++) {
			lamd_R = lamd_R + a[j] * r3ds[j];
		}
		//cout << "AA" << endl;
	
		g = g + n*lamd_R;
		h = (ob - *face[0])*n;
		rRds = Integral::Surface_rR1(ob, face, po.face_node_number[i]);
		Rds = Integral::Surface_R1(ob, face, po.face_node_number[i]);
		g.x = g.x -( a[4]*(n.z*rRds.z+0.25*h*Rds)+a[5]*(n.z*rRds.y)
			+a[6]*(n.y*rRds.y+0.25*h*Rds) + 2 * a[7] * (n.z*rRds.x)
			+2*a[8]*(n.y*rRds.x)+3*a[9]*(n.x*rRds.x+0.25*h*Rds) );
		g.y=g.y -( a[1] * (n.z*rRds.z + 0.25*h*Rds) + 2*a[2] * (n.z*rRds.y)
			+ 3*a[3] * (n.y*rRds.y + 0.25*h*Rds) + a[5] * (n.z*rRds.x)
			+ 2 * a[6] * (n.y*rRds.x) +  a[8] * (n.x*rRds.x + 0.25*h*Rds));
		g.z=g.z -( 3*a[0] * (n.z*rRds.z + 0.25*h*Rds) + 2*a[1] * (n.z*rRds.y)
			+ a[2] * (n.y*rRds.y + 0.25*h*Rds) + 2 * a[4] * (n.z*rRds.x)
			+ a[5] * (n.y*rRds.x) + a[7] * (n.x*rRds.x + 0.25*h*Rds));
	}
	g = -G*1e8*g;//unit:mGal
	return g;
}

Dyadic Gravity::tensor_const(const Polyhedral& phl, const Point& ob, const double& rho)
{
	Dyadic T;
	//double S_R_1;

	//vector<Point> face;
	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++)
		//	face[k] = *phl.node[phl.face_node[i][k]];
		T = T + Dyadic(Integral::Surface_rR_3(ob, phl.faces[i], phl.face_node_number[i]),phl.face_normal_vector[i]);
	}
	T = T*(-G*rho);//µ¥Î»s^{-2}
	return T;
}

Dyadic Gravity::tensor_1st(const Polyhedral& phl, const Point& ob, 
	                       const double& a, const double& b, const double& c)
{
	Dyadic T;
	T = tensor_1st0(phl, ob, a, b, c);
	T += tensor_const(phl, ob, a*ob.x + b*ob.y + c*ob.z);
	return T;
}

Dyadic Gravity::tensor_2nd(const Polyhedral& phl, const Point& ob, double a, double b, double c, double d, double e, double f)
{
	Dyadic T;
	double x0 = ob.x, y0 = ob.y, z0 = ob.z;
	T += tensor_2nd0(phl, ob, a, b, c, d, e, f);
	T += tensor_1st0(phl, ob, (2.0 * a*x0 + d*z0 + f*y0),
		(2.0 * b*y0 + e*z0 + f*x0),
		(2.0 * c*z0 + d*x0 + e*y0));
	T += tensor_const(phl, ob, a*x0*x0 + b*y0*y0 + c*z0*z0 + d*x0*z0 + e*y0*z0 + f*x0*y0);
	return T;
}

Dyadic Gravity::tensor_3rd(const Polyhedral& po, const Point& ob, double a[10])
{
	Dyadic T;
	double x0 = ob.x, y0 = ob.y, z0 = ob.z;
	//cout << "here";
	T = T + tensor_3rd0(po, ob, a);
	//cout << "here";
	T = T + tensor_2nd0(po, ob, a[7] * z0 + a[8] * y0 + 3 * a[9] * x0,
		a[2] * z0 + 3 * a[3] * y0 + a[6] * x0,
		3 * a[0] * z0 + a[1] * y0 + a[4] * x0,
		2 * a[4] * z0 + a[5] * y0 + 2 * a[7] * x0,
		2 * a[1] * z0 + 2 * a[2] * y0 + a[5] * x0,
		a[5] * z0 + 2 * a[6] * y0 + 2 * a[8] * x0);
	T = T + tensor_1st0(po, ob, a[4] * z0*z0 + a[5] * y0*z0 + a[6] * y0*y0 + 2 * a[7] * x0*z0 + 2 * a[8] * x0*y0 + 3 * a[9] * x0*x0,
		a[1] * z0*z0 + 2 * a[2] * y0*z0 + 3 * a[3] * y0*y0 + a[5] * x0*z0 + 2 * a[6] * x0*y0 + a[8] * x0*x0,
		3 * a[0] * z0*z0 + a[2] * y0*y0 + 2 * a[1] * y0*z0 + 2 * a[4] * x0*z0 + a[5] * x0*y0 + a[7] * x0*x0);
	T = T + tensor_const(po, ob, a[0] * z0*z0*z0 + a[1] * y0*z0*z0 + a[2] * y0*y0*z0 + a[3] * y0*y0*y0
		+ a[4] * x0*z0*z0 + a[5] * x0*y0*z0 + a[6] * x0*y0*y0
		+ a[7] * x0*x0*z0 + a[8] * x0*x0*y0 + a[9] * x0*x0*x0);
	return T;
}

Dyadic Gravity::tensor_1st0(const Polyhedral& phl, const Point& ob, double a, double b, double c)
{
	Dyadic T,T1,T2;
	double S_R_1 = 0.0,hi=0;
	Point ni(0., 0., 0.),A(a,b,c),mj(0.,0.,0.);
	Point S_lamda(0, 0, 0);
	//vector<Point> face;
	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++)
		//	face[k] = *phl.node[phl.face_node[i][k]];
		const vector<Node*>& face = phl.faces[i];
		S_R_1 = Integral::Surface_R_1(ob, face, phl.face_node_number[i]);
		ni = phl.face_normal_vector[i];
		
		T1 += Dyadic(ni, A*S_R_1);
		/*******************************************/
		
		hi = ni* (ob - *face[0]);
		S_lamda.zero();
		S_lamda += ni*(hi*(A*Integral::Surface_rR_3(ob, face, phl.face_node_number[i])));
		S_lamda += ((ni*A)*ni-A)*S_R_1;

		//face.push_back(face[0]);
		unsigned int next = 0;
		for (int k = 0; k < phl.face_node_number[i]; k++) {
			next = (k + 1) % (phl.face_node_number[i]);
			mj = unitCross((*face[next] - *face[k]), ni);
			S_lamda =S_lamda+mj*(A*(Integral::Line_rR_1(ob, *face[k], *face[next])));
		}
		T2 += Dyadic(S_lamda, ni);
	}
	T = T1 - T2;
	T = T*(-G);
	return T;
}

Dyadic Gravity::tensor_2nd0(const Polyhedral& phl, const Point& ob, double a, double b, double c, double d, double e, double f)
{
	Dyadic T, T1, T2;
	double S_R_1 = 0.0, hi = 0;
	Point ni(0., 0., 0.), mj(0., 0., 0.);
	Point S_lamda(0, 0, 0);
	Point S_rR_1(0, 0, 0);
	Dyadic A(2 * a, f, d, f, 2 * b, e, d, e, 2 * c);
	Point A1 = A[0];
	Point A2 = A[1];
	Point A3 = A[2];
	//Point A1(2 * a, f, d);
	//Point A2(f, 2 * b, e);
	//Point A3(d, e, 2 * c);
	double S_x2R_3 = 0.0, S_y2R_3 = 0.0, S_z2R_3 = 0.0, S_xzR_3 = 0.0, S_yzR_3 = 0.0, S_xyR_3 = 0.0;
	//vector<Point> face;
	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++)
			//face[k] = *phl.node[phl.face_node[i][k]];
		const vector<Node*>& face = phl.faces[i];
		S_R_1 = Integral::Surface_R_1(ob, face, phl.face_node_number[i]);

		ni = phl.face_normal_vector[i];
		hi = ni* (ob - *face[0]);


		S_rR_1 = Integral::Surface_rR_1(ob, face, phl.face_node_number[i]);
		T1 = T1 + Dyadic(ni, Point(A1*S_rR_1, A2*S_rR_1, A3*S_rR_1));
		T1 = T1 + A*(0.5*hi*S_R_1);
		
		S_x2R_3 = Integral::Surface_x2R_3(ob, face, phl.face_node_number[i]);
		S_y2R_3 = Integral::Surface_y2R_3(ob, face, phl.face_node_number[i]);
		S_z2R_3 = Integral::Surface_z2R_3(ob, face, phl.face_node_number[i]);
		S_xzR_3 = Integral::Surface_xzR_3(ob, face, phl.face_node_number[i]);
		S_yzR_3 = Integral::Surface_yzR_3(ob, face, phl.face_node_number[i]);
		S_xyR_3 = Integral::Surface_xyR_3(ob, face, phl.face_node_number[i]);
		S_lamda.zero();
		S_lamda = S_lamda + ni*hi*(a*S_x2R_3 + b*S_y2R_3 + c*S_z2R_3 + d*S_xzR_3 + e*S_yzR_3 + f*S_xyR_3);
		S_lamda += (Point(A1*S_rR_1, A2*S_rR_1, A3*S_rR_1)*ni)*ni - Point(A1*S_rR_1, A2*S_rR_1, A3*S_rR_1);
		//face.push_back(face[0]);

		unsigned int next = 0;
		for (int k = 0; k < phl.face_node_number[i]; k++) {
			next = (k + 1) % (phl.face_node_number[i]);
			mj = unitCross((*face[next] - *face[k]), ni);
			S_lamda = S_lamda + mj*(a*Integral::Line_x2R_1(ob, *face[k], *face[next])
				+b*Integral::Line_y2R_1(ob, *face[k], *face[next])
				+c*Integral::Line_z2R_1(ob, *face[k], *face[next])
				+d*Integral::Line_xzR_1(ob, *face[k], *face[next])
				+e*Integral::Line_yzR_1(ob, *face[k], *face[next])
				+f*Integral::Line_xyR_1(ob, *face[k], *face[next]));
		}
		T2 = T2 + Dyadic(S_lamda, ni);
	}
	T = T1 - T2;
	T = T*(-G);
	return T;
}

Dyadic Gravity::tensor_3rd0(const Polyhedral& phl, const Point& ob, double a[10])
{
	Dyadic T, T1, T2;
	Point ni(0, 0, 0), mj(0, 0, 0);
	Point rR_1dv(0, 0, 0);
	Point gragra_lam_xx(6 * a[9], 2 * a[8], 2 * a[7]);
	Point gragra_lam_yx(2 * a[8], 2 * a[6], a[5]);
	Point gragra_lam_zx(2 * a[7], a[5], 2 * a[4]);
	Point gragra_lam_yy(2 * a[6], 6 * a[3], 2 * a[2]);
	Point gragra_lam_zy(a[5], 2 * a[2], 2 * a[1]);
	Point gragra_lam_zz(2*a[4], 2 * a[1], 6 * a[0]);

	Point S_gra_lam(0, 0, 0);

	Point S_lamda(0, 0, 0);

	double z2R_1 = 0.0, yzR_1 = 0.0, y2R_1 = 0.0, xzR_1 = 0.0, xyR_1 = 0.0, x2R_1 = 0.0;
	Point r2R_1;

	//vector<Point> face;
	int num=0;
	Point graS;
	double hi = 0.0;
	for (int i = 0; i < phl.face_number; i++) {
		//face.resize(phl.face_node_number[i]);
		//for (int k = 0; k < phl.face_node_number[i]; k++)
		//	face[k] = *phl.node[phl.face_node[i][k]];
		const vector<Node*>& face = phl.faces[i];
		ni = phl.face_normal_vector[i];
		hi = ni* (ob - *face[0]);

		num = phl.face_node_number[i];
		r2R_1 = Integral::Surface_r2R_1(ob, face, num);
		x2R_1 = r2R_1.x;
		y2R_1 = r2R_1.y;
		z2R_1 = r2R_1.z;
		xzR_1 = Integral::Surface_xzR_1(ob, face, num);
		yzR_1 = Integral::Surface_yzR_1(ob, face, num);
		xyR_1 = Integral::Surface_xyR_1(ob, face, num);
		
		graS.setPoint(a[4] * z2R_1 + a[5] * yzR_1 + a[6] * y2R_1 + 2 * a[7] * xzR_1 + 2*a[8] * xyR_1 + 3 * a[9] * x2R_1,
				a[1] * z2R_1 + 2 * a[2] * yzR_1 + 3 * a[3] * y2R_1 + a[5] * xzR_1 + 2 * a[6] * xyR_1 + a[8] * x2R_1,
				3 * a[0] * z2R_1 + 2 * a[1] * yzR_1 + a[2] * y2R_1 + 2 * a[4] * xzR_1 + a[5] * xyR_1 + a[7] * x2R_1);
		T1 = T1 + Dyadic(ni, graS);
		rR_1dv += ni*Integral::Surface_R1(ob, face, num);

		S_lamda.zero();
		double x3R_3, y3R_3, z3R_3, x2yR_3, x2zR_3, y2xR_3, y2zR_3, z2xR_3, z2yR_3,xyzR_3;
		Integral::Surface_f3R_3(ob,face, num,
			x3R_3, y3R_3, z3R_3,
			x2yR_3,x2zR_3,
			y2xR_3,y2zR_3,
			z2xR_3,z2yR_3,
			xyzR_3);
		S_lamda = S_lamda + ni*hi*(a[0] * z3R_3 + a[1] * z2yR_3 + a[2] * y2zR_3 + a[3] * y3R_3 + a[4] * z2xR_3 + a[5] * xyzR_3 +
			a[6] * y2xR_3 + a[7] * x2zR_3 + a[8] * x2yR_3 + a[9] * x3R_3);
		S_lamda = S_lamda + (graS*ni)*ni - graS;
		//face.push_back(face[0]);
		unsigned int next = 0;
		for (int k = 0; k < phl.face_node_number[i]; k++) {
			next = (k + 1) % phl.face_node_number[i];
			mj = unitCross((*face[next] - *face[k]), ni);
			S_lamda = S_lamda + mj*(a[0] * Integral::Line_f3R_1(ob, *face[k], *face[next], "z3")
				+ a[1] * Integral::Line_f3R_1(ob, *face[k], *face[next], "z2y")
				+ a[2] * Integral::Line_f3R_1(ob, *face[k], *face[next], "y2z")
				+ a[3] * Integral::Line_f3R_1(ob, *face[k], *face[next], "y3")
				+ a[4] * Integral::Line_f3R_1(ob, *face[k], *face[next], "z2x")
				+ a[5] * Integral::Line_f3R_1(ob, *face[k], *face[next], "xyz")
				+ a[6] * Integral::Line_f3R_1(ob, *face[k], *face[next], "y2x")
				+ a[7] * Integral::Line_f3R_1(ob, *face[k], *face[next], "x2z")
				+ a[8] * Integral::Line_f3R_1(ob, *face[k], *face[next], "x2y")
				+ a[9] * Integral::Line_f3R_1(ob, *face[k], *face[next], "x3"));
		}
		T2 = T2 + Dyadic(S_lamda, ni);
	}
	T1 = T1 - Dyadic(gragra_lam_xx*rR_1dv, gragra_lam_yx*rR_1dv, gragra_lam_zx*rR_1dv,
		gragra_lam_yx*rR_1dv, gragra_lam_yy*rR_1dv, gragra_lam_zy*rR_1dv,
		gragra_lam_zx*rR_1dv, gragra_lam_zy*rR_1dv, gragra_lam_zz*rR_1dv);
	T = T1 - T2;
	T = T*(-G);
	return T;
}
