#include<fstream>
#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<cstdlib>
using namespace std;
#define EPS 1e-14
int main(int argc, char* argv[]){
	if(argc <5) {
		printf("Usage: %s input_paramters_filename\n", argv[0]);
		return 1;
	}
	//first string: x
	//second string: y
	//third string£º z

	typedef vector<double>::iterator iter;
	int num=0;
	string xstring(argv[1]);
	string ystring(argv[2]);
	string zstring(argv[3]);
	double start_x, xspace, end_x;
	double start_y, yspace, end_y;
	double start_z, zspace, end_z;

	if (xstring.find(':') == string::npos) {
		start_x = atof(argv[1]);
		end_x = atof(argv[1]);
		xspace = 1;
		cout << setw(30) << left << "x coordinate:" << start_x<<endl;
	}
	else {
		if (xstring.find(':') == xstring.rfind(':')) {
			start_x = atof(xstring.substr(0, xstring.find(':')).c_str());
			xspace = 1;
			end_x = atof(xstring.substr(xstring.rfind(':') + 1).c_str());
		}
		else {
			start_x = atof(xstring.substr(0, xstring.find(':')).c_str());
			xspace = atof(xstring.substr(xstring.find(':') + 1, xstring.rfind(':') - xstring.find(':')).c_str());
			end_x = atof(xstring.substr(xstring.rfind(':') + 1).c_str());		
		}
		cout << setw(30) << left << "starting x coordinate:" << start_x << endl;
		cout << setw(30) << left << "x interval:" << xspace << endl;
		cout << setw(30) << left << "endding x coordinate:" << end_x << endl;	
	}


	cout << endl;

	
	if (ystring.find(':') == string::npos) {
		start_y = atof(argv[2]);
		end_y = atof(argv[2]);
		yspace = 1;
		cout << setw(30) << left << "y coordinate:" << start_y<<endl;
	}
	else{
		if (ystring.find(':') == ystring.rfind(':')) {
			start_y = atof(ystring.substr(0, ystring.find(':')).c_str());
			yspace = 1;
			end_y = atof(ystring.substr(ystring.rfind(':') + 1).c_str());
		}
		else {
			start_y = atof(ystring.substr(0, ystring.find(':')).c_str());
			yspace = atof(ystring.substr(ystring.find(':') + 1, ystring.rfind(':') - ystring.find(':')).c_str());
			end_y = atof(ystring.substr(ystring.rfind(':') + 1).c_str());
		}
		cout << setw(30) << left << "starting y coordinate:" << start_y << endl;
		cout << setw(30) << left << "y interval:" << yspace << endl;
		cout << setw(30) << left << "endding y coordinate:" << end_y << endl;
	}

	cout << endl;

	if (zstring.find(':') == string::npos) {
		start_z = atof(argv[3]);
		end_z = atof(argv[3]);
		zspace = 1;
		cout << setw(30) << left << "z coordinate:" << start_y<<endl;
	}
	else{
		if (zstring.find(':') == zstring.rfind(':')) {
			start_z = atof(zstring.substr(0, zstring.find(':')).c_str());
			zspace = 1;
			end_z = atof(zstring.substr(zstring.rfind(':') + 1).c_str());
		}
		else {
			start_z = atof(zstring.substr(0, zstring.find(':')).c_str());
			zspace = atof(zstring.substr(zstring.find(':') + 1, zstring.rfind(':') - zstring.find(':')).c_str());
			end_z = atof(zstring.substr(zstring.rfind(':') + 1).c_str());
		}
		cout << setw(30) << left << "starting z coordinate:" << start_z << endl;
		cout << setw(30) << left << "z interval:" << zspace << endl;
		cout << setw(30) << left << "endding z coordinate:" << end_z << endl;
	}


	//double xspace=0.05,yspace=0.05,start_x=-3,start_y=-3.5,end_x=6,end_y=3.5;
//	double z=atof(argv[3]);
//	cout <<'\n'<< setw(30) << left << "z coordinate:" << z << endl;
	vector<double> x;
	vector<double> y;
	vector<double> z;

	int lab=1;
	ofstream os;
	string outfile(argv[4]);
	os.open(outfile.c_str());


	int cou=0;
	double temp=start_x+cou*xspace;
	do{
		x.push_back(temp);
		cou++;
		temp=start_x+cou*xspace;
	} while ((temp < end_x) || (abs(temp - end_x) < EPS));

	cou=0;
	temp=start_y+cou*yspace;
	do{
		y.push_back(temp);
		cou++;
		temp=start_y+cou*yspace;
	} while ((temp < end_y) || abs(temp - end_y) < EPS);

	cou=0;
	temp=start_z+cou*zspace;
	do{
		z.push_back(temp);
		cou++;
		temp=start_z+cou*zspace;
	} while ((temp < end_z) || abs(temp - end_z) < EPS);


	num=y.size()*x.size()*z.size();
	os<<num<<endl;
	for(iter iz= z.begin();iz!=z.end();++iz){
		for(iter iy=y.begin();iy!=y.end();++iy){
			for(iter ix=x.begin();ix!=x.end();++ix){
				os<<setw(15)<<left<<lab
				  <<setw(15)<<left<<(*ix)
				  <<setw(15)<<left<<(*iy)
				  <<setw(15)<<left<<(*iz)
				  <<'\n';
				  lab++;			  			  
			}
		}
	}
	return 0;
}
