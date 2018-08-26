#include<iostream>
#include"Point.h"
#include"Dyadic.h"
#include<vector>
#include <algorithm>
#include"Modeling.h"
#include<ctime>
#include"timer.h"
#include<string>
using namespace std;

int main(int argc, char* argv[]) {
	// clock_t t_start = clock();
	Timer t;
	t.start();
	if (argc <5) {
		printf("Usage: %s input_paramters_filename\n", argv[0]);
		return 1;
	}
	const char *model_file = argv[1];
	const char *coor_file = argv[2];
	const char *out_file = argv[3];
	string specifier = argv[4];
	Modeling m;
	
	m.read_model(model_file);
	m.read_site(coor_file);
	string field;
	if (specifier=="ggt") {
		field="gravity gradient tensor";
		cout << "Computing " <<field<<" ...\n";
		m.analytic_T();
		m.out_T(out_file);
	}
	else if(specifier=="g"){
		field="gravity field";
		cout << "Computing " <<field<<" ...\n";
		m.analytic_g();
		m.out_g(out_file);
	}
	t.stop();

	// clock_t t_end = clock();
	double time=t.getElapsedTimeInMilliSec();
	cout << "Total running time: "<<time<< "ms\n\n";

	return 0;
}


