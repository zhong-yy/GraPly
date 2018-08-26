#include<iostream>

#include<vector>
#include <algorithm>

#include<ctime>
#include<string>
#include"Point.h"
#include"Dyadic.h"
#include"Tetrahedron.h"
#include"tet_model.h"

using namespace std;

int main(int argc, char* argv[]) {
	clock_t t_start = clock();
	if (argc <2) {
		printf("Usage: %s input_paramters_filename\n", argv[0]);
		return 1;
	}
	/**************************command parameters*******************************************/
	//the first input parameter specify the model file
	const char *config = argv[1];
	tet_model model;
	const string out_g_file("g.dat");
	const string out_ggt_file("ggt.dat");

	model.read_config(config);
	model.compute_g();
	model.out_g(out_g_file.c_str());
	cout << "Calculated gravity data has been written into \"" << out_g_file <<"\""<< endl;

	

	model.compute_ggt();
	model.out_T(out_ggt_file.c_str());
	cout << "Calculated GGT data has been written into \"" << out_ggt_file << "\"" << endl;
	return 0;
}


