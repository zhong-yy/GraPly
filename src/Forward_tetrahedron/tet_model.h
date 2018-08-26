#pragma once
#include<map>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include"tet_mesh.h"
#include"Observation.h"
#include"timer.h"
class tet_model
{
public:
	tet_model();
	~tet_model();

	//void read_site(const char * name);
	void compute_g();
	void compute_ggt();
	void out_g(const char * name);
	void out_T(const char* name);

	void read_config(string config_file);

	void set_density();
private:
	Observation observation;
	tet_mesh mesh;
	vector<Point> g;
	vector<Dyadic> T;
	map<int, double > region_table_0;
	map<int, std::vector<double> > region_table_1;
	map<int, std::vector<double> > region_table_2;
	map<int, std::vector<double> > region_table_3;
	int density_order;//0 or 1 or 2 or 3

	//delete spaces at the begining and the end, and delete comments
	void line_process(string &line, const string comment_str = "#")
	{
		// if(line.size()==1&&(line[0]=='\n'||line[0]=='\r')){
		// 	cout<<"23333"<<endl;
		// }
		for (char &c : line)
			if (c == '\t' || c == ',' || c == ';'||c=='\r'||c=='\n')
				c = ' ';
		line.erase(0, line.find_first_not_of(" "));
		line.erase(line.find_last_not_of(" ") + 1);
		int n_comment_start = line.find_first_of(comment_str);
		if (n_comment_start != string::npos)
			line.erase(n_comment_start);
	}	
};

