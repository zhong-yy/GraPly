#pragma once
/*
*Observation class is to represent a list of observation points and read observation points from ASCII file
*/
#include<vector>
#include<cassert>
#include<fstream>
#include<sstream>
#include"Point.h"
class Observation
{
public:
	Observation();
	~Observation();
	void read_site(const string& name);
	const Point & operator()(unsigned i) {
		assert(i < n_obs);
		return obs[i];
	}
	void set_obs(unsigned i, double x,double y,double z) {
		assert(i<n_obs);
		this->obs[i].setPoint(x, y, z);
	}
	unsigned int get_n_obs() { return n_obs; }
private:
	vector<Point> obs;
	unsigned int n_obs;

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

