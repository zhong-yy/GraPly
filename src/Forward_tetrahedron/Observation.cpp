#include "Observation.h"

Observation::Observation()
{
}

Observation::~Observation()
{
}

void Observation::read_site(const string &name)
{
	ifstream in_stream;
	in_stream.open(name.c_str());
	assert(in_stream.good());
	string line;
	while (std::getline(in_stream, line))
	{
		line_process(line, "#");
		if (line.empty())
			continue;
		std::istringstream iss(line);
		iss >> n_obs;
		break;
	}
	cout << "\nReading observation sites ...\nThe total number of sites is " << n_obs << ".\n";
	unsigned lab;
	obs.resize(n_obs);
	double x, y, z;
	for (unsigned i = 0; i < n_obs; i++)
	{
		assert(in_stream.good());
		while (std::getline(in_stream, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			iss >> lab>> x >> y >> z;
			break;
		}
		obs[i].setPoint(x, y, z);
	}
	cout << "Reading sites done!\n\n";
}
