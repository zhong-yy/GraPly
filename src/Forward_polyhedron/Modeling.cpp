#include "Modeling.h"
#include <ctime>

Modeling::Modeling()
{
}

Modeling::~Modeling()
{
	assert(num_po == po.size());
	for (int i = 0; i < num_po; i++)
	{
		assert(po[i].node_number == po[i].node.size());
		for (int k = 0; k < po[i].node_number; k++)
		{
			delete po[i].node[k];
		}
	}
}

void Modeling::read_model(const char *name)
{
	ifstream is;
	is.open(name);
	unsigned lab;
	assert(is.good());
	std::string line;
	cout << "Reading model file: \"" << string(name) <<"\""<< endl;
	while (std::getline(is, line))
	{
		line_process(line, "#");
		if (line.empty())
			continue;
		std::istringstream iss(line);
		iss >> num_po;
		break;
	}
	po.resize(num_po);
	cout << "\nNumber of polyhedrons: " << num_po << "\n\n";

	for (unsigned i = 0; i < num_po; i++)
	{
		cout << "Polyhedral " << i << endl;
		// assert(is.good());
		while (std::getline(is, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			iss >> lab;
			// cout<<"lab="<<lab<<",i="<<i<<endl;
			assert(lab == i);
			iss >> po[i].node_number; //total number of vertices in the polyhedral
			po[i].node.resize(po[i].node_number);
			// cout << po[i].node.size() << endl;
			for (int kk = 0; kk < po[i].node_number; kk++)
			{
				po[i].node[kk] = NULL;
			}

			iss >> po[i].face_number;
			cout << "     Number of nodes : " << po[i].node_number << "\n"
				 << "     Number of facets: " << po[i].face_number << "\n\n";
			po[i].face_node_number.resize(po[i].face_number);
			po[i].face_node.resize(po[i].face_number);
			po[i].face_normal_vector.resize(po[i].face_number);
			break;
		}

		unsigned temp;

		for (int j = 0; j < po[i].node_number; j++)
		{
			assert(is.good());
			while (std::getline(is, line))
			{
				line_process(line, "#");
				if (line.empty())
					continue;

				std::istringstream iss(line);
				iss >> temp;
				//coordinate of each node
				double x, y, z;
				iss >> x >> y >> z;
				po[i].node[j] = new Node(x, y, z);
				break;
			}
		}
		// for (int j = 0; j < po[i].node_number; j++) {
		// 	cout << "node " << j << ": (x, y, z)=" << (*po[i].node[j]).x <<
		// 		", " << (*po[i].node[j]).y << ", " << (*po[i].node[j]).z << endl;
		// 	// (*po[i].node[j]).display();
		// }

		for (int j = 0; j < po[i].face_number; j++)
		{
			assert(is.good());
			while (std::getline(is, line))
			{
				line_process(line, "#");
				if (line.empty())
					continue;
				std::istringstream iss(line);
				iss >> temp;
				iss >> po[i].face_node_number[j];
				po[i].face_node[j].resize(po[i].face_node_number[j]);
				// cout << "Face " << j << ": " << "(";
				for (int k = 0; k < po[i].face_node_number[j]; k++)
				{
					iss >> po[i].face_node[j][k];
					// cout << po[i].face_node[j][k] << ", ";
				}

				// cout << ")\n";
				// for (int k = 0; k < po[i].face_node_number[j]; k++) {
				// 	cout << po[i].face_node[j][k] << ",";
				// 	(*po[i].node[po[i].face_node[j][k]]).display();
				// }
				break;
			}
		}

		assert(is.good());
		while (std::getline(is, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			iss >> po[i].const_density;
			break;
		}

		//cout << po[i].const_density<<endl;

		assert(is.good());
		while (std::getline(is, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			for (int j = 0; j < 3; j++)
			{
				iss >> po[i].lin[j];
				//cout << po[i].lin[j]<<"\t";
			}
			break;
		}

		//cout << endl;
		assert(is.good());
		while (std::getline(is, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			for (int j = 0; j < 6; j++)
			{
				iss >> po[i].qua[j];
				//cout << po[i].qua[j] << "\t";
			}
			break;
		}
		//cout << endl;

		while (std::getline(is, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			for (int j = 0; j < 10; j++)
			{
				iss >> po[i].cub[j];
				//cout << po[i].cub[j]<<"\t";
			}
			break;
		}

		//cout << endl;

		po[i].init();
		//po[i].sort_and_compute_normal_vector();//sort the order of face nodes

		// to check whether all vertice on a face are coplanar, if not, exit
		for (int j = 0; j < po[i].face_number; j++)
		{
			if (po[i].face_node_number[j] > 3)
			{
				for (int k = 3; k < po[i].face_node_number[j]; k++)
				{
					double flag = po[i].face_normal_vector[j] * (*po[i].node[po[i].face_node[j][k]] - *po[i].node[po[i].face_node[j][0]]);
					if (abs(flag) > TOLERANCE)
					{
						cerr << "The " << k << "-th vertex on the " << j << "-th facet of the " << i << "-th polyhedron is not coplanar with the first three point on face" << j << endl;
						exit(1);
					}
				}
			}
		}
	}
}

void Modeling::read_site(const char *name)
{
	ifstream is;
	is.open(name);
	string line;
	assert(is.good());
	cout << "Reading computation points in file: \"" << string(name) << "\"" << endl;
	while (std::getline(is, line))
	{
		line_process(line, "#");
		if (line.empty())
			continue;
		std::istringstream iss(line);
		iss >> num_site;
		break;
	}

	cout << "Number of computation points: " << num_site << endl
		 << endl;
	unsigned lab;
	coor.resize(num_site);
	g.resize(num_site);
	T.resize(num_site);
	//gx.resize(num_site);
	//gy.resize(num_site);
	//gz.resize(num_site);
	double x, y, z;
	for (unsigned i = 0; i < num_site; i++)
	{
		assert(is.good());
		while (std::getline(is, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			iss >> lab;
			iss >> x >> y >> z;
			break;
		}
		coor[i].setPoint(x, y, z);
	}
}

void Modeling::out_g(const char *name)
{
	ofstream os;
	os.open(name);
	assert(os.good());

	os << setw(25) << left << "x(km)"
	   << setw(25) << left << "y(km)"
	   << setw(25) << left << "z(km)"
	   << setw(25) << left << "gravity_x(mGal)"
	   << setw(25) << left << "gravity_y(mGal)"
	   << setw(25) << left << "gravity_z(mGal)";
	os << '\n';
	for (unsigned i = 0; i < num_site; i++)
	{
		assert(os.good());
		os << setprecision(7);
		os << fixed;
		os << setw(25) << left << coor[i].x
		   << setw(25) << left << coor[i].y
		   << setw(25) << left << coor[i].z;
		os << setprecision(14);
		os << scientific;
		os << setw(25) << left << g[i].x;
		os << setw(25) << left << g[i].y;
		os << setw(25) << left << g[i].z << '\n';
	}
}

void Modeling::out_T(const char *name)
{
	ofstream os;
	os.open(name);
	assert(os.good());

	os << setw(25) << left << "x(km)"
	   << setw(25) << left << "y(km)"
	   << setw(25) << left << "z(km)"
	   << setw(25) << left << "Txx(1/s^2)"
	   << setw(25) << left << "Tyx(1/s^2)"
	   << setw(25) << left << "Tzx(1/s^2)"
	   << setw(25) << left << "Txy(1/s^2)"
	   << setw(25) << left << "Tyy(1/s^2)"
	   << setw(25) << left << "Tzy(1/s^2)"
	   << setw(25) << left << "Txz(1/s^2)"
	   << setw(25) << left << "Tyz(1/s^2)"
	   << setw(25) << left << "Tzz(1/s^2)";
	os << '\n';
	for (unsigned i = 0; i < num_site; i++)
	{
		assert(os.good());
		os << setprecision(7);
		os << fixed;
		os << setw(25) << left << coor[i].x
		   << setw(25) << left << coor[i].y
		   << setw(25) << left << coor[i].z;
		os << setprecision(14);
		os << scientific;
		os << setw(25) << left << T[i].xx;
		os << setw(25) << left << T[i].yx;
		os << setw(25) << left << T[i].zx;
		os << setw(25) << left << T[i].xy;
		os << setw(25) << left << T[i].yy;
		os << setw(25) << left << T[i].zy;
		os << setw(25) << left << T[i].xz;
		os << setw(25) << left << T[i].yz;
		os << setw(25) << left << T[i].zz << '\n';
	}
}

void Modeling::analytic_g()
{
	vector<Gravity> gra;
	gra.resize(num_site);
#pragma omp parallel for
	for (int i = 0; i < num_site; i++)
	{
		g[i].setPoint(0., 0., 0.);
		for (int j = 0; j < num_po; j++)
		{
			//clock_t t_start = clock();
			g[i] += gra[i].g_const(po[j], coor[i], po[j].const_density);
			g[i] += gra[i].g_1st(po[j], coor[i], po[j].lin[0], po[j].lin[1], po[j].lin[2]);
			g[i] += gra[i].g_2nd(po[j], coor[i], po[j].qua[0], po[j].qua[1], po[j].qua[2], po[j].qua[3], po[j].qua[4], po[j].qua[5]);
			g[i] += gra[i].g_3rd(po[j], coor[i], po[j].cub);
		}
	}
}
void Modeling::analytic_T()
{
	vector<Gravity> gra;
	gra.resize(num_site);
#pragma omp parallel for
	for (int i = 0; i < num_site; i++)
	{
		T[i].set();
		for (int j = 0; j < num_po; j++)
		{
			//clock_t t_start = clock();
			T[i] = T[i] + gra[i].tensor_const(po[j], coor[i], po[j].const_density);
			T[i] = T[i] + gra[i].tensor_1st(po[j], coor[i], po[j].lin[0], po[j].lin[1], po[j].lin[2]);
			T[i] += gra[i].tensor_2nd(po[j], coor[i], po[j].qua[0], po[j].qua[1], po[j].qua[2], po[j].qua[3], po[j].qua[4], po[j].qua[5]);
			T[i] += gra[i].tensor_3rd(po[j], coor[i], po[j].cub);
		}
	}
}
