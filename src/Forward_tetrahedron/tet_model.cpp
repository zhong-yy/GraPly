#include "tet_model.h"

tet_model::tet_model()
{
}

tet_model::~tet_model()
{
}

void tet_model::compute_g()
{
	vector<Gravity> gra;
	unsigned int num_site = observation.get_n_obs();
	unsigned int num_tet = mesh.get_n_tets();
	const vector<Tetrahedron *> &tetrahedrons = mesh.get_tets();
	gra.resize(num_site);

	Timer ti;
	ti.start();
	//cout << (*tetrahedrons[1]).lin[0] << (*tetrahedrons[1]).lin[1] << (*tetrahedrons[1]).lin[2] << endl;
	if (density_order == 0)
	{
		cout << "\nThe density contrast is constant.\n";
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			g[i].setPoint(0., 0., 0.);
			for (unsigned j = 0; j < num_tet; j++)
			{
				//clock_t t_start = clock();
				g[i] += gra[i].g_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);
			}
		}
	}
	else if (density_order == 1)
	{
		cout << "\nThe density contrast is linear.\n";
		cout << "Computing gravity ..." << endl;
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			g[i].setPoint(0., 0., 0.);
			for (unsigned j = 0; j < num_tet; j++)
			{
				//clock_t t_start = clock();
				g[i] += gra[i].g_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);
				g[i] += gra[i].g_1st(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).lin[0],
									 (*tetrahedrons[j]).lin[1], (*tetrahedrons[j]).lin[2]);
			}
		}
	}
	else if (density_order == 2)
	{
		cout << "\nThe density contrast is quadratic.\n";
		cout << "Computing gravity ..." << endl;
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			g[i].setPoint(0., 0., 0.);
			for (unsigned j = 0; j < num_tet; j++)
			{
				//clock_t t_start = clock();
				g[i] += gra[i].g_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);

				g[i] += gra[i].g_1st(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).lin[0],
									 (*tetrahedrons[j]).lin[1], (*tetrahedrons[j]).lin[2]);

				g[i] += gra[i].g_2nd(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).qua[0],
									 (*tetrahedrons[j]).qua[1], (*tetrahedrons[j]).qua[2], (*tetrahedrons[j]).qua[3],
									 (*tetrahedrons[j]).qua[4], (*tetrahedrons[j]).qua[5]);
			}
		}
	}
	else if (density_order == 3)
	{
		cout << "\nThe density contrast is cubic.\n";
		cout << "Computing gravity ..." << endl;
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			g[i].setPoint(0., 0., 0.);
			for (unsigned j = 0; j < num_tet; j++)
			{
				//clock_t t_start = clock();
				g[i] += gra[i].g_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);

				g[i] += gra[i].g_1st(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).lin[0],
									 (*tetrahedrons[j]).lin[1], (*tetrahedrons[j]).lin[2]);

				g[i] += gra[i].g_2nd(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).qua[0],
									 (*tetrahedrons[j]).qua[1], (*tetrahedrons[j]).qua[2], (*tetrahedrons[j]).qua[3],
									 (*tetrahedrons[j]).qua[4], (*tetrahedrons[j]).qua[5]);

				g[i] += gra[i].g_3rd(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).cub);
			}
		}
	}
	ti.stop();
	cout << "Computing gravity vectors done! ";
	cout << "Computation time for gravity vectors is " << ti.getElapsedTimeInSec() << " second(s).\n\n";
}

void tet_model::compute_ggt()
{
	vector<Gravity> gra;
	unsigned int num_site = observation.get_n_obs();
	unsigned int num_tet = mesh.get_n_tets();
	const vector<Tetrahedron *> &tetrahedrons = mesh.get_tets();
	gra.resize(num_site);

	Timer ti;
	ti.start();
	if (density_order == 0)
	{
		cout << "\nThe density contrast is constant.\n";
		cout << "Computing gravity gradien tensor ..." << endl;
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			T[i].set();
			for (unsigned j = 0; j < num_tet; j++)
			{
				T[i] = T[i] + gra[i].tensor_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);
			}
		}
	}
	else if (density_order == 1)
	{
		cout << "\nThe density contrast is linear.\n";
		cout << "Computing gravity gradien tensor ..." << endl;
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			T[i].set();
			for (unsigned j = 0; j < num_tet; j++)
			{
				//clock_t t_start = clock();
				T[i] = T[i] + gra[i].tensor_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);
				T[i] = T[i] + gra[i].tensor_1st(*tetrahedrons[j], observation(i),
												(*tetrahedrons[j]).lin[0], (*tetrahedrons[j]).lin[1], (*tetrahedrons[j]).lin[2]);
			}
		}
	}
	else if (density_order == 2)
	{
		cout << "\nThe density contrast is quadratic.\n";
		cout << "Computing gravity gradien tensor ..." << endl;
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			T[i].set();
			for (unsigned j = 0; j < num_tet; j++)
			{
				//clock_t t_start = clock();
				T[i] = T[i] + gra[i].tensor_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);
				T[i] = T[i] + gra[i].tensor_1st(*tetrahedrons[j], observation(i),
												(*tetrahedrons[j]).lin[0], (*tetrahedrons[j]).lin[1], (*tetrahedrons[j]).lin[2]);
				T[i] += gra[i].tensor_2nd(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).qua[0],
										  (*tetrahedrons[j]).qua[1], (*tetrahedrons[j]).qua[2], (*tetrahedrons[j]).qua[3],
										  (*tetrahedrons[j]).qua[4], (*tetrahedrons[j]).qua[5]);
			}
		}
	}
	else if (density_order == 3)
	{
		cout << "\nThe density contrast is cubic.\n";
		cout << "Computing gravity gradien tensor ..." << endl;
#pragma omp parallel for
		for (int i = 0; i < num_site; i++)
		{
			T[i].set();
			for (unsigned j = 0; j < num_tet; j++)
			{
				T[i] = T[i] + gra[i].tensor_const(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).const_density);
				T[i] = T[i] + gra[i].tensor_1st(*tetrahedrons[j], observation(i),
												(*tetrahedrons[j]).lin[0], (*tetrahedrons[j]).lin[1], (*tetrahedrons[j]).lin[2]);
				T[i] += gra[i].tensor_2nd(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).qua[0],
										  (*tetrahedrons[j]).qua[1], (*tetrahedrons[j]).qua[2], (*tetrahedrons[j]).qua[3],
										  (*tetrahedrons[j]).qua[4], (*tetrahedrons[j]).qua[5]);
				T[i] += gra[i].tensor_3rd(*tetrahedrons[j], observation(i), (*tetrahedrons[j]).cub);
			}
		}
	}
	ti.stop();
	cout << "Computing gravity gradient tensors done!\n";
	cout << " Computaion time for GGTs is " << ti.getElapsedTimeInSec() << " second(s)\n\n";
}

void tet_model::read_config(string config_file)
{
	ifstream in_stream(config_file.c_str());
	string model_mesh;

	std::string line;
	while (std::getline(in_stream, line))
	{
		line_process(line, "#");
		if (line.empty())
			continue;
		std::istringstream iss(line);
		iss >> model_mesh;
		break;
	}
	std::cout << model_mesh << "\n";

	string station_file;
	while (std::getline(in_stream, line))
	{
		line_process(line, "#");
		if (line.empty())
			continue;
		std::istringstream iss(line);
		iss >> station_file;
		break;
	}
	std::cout << station_file << "\n";

	while (std::getline(in_stream, line))
	{
		line_process(line, "#");
		if (line.empty())
			continue;
		std::istringstream iss(line);
		iss >> density_order;
		break;
	}
	std::cout << density_order << "\n";

	unsigned int n_regions;
	while (std::getline(in_stream, line))
	{
		line_process(line, "#");
		if (line.empty())
			continue;
		std::istringstream iss(line);
		iss >> n_regions;
		break;
	}

	std::cout << n_regions << "\n";
	for (unsigned i = 0; i < n_regions; i++)
	{
		std::vector<double> temp;
		int marker;
		while (std::getline(in_stream, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			iss >> marker;
			break;
		}
		std::cout << marker << "\n";

		double rho_const;
		while (std::getline(in_stream, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			iss >> rho_const;
			break;
		}
		std::cout << "\t" << rho_const << "\n";
		region_table_0[marker] = rho_const;

		temp.resize(3);
		cout << "\t";
		while (std::getline(in_stream, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			for (int j = 0; j < 3; j++)
			{
				iss >> temp[j];
				std::cout << temp[j] << "\t";
			}
			break;
		}
		std::cout << "\n";
		region_table_1[marker] = temp;

		temp.resize(6);
		cout << "\t";
		while (std::getline(in_stream, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			for (int j = 0; j < 6; j++)
			{
				iss >> temp[j];
				std::cout << temp[j] << "\t";
			}
			break;
		}
		std::cout << "\n";
		region_table_2[marker] = temp;

		temp.resize(10);
		cout << "\t";
		while (std::getline(in_stream, line))
		{
			line_process(line, "#");
			if (line.empty())
				continue;
			std::istringstream iss(line);
			for (int j = 0; j < 10; j++)
			{
				iss >> temp[j];
				std::cout << temp[j] << "\t";
			}
			break;
		}
		std::cout << "\n";
		region_table_3[marker] = temp;
	}

	mesh.read_mesh(model_mesh);
	this->set_density();

	observation.read_site(station_file);
	int n_data = observation.get_n_obs();
	this->g.resize(n_data);
	this->T.resize(n_data);
}

void tet_model::set_density()
{
	typedef std::map<int, double>::iterator IT0;
	typedef std::map<int, std::vector<double>>::iterator IT;
	const vector<Tetrahedron *> &tetrahedrons = mesh.get_tets();
	for (int i = 0; i < mesh.get_n_tets(); i++)
	{
		int marker;
		marker = (*tetrahedrons[i]).get_marker();
		IT0 it0 = region_table_0.find(marker);
		IT it1 = region_table_1.find(marker);
		IT it2 = region_table_2.find(marker);
		IT it3 = region_table_3.find(marker);
		mesh.set_tet_density_0(i, (*it0).second);
		mesh.set_tet_density_1(i, (*it1).second);
		mesh.set_tet_density_2(i, (*it2).second);
		mesh.set_tet_density_3(i, (*it3).second);
	}
}

void tet_model::out_g(const char *name)
{
	unsigned num_site = observation.get_n_obs();
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
		os << setw(25) << left << observation(i).x
		   << setw(25) << left << observation(i).y
		   << setw(25) << left << observation(i).z;
		os << setprecision(14);
		os << scientific;
		os << setw(25) << left << g[i].x;
		os << setw(25) << left << g[i].y;
		os << setw(25) << left << g[i].z << '\n';
	}
}

void tet_model::out_T(const char *name)
{

	unsigned num_site = observation.get_n_obs();
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
		os << setw(25) << left << observation(i).x
		   << setw(25) << left << observation(i).y
		   << setw(25) << left << observation(i).z;
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