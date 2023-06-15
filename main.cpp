#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include "Mesh/MeshData.h"
#include "ProjectionCones.h"
#include "Prune.h"
#include "Opt.h"
#include "Simplifier/GCLF-Simplification/interface.h"
#include "QEMSimplification.h"


void printHelp(const std::string& programName)
{
	std::cout << "Help: "
		<< programName
		<< " input_obj_file"
		<< " output_file"
		<<" [--normBound=sigma] "
		<< std::endl;
}

bool doesArgExist(const std::string& arg, const std::string& searchStr)
{
	return arg.find(searchStr) != std::string::npos;
}

bool parseArg(const std::string& arg, const std::string& searchStr, std::string& value)
{
	if (doesArgExist(arg, searchStr)) {
		value = arg.substr(arg.find_first_of(searchStr[searchStr.size() - 1]) + 1);
		return true;
	}
	return false;
}

void parseArgs(int argc, const char* argv[], std::string& objPath, std::string& outpath, double& sigma)
{
	if (argc < 3) {
		printHelp(argv[0]);
		exit(EXIT_FAILURE);
	}
	else {
		objPath = argv[1];
		outpath = argv[2];

		std::string tmp;
		for (int i = 3; i < argc; i++) {
			if (parseArg(argv[i], "--normBound=", tmp)) sigma = std::stod(tmp);
		}
	}
}

void saveCones(const VectorXd& conesK, std::string conesPath, double eps = 1e-3)
{
	std::ofstream conesFile(conesPath);
	if (conesFile.fail())
	{
		std::cout << "Open " << conesPath << "failed\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < conesK.size(); ++i)
	{
		if (conesK[i] > -eps && conesK[i] < eps) continue;
		conesFile << i + 1 << " " << conesK[i] * 0.5 * M_PI << std::endl;
	}
	conesFile.close();
}

void saveU(const VectorXd& u, std::string uPath)
{
	std::ofstream uFile(uPath);
	if (uFile.fail())
	{
		std::cout << "Open " << uPath << "failed\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < u.size(); ++i)
	{
		uFile << u[i] << std::endl;
	}
	uFile.close();
}

void saveINFO(double L2norm, double parameter1, double parameter2,double parameter3,double parameter4,double parameter5, std::string infoPath)
{
	std::ofstream infoFile(infoPath);
	if (infoFile.fail())
	{
		std::cout << "Open " << infoPath << "failed\n";
		exit(EXIT_FAILURE);
	}

	infoFile << "distortion : \t" << L2norm << std::endl;
	//infoFile << "otherInfo : \t" << parameter1 <<" "<< parameter3 <<" " << parameter4 <<" "<< parameter5 << std::endl;
	infoFile << "cost time : \t" << parameter2 << " s\n";
	infoFile.close();
}

int main(int argc, const char* argv[])
{
	std::string objPath = "";
	std::string outPath = "";
	double sigma = 0.2;// default distortion bound

	parseArgs(argc, argv, objPath, outPath, sigma);

	// load mesh
	Mesh mesh;
	Mesh mesh_ori;
	std::cout << "load mesh from " << objPath << std::endl;
	if (!MeshTools::ReadMesh(mesh, objPath))
	{
		std::cout << "load failed!\n";
		exit(EXIT_FAILURE);
	}
	
	double ori_sigma = sigma;
	double ori_sigma_bound = 0.8*sigma;

	// run our method
	if (mesh.n_vertices() > 5000)
	{	
		clock_t start, end_smiplify, end,end_prune;
		start = clock();
		std::cout << "QEM: mesh vertex number above 5000!" << endl;

		QEMSimplification qs(mesh);
		qs.Simplify(5000, 1);
		Mesh mesh_simplify0;
		MeshTools::Reassign(qs.GetMesh(), mesh_simplify0);
		OpenMesh::IO::write_mesh(mesh_simplify0, ".//simplify_mesh_qem.obj");

		int mesh_simlify_num = mesh_simplify0.n_vertices();
		std::cout << "origin mesh vertex number: " << mesh.n_vertices() << " simplified mesh vertex number: " << mesh_simplify0.n_vertices() << endl;
		end_smiplify = clock();
		double smiplify_time = (double)(end_smiplify - start) / CLOCKS_PER_SEC;
		std::cout << "************ IRLS Start ************" << endl;
		Prune pr0(mesh_simplify0);
		int try_wrong_num_geq_m3 = 0;
		int try_wrong_num_leq_m3 = 0;
		flag2:
		try {
			pr0.socp(mesh_simplify0.n_vertices(), 0.8 * sigma,20);
		}
		catch (...)
		{
			int signal_=round(std::accumulate(pr0.socp_rho.begin(), pr0.socp_rho.end(), 0.0));
			if (pr0.test_socp_num > 0&& signal_==8)
			{ 
				std::cout << "***** mode 2 *****" << endl;
				try {
					pr0.round();
				}
				catch (...)
				{
					Prune pr1(mesh_simplify0);
					int try_wrong_num_geq = 0;
					int try_wrong_num_leq = 0;
					flag_1:
					try {
						pr1.socp(mesh_simplify0.n_vertices(), 0.8*sigma, 20);
					}
					catch (...)
					{
						std::cout << "socp out" << endl;
					}

					try {
						pr1.round();
					}
					catch (...)
					{
						if (sigma > 0.1)
						{
							if (try_wrong_num_leq==0)
							{
								sigma = sigma - 0.01;
								pr1.prune_clear();
								try_wrong_num_geq++;
								goto flag_1;
							}
							else
							{
								sigma = sigma + 0.01;
								pr1.prune_clear();
								try_wrong_num_leq++;
								goto flag_1;
							}
						}
						else
						{
							if (try_wrong_num_geq !=0)
							{
								sigma = ori_sigma + 0.01;
								ori_sigma = sigma;
								try_wrong_num_leq++;
								pr1.prune_clear();
								goto flag_1;
							}
							else
							{
								sigma = sigma + 0.01;
								pr1.prune_clear();
								try_wrong_num_leq++;
								goto flag_1;
							}
						}
					}

					end_prune = clock();
					double prune_time = (double)(end_prune - start) / CLOCKS_PER_SEC;

					Opt intcones1(mesh_simplify0);
					intcones1.sigma_bound = 0.8 * sigma;
					intcones1.search_rings_num = 10;

					intcones1.calc_cones_position(pr1.socp_rho, pr1.round_rho, prune_time);
					std::cout << std::endl;
					std::cout << "projection simplify to ori and opt" << endl;
					Opt intcones_ori1(mesh);
					ProjectionCones pro_cones1;
					pro_cones1.projection_onto_ori_mesh(mesh_simplify0, mesh, intcones1.vinfo);
					for (int i = 0; i < intcones1.vinfo.size(); i++)
					{
						intcones1.vinfo[i].idx_i = pro_cones1.after_projection[i];
					}
					intcones_ori1.ori_mesh_local_search(intcones1.vinfo, 3);

					std::cout << endl;
					end = clock();
					std::cout << "total time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
					double total_time = (double)(end - start) / CLOCKS_PER_SEC;
					// output
					std::string infoPath = outPath + "-info.txt";
					std::string uPath = outPath + "-u.txt";
					std::string conesPath = outPath + "-cones.txt";
					saveINFO(intcones_ori1.ori_mesh_distortion, intcones1.search_merge_info[0], total_time, pr1.interval_round_time, pr1.total_round_time, total_time- prune_time, infoPath);
					saveU(intcones_ori1.final_U, uPath);
					saveCones(intcones_ori1.rho_vec, conesPath);
					return 0;
				}
				end_prune= clock();
				double prune_time= (double)(end_prune - start) / CLOCKS_PER_SEC;

				Opt intcones0(mesh_simplify0);
				intcones0.sigma_bound = 0.8 * sigma;
				intcones0.search_rings_num = 10;
				intcones0.calc_cones_position(pr0.socp_rho, pr0.round_rho, prune_time);
				std::cout << std::endl;
				std::cout << "projection simplify to ori and opt" << endl;
				Opt intcones_ori0(mesh);
				ProjectionCones pro_cones0;
				pro_cones0.projection_onto_ori_mesh(mesh_simplify0, mesh, intcones0.vinfo);
				for (int i = 0; i < intcones0.vinfo.size(); i++)
				{
					intcones0.vinfo[i].idx_i = pro_cones0.after_projection[i];
				}
				intcones_ori0.ori_mesh_local_search(intcones0.vinfo, 3);

				std::cout << endl;
				end = clock();
				std::cout << "total time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
				double total_time = (double)(end - start) / CLOCKS_PER_SEC;
				// output
				std::string infoPath = outPath + "-info.txt";
				std::string uPath = outPath + "-u.txt";
				std::string conesPath = outPath + "-cones.txt";
				saveINFO(intcones_ori0.ori_mesh_distortion, intcones0.search_merge_info[0], total_time, pr0.interval_round_time, pr0.total_round_time, total_time - prune_time, infoPath);
				saveU(intcones_ori0.final_U, uPath);
				saveCones(intcones_ori0.rho_vec, conesPath);
				return 0;
			}
			else
			{
				std::cout << "***** mode 3 *****" << endl;
				clock_t start, end_smiplify, end;
				start = clock();
				std::cout << "NonQEM: mesh vertex number above 5000!" << endl;

				GCLF::Simplification::SMeshT smesh_ori, smesh_simlify;
				GCLF::Simplification::SMeshT smesh;
				Mesh mesh_simlify;
				convertMeshTo(mesh, smesh);
				GCLF::Simplification::simplify(smesh, smesh_simlify, 5000);
				convertMeshFrom(smesh_simlify, mesh_simlify);
				OpenMesh::IO::write_mesh(mesh_simlify, ".//simplify_mesh_noqem.obj");
				int mesh_simlify_num = mesh_simlify.n_vertices();
				std::cout << "origin mesh vertex number: " << mesh.n_vertices() << " simplified mesh vertex number: " << mesh_simlify.n_vertices() << endl;
				end_smiplify = clock();
				double smiplify_time = (double)(end_smiplify - start) / CLOCKS_PER_SEC;

				Prune pr(mesh_simlify);

				int try_wrong_num_geq = 0;
				int try_wrong_num_leq = 0;
				flag_0:
				try {
					pr.socp(mesh_simlify.n_vertices(), 0.8 * sigma, 20);
				}
				catch (...)
				{
					std::cout << "socp exit ahead!" << endl;

				}

				try {
					pr.round();
				}
				catch (...)
				{
					if (sigma > 0.1)
					{
						if (try_wrong_num_leq == 0)
						{
							sigma = sigma - 0.01;
							pr.prune_clear();
							try_wrong_num_geq++;
							goto flag_0;
						}
						else
						{
							sigma = sigma + 0.01;
							pr.prune_clear();
							try_wrong_num_leq++;
							goto flag_0;
						}
					}
					else
					{
						if (try_wrong_num_geq != 0)
						{
							sigma = ori_sigma + 0.01;
							ori_sigma = sigma;
							try_wrong_num_leq++;
							pr.prune_clear();
							goto flag_0;
						}
						else
						{
							sigma = sigma + 0.01;
							pr.prune_clear();
							try_wrong_num_leq++;
							goto flag_0;
						}
					}
				}

				end_prune = clock();
				double prune_time = (double)(end_prune - start) / CLOCKS_PER_SEC;

				Opt intcones(mesh_simlify);
				intcones.sigma_bound = 0.8 * sigma;
				intcones.search_rings_num = 10;
				intcones.calc_cones_position(pr.socp_rho, pr.round_rho, prune_time);
				std::cout << std::endl;
				std::cout << "projection simplify to ori and opt" << endl;
				Opt intcones_ori(mesh);
				ProjectionCones pro_cones;
				pro_cones.projection_onto_ori_mesh(mesh_simlify, mesh, intcones.vinfo);
				for (int i = 0; i < intcones.vinfo.size(); i++)
				{
					intcones.vinfo[i].idx_i = pro_cones.after_projection[i];
				}
				intcones_ori.ori_mesh_local_search(intcones.vinfo, 3);

				std::cout << endl;
				end = clock();
				std::cout << "total time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
				double total_time = (double)(end - start) / CLOCKS_PER_SEC;
				// output
				std::string infoPath = outPath + "-info.txt";
				std::string uPath = outPath + "-u.txt";
				std::string conesPath = outPath + "-cones.txt";
				saveINFO(intcones_ori.ori_mesh_distortion, intcones.search_merge_info[0], total_time, pr.interval_round_time, pr.total_round_time, total_time - prune_time, infoPath);
				saveU(intcones_ori.final_U, uPath);
				saveCones(intcones_ori.rho_vec, conesPath);
				return 0;
			}
		}

		std::cout << "***** mode 1 *****" << endl;
		try {
			pr0.round();
		}
		catch (...)
		{
			if (sigma > 0.1)
			{
				if (try_wrong_num_leq_m3 == 0)
				{
					sigma = sigma - 0.01;
					pr0.prune_clear();
					try_wrong_num_geq_m3++;
					goto flag2;
				}
				else
				{
					sigma = sigma + 0.01;
					pr0.prune_clear();
					try_wrong_num_leq_m3++;
					goto flag2;
				}
			}
			else
			{
				if (try_wrong_num_geq_m3 != 0)
				{
					sigma = ori_sigma + 0.01;
					ori_sigma = sigma;
					try_wrong_num_leq_m3++;
					pr0.prune_clear();
					goto flag2;
				}
				else
				{
					sigma = sigma + 0.01;
					pr0.prune_clear();
					try_wrong_num_leq_m3++;
					goto flag2;
				}
			}
		}

		end_prune = clock();
		double prune_time = (double)(end_prune - start) / CLOCKS_PER_SEC;

		Opt intcones0(mesh_simplify0);
		intcones0.sigma_bound = 0.8 * sigma;
		intcones0.search_rings_num =10;
		intcones0.calc_cones_position(pr0.socp_rho, pr0.round_rho, prune_time);

		std::cout << std::endl;
		std::cout << "projection simplify to ori and opt" << endl;
		Opt intcones_ori0(mesh);
		ProjectionCones pro_cones0;
		pro_cones0.projection_onto_ori_mesh(mesh_simplify0, mesh, intcones0.vinfo);
		for (int i = 0; i < intcones0.vinfo.size(); i++)
		{
			intcones0.vinfo[i].idx_i = pro_cones0.after_projection[i];
		}
		intcones_ori0.ori_mesh_local_search(intcones0.vinfo, 3);

		std::cout << endl;
		end = clock();
		std::cout << "total time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
		double total_time = (double)(end - start) / CLOCKS_PER_SEC;
		// output
		std::string infoPath = outPath + "-info.txt";
		std::string uPath = outPath + "-u.txt";
		std::string conesPath = outPath + "-cones.txt";
		saveINFO(intcones_ori0.ori_mesh_distortion, intcones0.search_merge_info[0], total_time, pr0.interval_round_time, pr0.total_round_time, total_time - prune_time, infoPath);
		saveU(intcones_ori0.final_U, uPath);
		saveCones(intcones_ori0.rho_vec, conesPath);
	}
	else
	{
		std::cout << "************ IRLS Start ************" << endl;
		clock_t start, end, end_prune;
		start = clock();
		std::cout << "mesh vertex number less 5000!" << endl;
		Prune pr(mesh);
		int try_wrong_num_geq_m5= 0;
		int try_wrong_num_leq_m5 = 0;
		flag_5:
		try {
			pr.socp(mesh.n_vertices(), 0.8 * sigma,20);
		}
		catch (...)
		{
			std::cout << "socp exit ahead!" << endl;
			int signal_1 = round(std::accumulate(pr.socp_rho.begin(), pr.socp_rho.end(), 0.0));
			if (pr.test_socp_num > 0 && signal_1 == 8)
			{
				std::cout << "***** model 4 ******" << endl;
				int try_wrong_num_geq = 0;
				int try_wrong_num_leq = 0;
				flag_4:
				try{
					pr.socp(mesh.n_vertices(), 0.8*sigma, 20);
				}
				catch (...)
				{
					std::cout << "socp exit ahead again!" << endl;
				}

				try {
					pr.round();
				}
				catch (...)
				{
					if (sigma > 0.1)
					{
						if (try_wrong_num_leq == 0)
						{
							sigma = sigma - 0.01;
							pr.prune_clear();
							try_wrong_num_geq++;
							goto flag_4;
						}
						else
						{
							sigma = sigma + 0.01;
							pr.prune_clear();
							try_wrong_num_leq++;
							goto flag_4;
						}
					}
					else
					{
						if (try_wrong_num_geq != 0)
						{
							sigma = ori_sigma + 0.01;
							ori_sigma = sigma;
							try_wrong_num_leq++;
							pr.prune_clear();
							goto flag_4;
						}
						else
						{
							sigma = sigma + 0.01;
							pr.prune_clear();
							try_wrong_num_leq++;
							goto flag_4;
						}
					}
				}

				end_prune = clock();
				double prune_time = (double)(end_prune - start) / CLOCKS_PER_SEC;
				Opt intcones(mesh);
				intcones.sigma_bound = 0.8 * sigma;
				intcones.search_rings_num = 10;
				intcones.calc_cones_position(pr.socp_rho, pr.round_rho, prune_time);

				std::cout << endl;
				end = clock();
				std::cout << "total time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
				double total_time = (double)(end - start) / CLOCKS_PER_SEC;
				// output
				std::string infoPath = outPath + "-info.txt";
				std::string uPath = outPath + "-u.txt";
				std::string conesPath = outPath + "-cones.txt";
				saveINFO(intcones.search_merge_info[1], intcones.search_merge_info[0], total_time, pr.interval_round_time, pr.total_round_time, total_time - prune_time, infoPath);
				saveU(intcones.final_U, uPath);
				saveCones(intcones.rho_vec, conesPath);
				return 0;
			}
			else
			{
				pr.prune_clear();
				std::cout << "***** model 6 ******" << endl;
				Prune pr6(mesh);
				int try_wrong_num_geq_m6 = 0;
				int try_wrong_num_leq_m6 = 0;
				int try_wrong_num_geq_m6_ = 0;
				int try_wrong_num_leq_m6_ = 0;
				flag_6:
				try {
					pr6.socp(mesh.n_vertices(), 0.8 * sigma, 20);
				}
				catch (...)
				{
					if (sigma > 0.1)
					{
						if (try_wrong_num_geq_m6 == 0)
						{
							sigma = sigma - 0.01;
							pr6.prune_clear();
							try_wrong_num_geq_m6++;
							goto flag_6;
						}
						else
						{
							sigma = sigma + 0.01;
							pr6.prune_clear();
							try_wrong_num_leq_m6++;
							goto flag_6;
						}
					}
					else
					{
						if (try_wrong_num_geq_m6 != 0)
						{
							sigma = ori_sigma + 0.01;
							ori_sigma = sigma;
							try_wrong_num_leq_m6++;
							pr6.prune_clear();
							goto flag_6;
						}
						else
						{
							sigma = sigma + 0.01;
							pr6.prune_clear();
							try_wrong_num_leq_m6++;
							goto flag_6;
						}
					}
				}

				try {
					pr6.round();
				}
				catch (...)
				{
					if (sigma > 0.1)
					{
						if (try_wrong_num_leq_m6_ == 0)
						{
							sigma = sigma - 0.01;
							pr6.prune_clear();
							try_wrong_num_geq_m6_++;
							goto flag_6;
						}
						else
						{
							sigma = sigma + 0.01;
							pr6.prune_clear();
							try_wrong_num_leq_m6_++;
							goto flag_6;
						}
					}
					else
					{
						if (try_wrong_num_geq_m6_ != 0)
						{
							sigma = ori_sigma + 0.01;
							ori_sigma = sigma;
							try_wrong_num_leq_m6_++;
							pr6.prune_clear();
							goto flag_6;
						}
						else
						{
							sigma = sigma + 0.01;
							pr6.prune_clear();
							try_wrong_num_leq_m6_++;
							goto flag_6;
						}
					}
				}

				end_prune = clock();
				double prune_time = (double)(end_prune - start) / CLOCKS_PER_SEC;

				Opt intcones(mesh);
				intcones.sigma_bound = 0.8 * sigma;
				intcones.search_rings_num = 10;
				intcones.calc_cones_position(pr6.socp_rho, pr6.round_rho, prune_time);

				std::cout << endl;
				end = clock();
				std::cout << "total time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
				double total_time = (double)(end - start) / CLOCKS_PER_SEC;
				// output
				std::string infoPath = outPath + "-info.txt";
				std::string uPath = outPath + "-u.txt";
				std::string conesPath = outPath + "-cones.txt";
				saveINFO(intcones.search_merge_info[1], intcones.search_merge_info[0], total_time, pr.interval_round_time, pr.total_round_time, total_time - prune_time, infoPath);
				saveU(intcones.final_U, uPath);
				saveCones(intcones.rho_vec, conesPath);
				return 0;
			}
		}

		std::cout << "***** model 5 ******" << endl;
		try {
			pr.round();
		}
		catch (...)
		{ 
			if (sigma > 0.1)
			{
				if (try_wrong_num_leq_m5 == 0)
				{
					sigma = sigma - 0.01;
					pr.prune_clear();
					try_wrong_num_geq_m5++;
					goto flag_5;
				}
				else
				{
					sigma = sigma + 0.01;
					pr.prune_clear();
					try_wrong_num_leq_m5++;
					goto flag_5;
				}
			}
			else
			{
				if (try_wrong_num_geq_m5 != 0)
				{
					sigma = ori_sigma + 0.01;
					ori_sigma = sigma;
					try_wrong_num_leq_m5++;
					pr.prune_clear();
					goto flag_5;
				}
				else
				{
					sigma = sigma + 0.01;
					pr.prune_clear();
					try_wrong_num_leq_m5++;
					goto flag_5;
				}
			}
		}

		end_prune = clock();
		double prune_time = (double)(end_prune - start) / CLOCKS_PER_SEC;

		Opt intcones(mesh);
		intcones.sigma_bound = 0.8*sigma;
		intcones.search_rings_num = 10;
		intcones.calc_cones_position(pr.socp_rho, pr.round_rho, prune_time);

		std::cout << endl;
		end = clock();
		std::cout << "total time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
		double total_time = (double)(end - start) / CLOCKS_PER_SEC;
		// out info
		std::string infoPath = outPath + "-info.txt";
		std::string uPath = outPath + "-u.txt";
		std::string conesPath = outPath + "-cones.txt";
		saveINFO(intcones.search_merge_info[1], intcones.search_merge_info[0], total_time, pr.interval_round_time, pr.total_round_time, total_time - prune_time, infoPath);
		saveU(intcones.final_U, uPath);
		saveCones(intcones.rho_vec, conesPath);
	}
	return 0;
}

