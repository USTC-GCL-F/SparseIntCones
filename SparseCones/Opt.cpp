#include "Opt.h"

Opt::Opt(Mesh& mesh) :mesh(mesh)
{

}

void Opt::VertArea()
{
	vec_A.setConstant(mesh.n_vertices(), 0);
	double fA;
	for (const FH& f_h : mesh.faces())
	{
		OpenMesh::Vec3d p[3];
		auto cfv_it = mesh.cfv_begin(f_h);
		p[0] = mesh.point(*cfv_it);		++cfv_it;
		p[1] = mesh.point(*cfv_it);		++cfv_it;
		p[2] = mesh.point(*cfv_it);

		fA = ((p[0] - p[1]) % (p[1] - p[2])).norm() / 3.0;

		for (const VH& fv_h : mesh.fv_range(f_h))
		{
			vec_A[fv_h.idx()] += fA;
		}
	}

	vec_A = vec_A / vec_A.sum();
}

void Opt::FaceAngle()
{
	f_angle.setConstant(mesh.n_faces() * 3, 0);
	for (const FH& f_h : mesh.faces())
	{
		std::vector<OpenMesh::Vec3d> he_v(3);
		int cout = 0;
		for (const HEH& he_h : mesh.fh_range(f_h))
		{
			he_v[cout] = mesh.calc_edge_vector(he_h).normalized();
			cout++;
		}

		double cos_0 = -dot(he_v[1], he_v[2]);
		double cos_1 = -dot(he_v[2], he_v[0]);
		double cos_2 = -dot(he_v[0], he_v[1]);

		f_angle[3 * f_h.idx() + 0] = acos(cos_0);
		f_angle[3 * f_h.idx() + 1] = acos(cos_1);
		f_angle[3 * f_h.idx() + 2] = acos(cos_2);
	}
}

void Opt::VertGauss()
{
	K_ori.setConstant(mesh.n_vertices(), 2 * M_PI);

	for (const VH& v_h : mesh.vertices())
	{
		if (mesh.is_boundary(v_h)) {
			K_ori[v_h.idx()] = M_PI;
		}
	}

	for (const FH& f_h : mesh.faces())
	{
		std::vector<int> v_ids(3);
		int cout = 0;
		for (const HEH& he_h : mesh.fh_range(f_h))
		{
			v_ids[cout] = mesh.to_vertex_handle(mesh.next_halfedge_handle(he_h)).idx();
			cout++;
		}

		K_ori[v_ids[0]] -= f_angle[3 * f_h.idx() + 0];
		K_ori[v_ids[1]] -= f_angle[3 * f_h.idx() + 1];
		K_ori[v_ids[2]] -= f_angle[3 * f_h.idx() + 2];
	}

}


void Opt::MatrixL()
{
	double cot_weight;
	VH to_v, from_v;
	std::vector<T> trip;
	trip.reserve(12 * mesh.n_faces());
	for (const FH& f_h : mesh.faces())
	{
		int cout = 0;
		for (const HEH& fh_h : mesh.fh_range(f_h))
		{
			to_v = mesh.to_vertex_handle(fh_h);
			from_v = mesh.from_vertex_handle(fh_h);

			cot_weight = 0.5 / tan(f_angle[3 * f_h.idx() + cout]);
			cout++;

			trip.emplace_back(to_v.idx(), from_v.idx(), -cot_weight);
			trip.emplace_back(from_v.idx(), to_v.idx(), -cot_weight);
			trip.emplace_back(from_v.idx(), from_v.idx(), cot_weight);
			trip.emplace_back(to_v.idx(), to_v.idx(), cot_weight);
		}
	}

	m_L.resize(mesh.n_vertices(), mesh.n_vertices());
	m_L.setFromTriplets(trip.begin(), trip.end());
}

void Opt::Init()
{
	VertArea();
	FaceAngle();
	VertGauss();
	MatrixL();


	P.resize(mesh.n_vertices(), mesh.n_vertices() - 1);
	std::vector<T> element;
	for (int i = 0; i < mesh.n_vertices() - 1; i++)
	{
		element.emplace_back(i, i, 1);
	}

	P.setFromTriplets(element.begin(), element.end());
	PT = P.transpose();

}

void Opt::init_vinfo_noncone(vector<double>& a, vector<double>& b)
{
	non_cone.clear();
	for (int i = 0; i < a.size(); i++)
	{
		if (fabs(a[i]) > 1e-4)
		{
			VInfo vif;
			vif.after_reweighed_rho = b[i];
			vif.after_round = a[i];
			vif.idx_i = i;
			vif.gap = fabs(b[i] - a[i]);
			vinfo.push_back(vif);
		}
		else
		{
			non_cone.push_back(i);
		}
	}

	// Ensure that the first two are one positive and one negative
	std::sort(vinfo.begin(), vinfo.end(), com_smaller_gap);
	double first_element_rho = vinfo[0].after_round;
	VInfo position2th;
	for (int i = 1; i < vinfo.size(); i++)
	{
		double element_rho = vinfo[i].after_round;
		if (first_element_rho * element_rho < 0)
		{
			position2th.idx_i = vinfo[i].idx_i;
			position2th.is_cone = vinfo[i].is_cone;
			position2th.after_reweighed_rho = vinfo[i].after_reweighed_rho;
			position2th.after_round = vinfo[i].after_round;
			position2th.gap = vinfo[i].gap;
			vinfo.erase(vinfo.begin() + i);
			vinfo.insert(vinfo.begin() + 1, position2th);
			break;
		}
	}
}

void Opt::init_equatuion_data()
{
	Init();
	E.setOnes(mesh.n_vertices());
	solver.compute(PT * m_L * P);
	vector<pair<int, double>> compare;

	B = 0.5 * M_PI * rho_vec - K_ori;
	U = solver.solve(PT * B);
	Eigen::VectorXd U2;
	U2.resize(mesh.n_vertices());
	for (int k = 0; k < U.size(); k++)
	{
		U2[k] = U[k];
	}
	U2[mesh.n_vertices() - 1] = 0;

	double xishu1 = U2.transpose() * vec_A.asDiagonal() * E;
	double xishi2 = E.transpose() * vec_A.asDiagonal() * E;
	Eigen::VectorXd U1 = U2 - ((double)(xishu1 / xishi2)) * E;
	distortion = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
	distortion_last = distortion;
	std::cout << "Initial distortion of the optimization£º " << distortion << endl;
	final_U = U1;
	//cout <<"error: " <<(m_L * U1 - B).norm() << endl;
	search_merge_info.push_back(distortion);
}

double Opt::calc_mesh_distortion(Eigen::VectorXd rho_vec)
{
	B = 0.5 * M_PI * rho_vec - K_ori;
	U = solver.solve(PT * B);

	Eigen::VectorXd U2;
	U2.resize(mesh.n_vertices());
	for (int k = 0; k < U.size(); k++)
	{
		U2[k] = U[k];
	}
	U2[mesh.n_vertices() - 1] = 0;

	double xishu1 = U2.transpose() * vec_A.asDiagonal() * E;
	double xishi2 = E.transpose() * vec_A.asDiagonal() * E;
	Eigen::VectorXd U1 = U2 - ((double)(xishu1 / xishi2)) * E;
	double d = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
	final_U = U1;

	return d;
}

void Opt::global_search()
{
	for (int i = 0; i < non_cone.size(); i++)
	{
		rho_vec[non_cone[i]] = rj;
		distortion = calc_mesh_distortion(rho_vec);
		compare.push_back(make_pair(non_cone[i], distortion));
		rho_vec[non_cone[i]] = 0;
	}
}

void Opt::local_search(vector<int>& a)
{
	for (int i = 0; i < a.size(); i++)
	{
		rho_vec[a[i]] = rj;
		double distortion = calc_mesh_distortion(rho_vec);
		compare.push_back(make_pair(a[i], distortion));
		rho_vec[a[i]] = 0;
	}
}

void Opt::calc_cones_position(vector<double> socp_res, vector<double> round_res,double p_time)
{
	std::cout << std::endl;
	std::cout << "************** The Optimization Procedure Start **************" << endl;
	// 1.load data
	clock_t startTime, endTime,end_opt;
	startTime = clock();
	vector<double> vec_rho = round_res;
	vector<double> after_reweighed_rho = socp_res;

	// 2.init vinfo and non_cone
	init_vinfo_noncone(vec_rho, after_reweighed_rho);

	// 3.init data and calc init distortion ( after around ) 
	rho_vec = Eigen::Map<VectorXd>(vec_rho.data(), mesh.n_vertices());
	init_equatuion_data();

	if (p_time>1800)
	{
		search_merge_info.push_back(distortion);
		int cones_number = 0;
		for (int i = 0; i < rho_vec.size(); i++)
		{
			if (rho_vec[i] != 0)
			{
				cones_number++;
			}
		}
		search_merge_info.push_back(cones_number);
		goto flag1;
	}

	// 4.loop {global, local(break), merge}_{loop=3}
	for (int y = 0; y < loop_mix_num; y++)
	{
		int global_iter = 1;
		double delta_distortion = 1;
		distortion_last = calc_mesh_distortion(rho_vec);
		int decrease_num = 1;
		while (delta_distortion > 1e-2 && decrease_num <= 3)
		{
			global_search_num = 0;
			for (int j = 0; j < vinfo.size(); j++)
			{
				if (global_iter <= 2)
				{
					rho_vec[vinfo[j].idx_i] = 0;
					//std::cout << "before moving idx: " << vinfo[j].idx_i << "  **  ";
					non_cone.push_back(vinfo[j].idx_i);
					rj = vinfo[j].after_round;

					global_search();

					std::sort(compare.begin(), compare.end(), cmp_bigger);
					rho_vec[compare[0].first] = rj;

					vector<int>::iterator it = find(non_cone.begin(), non_cone.end(), compare[0].first);
					non_cone.erase(it);
					non_cone.push_back(vinfo[j].idx_i);

					vinfo[j].idx_i = compare[0].first;//update cone idx

					//std::cout << " global search  ***  after moving idx: " << compare[0].first << " " << "distortion: " << compare[0].second << "  round value: " << " " << rj << endl;

					distortion = compare[0].second;

					compare.clear();
					non_cone.pop_back();

					global_iter++;
					global_search_num++;

					end_opt = clock();
					double opt_time = (double)(end_opt - startTime) / CLOCKS_PER_SEC;
					if (opt_time+p_time>1800)
					{
						goto flag2;
					}
				}
				else
				{
					if (distortion > ratio_sigma * sigma_bound) //ratio_sigma*0.2
					{
						rho_vec[vinfo[j].idx_i] = 0;
						//std::cout << "before moving idx: " << vinfo[j].idx_i << "  **  ";
						non_cone.push_back(vinfo[j].idx_i);
						rj = vinfo[j].after_round;

						global_search();

						std::sort(compare.begin(), compare.end(), cmp_bigger);
						rho_vec[compare[0].first] = rj;

						vector<int>::iterator it = find(non_cone.begin(), non_cone.end(), compare[0].first);
						non_cone.erase(it);
						non_cone.push_back(vinfo[j].idx_i);

						vinfo[j].idx_i = compare[0].first;//update cone idx

						//std::cout << " global search  ***  after moving idx: " << compare[0].first << " " << "distortion: " << compare[0].second << "  round value: " << " " << rj << endl;

						distortion = compare[0].second;

						compare.clear();
						non_cone.pop_back();
						global_search_num++;

						end_opt = clock();
						double opt_time = (double)(end_opt - startTime) / CLOCKS_PER_SEC;
						if (opt_time + p_time > 1800)
						{
							goto flag2;
						}
					}
					else
					{
						set<int> n_rings_set;
						vector<int> n_rings_vec;
						n_rings_set.clear();
						n_rings_vec.clear();
						auto v = mesh.vertex_handle(vinfo[j].idx_i);
						n_rings_set.insert(vinfo[j].idx_i);

						for (auto vv : mesh.vv_range(v))
						{
							n_rings_set.insert(vv.idx());
							for (auto vvv : mesh.vv_range(vv))
							{
								n_rings_set.insert(vvv.idx());
							}
						}

						n_rings_vec.assign(n_rings_set.begin(), n_rings_set.end());

						for (int k = 0; k < vinfo.size(); k++)
						{
							vector<int>::iterator it = find(n_rings_vec.begin(), n_rings_vec.end(), vinfo[k].idx_i);
							if (it != n_rings_vec.end())
							{
								n_rings_vec.erase(it);
							}
						}

						rho_vec[vinfo[j].idx_i] = 0;
						//std::cout << "before moving idx: " << vinfo[j].idx_i;
						n_rings_vec.push_back(vinfo[j].idx_i);
						rj = vinfo[j].after_round;

						local_search(n_rings_vec);

						std::sort(compare.begin(), compare.end(), cmp_bigger);
						rho_vec[compare[0].first] = rj;
						vinfo[j].idx_i = compare[0].first;

						//std::cout << " ** local search  ***  after moving idx: " << compare[0].first << " " << " distortion: " << compare[0].second << "  round value: " << " " << rj << endl;

						compare.clear();
						n_rings_set.clear();
						n_rings_vec.clear();

						end_opt = clock();
						double opt_time = (double)(end_opt - startTime) / CLOCKS_PER_SEC;
						if (opt_time + p_time > 1800)
						{
							goto flag2;
						}
					}
				}
			}
			double distortion_new = calc_mesh_distortion(rho_vec);
			delta_distortion = distortion_last - distortion_new;
			distortion_last = distortion_new;
	
			vector<VInfo> replace_vinfo;
			replace_vinfo.resize(vinfo.size());
			std::copy(vinfo.begin() + global_search_num, vinfo.end(), replace_vinfo.begin());
			std::copy(vinfo.begin(), vinfo.begin() + global_search_num, replace_vinfo.begin() + vinfo.size() - global_search_num);
			vinfo.clear();
			vinfo = replace_vinfo;
			replace_vinfo.clear();
			decrease_num++;
		}

		if (y == loop_mix_num - 1)
		{
			flag2:
			double stop_dist = calc_mesh_distortion(rho_vec);
			search_merge_info.push_back(stop_dist);
			int cones_number = 0;
			for (int i = 0; i < rho_vec.size(); i++)
			{
				if (rho_vec[i] != 0)
				{
					cones_number++;
				}
			}
			std::cout << "Iters: " << y + 1 << " distortion: " << stop_dist << " cone number: " << cones_number << endl;
			search_merge_info.push_back(cones_number);
			break;
		}

		merge_cones_pair(search_rings_num);

		double final_distortion = calc_mesh_distortion(rho_vec);

		int c_num = 0;
		for (int i = 0; i < mesh.n_vertices(); i++)
		{
			double a = rho_vec[i];
			if (abs(a) > 1e-8)
			{
				c_num++;
			}
		}
		std::cout << "Iters: " << y + 1 << " distortion: "<< final_distortion <<" cone number: " << c_num << endl;
	}
	flag1:
	
	endTime = clock();
	std::cout << "************** The Optimization Procedure End **************" << endl;
	std::cout << "The optimization process time: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	search_merge_info.push_back((double)(endTime - startTime) / CLOCKS_PER_SEC);
}

void Opt::ori_mesh_local_search(vector<VInfo> vinfo, int loop_num)
{
	clock_t start, end;
	start = clock();
	Init();

	E.setOnes(mesh.n_vertices());
	rho_vec.setZero(mesh.n_vertices());
	for (int i = 0; i < vinfo.size(); i++)
	{
		rho_vec[vinfo[i].idx_i] = vinfo[i].after_round;
	}
	solver.compute(PT * m_L * P);

	B = 0.5 * M_PI * rho_vec - K_ori;
	U = solver.solve(PT * B);

	Eigen::VectorXd U2;
	U2.resize(mesh.n_vertices());
	for (int k = 0; k < U.size(); k++)
	{
		U2[k] = U[k];
	}
	U2[mesh.n_vertices() - 1] = 0;

	double xishu1 = U2.transpose() * vec_A.asDiagonal() * E;
	double xishi2 = E.transpose() * vec_A.asDiagonal() * E;
	Eigen::VectorXd U1 = U2 - ((double)(xishu1 / xishi2)) * E;
	distortion = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
	std::cout << "Original mesh initial distortion£º " << distortion << endl;
	ori_mesh_local_info.clear();
	ori_mesh_local_info.push_back(distortion);
	ori_mesh_local_info.push_back(vinfo.size());
	double final_local_search_distortion = 0;

	for (int m = 0; m < loop_num; m++)
	{
		std::cout << "Iter: " << m + 1 << " local search on original mesh" << endl;
		for (int j = 0; j < vinfo.size(); j++)
		{
			set<int> n_rings_set;
			vector<int> n_rings_vec;
			n_rings_set.clear();
			n_rings_vec.clear();
			auto v = mesh.vertex_handle(vinfo[j].idx_i);
			n_rings_set.insert(vinfo[j].idx_i);
			
			for (auto vv : mesh.vv_range(v))
			{
				n_rings_set.insert(vv.idx());
				for (auto vvv : mesh.vv_range(vv))
				{
					n_rings_set.insert(vvv.idx());
				}
			}

			n_rings_vec.assign(n_rings_set.begin(), n_rings_set.end());

			for (int k = 0; k < vinfo.size(); k++)
			{
				vector<int>::iterator it = find(n_rings_vec.begin(), n_rings_vec.end(), vinfo[k].idx_i);
				if (it != n_rings_vec.end())
				{
					n_rings_vec.erase(it);
				}
			}

			rho_vec[vinfo[j].idx_i] = 0;
			//std::cout << "before moving idx: " << vinfo[j].idx_i;
			n_rings_vec.push_back(vinfo[j].idx_i);
			rj = vinfo[j].after_round;

			local_search(n_rings_vec);

			std::sort(compare.begin(), compare.end(), cmp_bigger);
			rho_vec[compare[0].first] = rj;
			vinfo[j].idx_i = compare[0].first;//update cone placement

			//std::cout << " ** After projection original mesh ** local search  ***  after moving idx: " << compare[0].first << " " << " distortion: " << compare[0].second << "  round value: " << " " << rj << endl;
			final_local_search_distortion = compare[0].second;
			compare.clear();
			n_rings_set.clear();
			n_rings_vec.clear();
		}
	}
	end = clock();
	ori_mesh_local_info.push_back((double)(end - start) / CLOCKS_PER_SEC);
	ori_mesh_local_info.push_back(final_local_search_distortion);

	double ori_dist = calc_mesh_distortion(rho_vec);
	ori_mesh_distortion = ori_dist;
	std::cout << "***** original mesh *****"<<" distortion: "<< ori_dist <<endl;
	std::cout << "The time of original mesh optimization: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
}


double Opt::local_search2(VectorXd merge_rho, vector<VInfo> vv_info)
{
	vector<pair<int, double>> compare1;
	double rj_0 = 0.0;
	for (int j = 0; j < vv_info.size(); j++)
	{
		set<int> n_rings_set;
		vector<int> n_rings_vec;
		n_rings_set.clear();
		n_rings_vec.clear();
		auto v = mesh.vertex_handle(vv_info[j].idx_i);
		n_rings_set.insert(vv_info[j].idx_i);
		for (auto vv : mesh.vv_range(v))
		{
			n_rings_set.insert(vv.idx());
			for (auto vvv : mesh.vv_range(vv))
			{
				n_rings_set.insert(vvv.idx());
			}
		}

		n_rings_vec.assign(n_rings_set.begin(), n_rings_set.end());

		for (int k = 0; k < vv_info.size(); k++)
		{
			vector<int>::iterator it = find(n_rings_vec.begin(), n_rings_vec.end(), vv_info[k].idx_i);
			if (it != n_rings_vec.end())
			{
				n_rings_vec.erase(it);
			}
		}

		merge_rho[vv_info[j].idx_i] = 0;
		//std::cout << "before moving idx: " << vv_info[j].idx_i;
		n_rings_vec.push_back(vv_info[j].idx_i);
		rj_0 = vv_info[j].after_round;


		for (int i = 0; i < n_rings_vec.size(); i++)
		{
			merge_rho[n_rings_vec[i]] = rj_0;
			double distortion = calc_mesh_distortion(merge_rho);

			compare1.push_back(make_pair(n_rings_vec[i], distortion));
			//std::cout << "idx: " << compare[i].first << " " << "distortion: " << compare[i].second << endl;
			merge_rho[n_rings_vec[i]] = 0;
		}

		std::sort(compare1.begin(), compare1.end(), cmp_bigger);
		merge_rho[compare1[0].first] = rj_0;
		vv_info[j].idx_i = compare1[0].first;

		//std::cout << " ** local search  ***  after moving idx: " << compare1[0].first << " " << " distortion: " << compare1[0].second << "  round value: " << " " << rj_0 << endl;

		compare1.clear();
		n_rings_set.clear();
		n_rings_vec.clear();
	}
	vec_vec_vinfo.push_back(vv_info);
	vv_info.clear();
	compare1.clear();
	double dist_merge = calc_mesh_distortion(merge_rho);
	//cout << dist_merge << endl;
	return dist_merge;
}

void Opt::merge_cones_pair(int rings_num)
{
	//cout << "***** merge start *****" << endl;
	/*   update many variables
	//	 1): each update vinfo
	//   2): update rho_vec
	//   3): update non_cone; used for search
	//   4): posi_cone¡¢unposi_cone
	//   5): posi_cone_status unposi_cone_status
	*/

	// 1.initial
	clock_t start_mergy_time,now_mergy_time;
	start_mergy_time = clock();
	
	global_search_num = 0;

	vector<int> posi_cone;
	vector<int> unposi_cone;
	
	for (int i = 0; i < vinfo.size(); i++)
	{
		if (vinfo[i].after_round > 0)
		{
			posi_cone.push_back(vinfo[i].idx_i);
		}
		else
		{
			unposi_cone.push_back(vinfo[i].idx_i);
		}
	}
	vector<bool> unposi_cone_status;
	vector<bool> posi_cone_status;
	posi_cone_status.resize(posi_cone.size(), false);
	unposi_cone_status.resize(unposi_cone.size(), false);

	MeshCache MC(mesh);
	std::vector<int> path;
	std::vector<int> targets(MC.NVertices, 0);

	// 2.merge cone pair
	flag0:
	for (int i = 0; i < posi_cone.size(); i++)
	{
		vector<int> return_cones_in_vinfo_pos;
		return_cones_in_vinfo_pos.clear();
		return_cones_in_vinfo_pos.resize(mesh.n_vertices(),-1);
		for (int i = 0; i < vinfo.size(); i++)
		{
			return_cones_in_vinfo_pos[vinfo[i].idx_i] = i;
		}

		if (!posi_cone_status[i])
		{
			vector<int> candidate_delete;
			vector<DeleteInfo> deleinfo_vec;
			double before_dele_dist =calc_mesh_distortion(rho_vec);
			for (int j = 0; j < unposi_cone.size(); j++)
			{
				if (!unposi_cone_status[j])
				{
					targets[unposi_cone[j]] = 1;
					Algorithm::Dijkstra_with_nearest2(MC, posi_cone[i], targets, path);
					double distance = path.size() - 1;
					targets[unposi_cone[j]] = 0;
					path.clear();
					if (distance <= rings_num)
					{
						candidate_delete.push_back(unposi_cone[j]);
					}
				}
			}

			if (candidate_delete.size() > 0) 
			{
				vector<int> candidate_delete_pos;
				candidate_delete_pos.clear();
				candidate_delete_pos.resize(mesh.n_vertices(), -1);
				vec_vec_vinfo.clear();
				//vec_vec_vinfo.resize(candidate_delete.size());
				for (int k = 0; k < candidate_delete.size(); k++)
				{

					rho_vec[posi_cone[i]] = 0;
					rho_vec[candidate_delete[k]] = 0;
					VectorXd rho_vec_merge = rho_vec;
					vector<VInfo>judge_merge_vinfo;
					
					judge_merge_vinfo.clear();
					for (int q1 = 0; q1 < vinfo.size(); q1++)
					{
						if (q1 != return_cones_in_vinfo_pos[posi_cone[i]] && q1 != return_cones_in_vinfo_pos[candidate_delete[k]])
						{
							judge_merge_vinfo.push_back(vinfo[q1]);
						}
					}
					
					double dist = local_search2(rho_vec_merge, judge_merge_vinfo);
					double after_dele_dist = dist;
					double delta_dist = before_dele_dist - after_dele_dist;

					DeleteInfo dinfo;
					dinfo.posi_idx_i = posi_cone[i];
					dinfo.unposi_idx_i = candidate_delete[k];
					dinfo.delta_dist = delta_dist;
					dinfo.before_dele_dist = before_dele_dist;
					dinfo.after_dele_dist = after_dele_dist;
					deleinfo_vec.push_back(dinfo);

					candidate_delete_pos[candidate_delete[k]] = k;
					rho_vec[posi_cone[i]] = 1;
					rho_vec[candidate_delete[k]] = -1;
					
					now_mergy_time = clock();
					double now_time = (double)(now_mergy_time - start_mergy_time) / CLOCKS_PER_SEC;
					if (now_time>300)
					{
						goto flag_merge;
					}
				}

				std::sort(deleinfo_vec.begin(), deleinfo_vec.end(), com_smaller_delta_dist);

				if (deleinfo_vec[0].delta_dist >= 0)
				{
					rho_vec[posi_cone[i]] = 0;
					posi_cone_status[i] = true;
					rho_vec[deleinfo_vec[0].unposi_idx_i] = 0;

					auto it = find(unposi_cone.begin(), unposi_cone.end(), deleinfo_vec[0].unposi_idx_i);
					int order_idx = std::distance(unposi_cone.begin(), it);
					unposi_cone_status[order_idx] = true;

					//std::cout << "merge: " << posi_cone[i] << " and " << deleinfo_vec[0].unposi_idx_i << endl;

					vinfo = vec_vec_vinfo[candidate_delete_pos[deleinfo_vec[0].unposi_idx_i]];
					
					rho_vec.setZero(mesh.n_vertices());
					for (int q2 = 0; q2 < vinfo.size(); q2++)
					{
						rho_vec[vinfo[q2].idx_i] = vinfo[q2].after_round;
					}
					//cout << "after merging distortion: " << calc_mesh_distortion(rho_vec) << endl;
					double first_element_rho = vinfo[0].after_round;
					VInfo position2th;
					for (int q3 = 1; q3 < vinfo.size(); q3++)
					{
						double element_rho = vinfo[q3].after_round;
						if (first_element_rho * element_rho < 0)
						{
							position2th.idx_i = vinfo[q3].idx_i;
							position2th.is_cone = vinfo[q3].is_cone;
							position2th.after_reweighed_rho = vinfo[q3].after_reweighed_rho;
							position2th.after_round = vinfo[q3].after_round;
							position2th.gap = vinfo[q3].gap;
							vinfo.erase(vinfo.begin() + q3);
							vinfo.insert(vinfo.begin() + 1, position2th);
							break;
						}
					}

					//undata posi_vec
					posi_cone.clear();
					unposi_cone.clear();
					for (int i = 0; i < vinfo.size(); i++)
					{
						if (vinfo[i].after_round > 0)
						{
							posi_cone.push_back(vinfo[i].idx_i);
						}
						else
						{
							unposi_cone.push_back(vinfo[i].idx_i);
						}
					}
					posi_cone_status.resize(posi_cone.size(), false);
					unposi_cone_status.resize(unposi_cone.size(), false);

					non_cone.clear();
					for (int i = 0; i < rho_vec.size(); i++)
					{
						if (abs(rho_vec[i]) < 1e-4)
						{
							non_cone.push_back(i);
						}
					}
					goto flag0;
				}
				else
				{
					if (deleinfo_vec[0].after_dele_dist <= merge_bound * sigma_bound * 1.24)
					{
						if (deleinfo_vec[0].delta_dist >= -delta_sigma * sigma_bound * 1.24)
						{
							//cout << "after merging distortion is not greater than bound: " << posi_cone[i] << " " << deleinfo_vec[0].unposi_idx_i << endl;
							rho_vec[posi_cone[i]] = 0;
							posi_cone_status[i] = true;
							rho_vec[deleinfo_vec[0].unposi_idx_i] = 0;
							auto it = find(unposi_cone.begin(), unposi_cone.end(), deleinfo_vec[0].unposi_idx_i);
							int order_idx = std::distance(unposi_cone.begin(), it);
							unposi_cone_status[order_idx] = true;
							//std::cout << "merge: " << posi_cone[i] << " and " << deleinfo_vec[0].unposi_idx_i << endl;

							vinfo.clear();
							vinfo = vec_vec_vinfo[candidate_delete_pos[deleinfo_vec[0].unposi_idx_i]];
							rho_vec.setZero(mesh.n_vertices());
							
							for (int q2 = 0; q2 < vinfo.size(); q2++)
							{
								rho_vec[vinfo[q2].idx_i] = vinfo[q2].after_round;
							}
							//cout << "after merging distortion: " << calc_mesh_distortion(rho_vec) << endl;
							double first_element_rho = vinfo[0].after_round;
							VInfo position2th;
							for (int q3 = 1; q3 < vinfo.size(); q3++)
							{
								double element_rho = vinfo[q3].after_round;
								if (first_element_rho * element_rho < 0)
								{
									position2th.idx_i = vinfo[q3].idx_i;
									position2th.is_cone = vinfo[q3].is_cone;
									position2th.after_reweighed_rho = vinfo[q3].after_reweighed_rho;
									position2th.after_round = vinfo[q3].after_round;
									position2th.gap = vinfo[q3].gap;
									vinfo.erase(vinfo.begin() + q3);
									vinfo.insert(vinfo.begin() + 1, position2th);
									break;
								}
							}

							posi_cone.clear();
							unposi_cone.clear();
							for (int i = 0; i < vinfo.size(); i++)
							{
								if (vinfo[i].after_round > 0)
								{
									posi_cone.push_back(vinfo[i].idx_i);
								}
								else
								{
									unposi_cone.push_back(vinfo[i].idx_i);
								}
							}
							posi_cone_status.resize(posi_cone.size(), false);
							unposi_cone_status.resize(unposi_cone.size(), false);

							non_cone.clear();
							for (int i = 0; i < rho_vec.size(); i++)
							{
								if (abs(rho_vec[i]) < 1e-4)
								{
									non_cone.push_back(i);
								}
							}
							goto flag0;
						}
					}
				}
				candidate_delete.clear();
				deleinfo_vec.clear();
			}
		}
	}

	flag_merge:
	non_cone.clear();
	for (int i = 0; i < rho_vec.size(); i++)
	{
		if (abs(rho_vec[i]) < 1e-4)
		{
			non_cone.push_back(i);
		}
	}
}

