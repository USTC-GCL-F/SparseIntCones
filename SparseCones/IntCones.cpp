#include "IntCones.h"
#pragma comment (lib,"ReweighedSOCP.lib")

IntCones::IntCones(Mesh& mesh) :mesh(mesh)
{

}

void IntCones::VertArea()
{
	vec_A.setConstant(mesh.n_vertices(), 0);
	double fA;
	for (const FH& f_h : mesh.faces())
	{
		fA = mesh.calc_face_area(f_h) / 3;
		for (const VH& fv_h : mesh.fv_range(f_h))
		{
			vec_A[fv_h.idx()] += fA;
		}
	}

	vec_A = vec_A / vec_A.sum();

	/*
	std::ofstream _out("A_dog.txt");
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		_out << vec_A[i] << std::endl;
	}
	_out.close();*/
}

void IntCones::FaceAngle()
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

void IntCones::VertGauss()
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

	/*
	std::ofstream _out("K_ori_dog.txt");
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		_out << K_ori[i] << std::endl;
	}
	_out.close();
	std::cout << " " << endl;*/
}


void IntCones::MatrixL()
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

	//cout << "������˹���������ʽ: "<<(m_L.toDense()).determinant() << endl;
	/*
	std::ofstream _out("L_dog.txt");
	for (int k = 0; k < m_L.outerSize(); ++k)
	{
		for (SparseMatrix<double>::InnerIterator it(m_L, k); it; ++it)
		{


			_out << it.row()+1<<" "<< it.col()+1<<" "<< it.value() << std::endl;
			//std::cout << it.value() << std::endl;
			//it.row();   // row index
			//it.col();   // col index (here it is equal to k)
			//it.index(); // inner index, here it is equal to it.row()
		}
	}*/
	//m_C = m_L * m_M.cwiseInverse().asDiagonal() * m_L.transpose();
}

void IntCones::Init()
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

// 1.load txt
vector<int> IntCones::load_rho_int(string str)
{
	vector<int> rho_round;
	rho_round.clear();
	ifstream ifs;
	ifs.open(str, ios::in);
	vector<string> item;
	string temp;
	while (getline(ifs, temp))
	{
		item.push_back(temp);
	}
	for (auto it = item.begin(); it != item.end(); it++)
	{
		istringstream istr(*it);
		string str;
		int count = 0;
		while (istr >> str)
		{
			if (count == 0)
			{
				int r = atoi(str.c_str());
				rho_round.push_back(r);
			}
			count++;
		}
	}
	return rho_round;

}

vector<double> IntCones::load_rho_double(string str)
{
	vector<double> rho_socp;
	rho_socp.clear();
	ifstream ifs2;
	ifs2.open(str, ios::in);
	vector<string> item2;
	string temp2;
	while (getline(ifs2, temp2))
	{
		item2.push_back(temp2);
	}
	for (auto it = item2.begin(); it != item2.end(); it++)
	{
		istringstream istr(*it);
		string str;
		int count = 0;
		while (istr >> str)
		{
			if (count == 0)
			{
				double r = atof(str.c_str());
				rho_socp.push_back(r);
			}
			count++;
		}
	}
	return rho_socp;
}

void IntCones::init_vinfo_noncone(vector<double>& a, vector<double>& b)
{
	//vector<VInfo> vinfo;
	//vector<int> non_cone;
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

	// (1).��֤ǰ����Ԫ����һ��һ��
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
	//vinfo.insert(vinfo.begin() + 1, position2th);
}

void IntCones::init_equatuion_data()
{
	Init();
	E.setOnes(mesh.n_vertices());
	rho_vec.resize(mesh.n_vertices());

	//std::cout << "gauss_bonet: " << rho_vec.sum() << endl;
	solver.compute(PT * m_L * P);
	vector<pair<int, double>> compare;

	// 3.
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
	//Eigen::VectorXd U1 = U - ((U.transpose() * vec_A.asDiagonal() * E) / (E.transpose() * vec_A.asDiagonal() * E)) * E;
	distortion = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
	std::cout << "��ʼ�����ʼ (�ڶ��� round ��) ��Ť���� " << distortion << endl;
	 
	search_merge_info.push_back(distortion);

	cout << endl;
}

double IntCones::calc_mesh_distortion(Eigen::VectorXd rho_vec)
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
	//Eigen::VectorXd U1 = U - ((U.transpose() * vec_A.asDiagonal() * E) / (E.transpose() * vec_A.asDiagonal() * E)) * E;
	double d = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);

	final_U = U1;


	return d;
}

void IntCones::global_search()
{
	for (int i = 0; i < non_cone.size(); i++)
	{
		rho_vec[non_cone[i]] = rj;
		//std::cout << " - - - - - - - - - - " << endl;
		//std::cout << "gauss_bonet: " << rho_vec.sum() <<" "<<rj << endl;

		distortion = calc_mesh_distortion(rho_vec);
		compare.push_back(make_pair(non_cone[i], distortion));

		//std::cout << "����ֵ: " << compare[i].first << " " << "Ť��ֵ: " << compare[i].second << endl;
		rho_vec[non_cone[i]] = 0;
		//std::cout << " - - - - - - - - - - " << endl;
		//cout << endl;
	}
}

void IntCones::local_search(vector<int>& a)
{


	for (int i = 0; i < a.size(); i++)
	{
		rho_vec[a[i]] = rj;
		//std::cout << " - - - - - - - - - - " << endl;
		//std::cout << "gauss_bonet: " << rho_vec.sum() <<" "<<rj << endl;

		double distortion = calc_mesh_distortion(rho_vec);//*

		compare.push_back(make_pair(a[i], distortion));
		//std::cout << "����ֵ: " << compare[i].first << " " << "Ť��ֵ: " << compare[i].second << endl;
		rho_vec[a[i]] = 0;
		//std::cout << " - - - - - - - - - - " << endl;
		//cout << endl;
	}


}

void IntCones::merge_and_updata_vinfo(int rings_num)
{
	// 3.merge ��ŵ�
		// 1.�� 1 �� 1
	vector<int> posi_cone;
	vector<bool> posi_cone_status;
	vector<int> unposi_cone;
	vector<bool> unposi_cone_status;
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

	MeshCache MC(mesh);
	std::vector<int> path;
	std::vector<int> targets(MC.NVertices, 0);
	//targets[100] = 1;

	// 2.distance
	cout << endl;
	for (int i = 0; i < posi_cone.size(); i++)
	{
		if (!posi_cone_status[i])
		{
			vector<int> candidate_delete;
			vector<DeleteInfo> deleinfo_vec;
			//double before_dele_dist = calc_mesh_distortion(rho_vec);

			// ***********************************
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
			//Eigen::VectorXd U1 = U - ((U.transpose() * vec_A.asDiagonal() * E) / (E.transpose() * vec_A.asDiagonal() * E)) * E;
			double before_dele_dist = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
			//*************************************

			for (int j = 0; j < unposi_cone.size(); j++)
			{
				if (!unposi_cone_status[j])
				{
					// distance > n_rings
					// 1. calc distance
					targets[unposi_cone[j]] = 1;
					Algorithm::Dijkstra_with_nearest2(MC, posi_cone[i], targets, path);
					double distance = path.size() - 1;
					targets[unposi_cone[j]] = 0;
					path.clear();
					// 2. 
					if (distance <= rings_num)
					{
						candidate_delete.push_back(unposi_cone[j]);
					}
				}
			}

			if (candidate_delete.size() > 0) // ������������cone
			{
				for (int k = 0; k < candidate_delete.size(); k++)
				{
					rho_vec[posi_cone[i]] = 0;
					rho_vec[candidate_delete[k]] = 0;

					//*****************************************
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
					//Eigen::VectorXd U1 = U - ((U.transpose() * vec_A.asDiagonal() * E) / (E.transpose() * vec_A.asDiagonal() * E)) * E;
					double after_dele_dist = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
					//********************************************

					double delta_dist = before_dele_dist - after_dele_dist;

					DeleteInfo dinfo;
					dinfo.posi_idx_i = posi_cone[i];
					dinfo.unposi_idx_i = candidate_delete[k];
					dinfo.delta_dist = delta_dist;
					dinfo.before_dele_dist = before_dele_dist;
					dinfo.after_dele_dist = after_dele_dist;
					deleinfo_vec.push_back(dinfo);

					rho_vec[posi_cone[i]] = 1;
					rho_vec[candidate_delete[k]] = -1;

				}


				std::sort(deleinfo_vec.begin(), deleinfo_vec.end(), com_smaller_delta_dist);

				// zhzh 1029
				if (deleinfo_vec[0].delta_dist >=0)
				{
					rho_vec[posi_cone[i]] = 0;
					posi_cone_status[i] = true;
					rho_vec[deleinfo_vec[0].unposi_idx_i] = 0;
					auto it = find(unposi_cone.begin(), unposi_cone.end(), deleinfo_vec[0].unposi_idx_i);
					int order_idx = std::distance(unposi_cone.begin(), it);
					unposi_cone_status[order_idx] = true;
					//std::cout << endl;
					std::cout << "�ϲ�: " << posi_cone[i] << " �� " << deleinfo_vec[0].unposi_idx_i << endl;
				}
				else
				{
					if (deleinfo_vec[0].after_dele_dist<= merge_bound* sigma_bound* ratio_sigma)
					{
						if (deleinfo_vec[0].delta_dist >= -delta_sigma * sigma_bound)
						{
							rho_vec[posi_cone[i]] = 0;
							posi_cone_status[i] = true;
							rho_vec[deleinfo_vec[0].unposi_idx_i] = 0;
							auto it = find(unposi_cone.begin(), unposi_cone.end(), deleinfo_vec[0].unposi_idx_i);
							int order_idx = std::distance(unposi_cone.begin(), it);
							unposi_cone_status[order_idx] = true;
							//std::cout << endl;
							std::cout << "�ϲ�: " << posi_cone[i] << " �� " << deleinfo_vec[0].unposi_idx_i << endl;
						}
						
					}
				}
				candidate_delete.clear();
				deleinfo_vec.clear();
			}

		}
	}

	// zhzh 1106
	merge_bound *= 0.9;
	/*
	std::ofstream _out2("vinfo_before" + to_string(y) + ".csv");
	for (int i = 0; i < vinfo.size(); i++)
	{
		_out2 << vinfo[i].idx_i<< std::endl;
	}
	_out2.close();*/

	// ���� vinfo ������ true ����ɾ��
	
	// 1029 ************
	/*
	vector<VInfo> replace_vinfo;
	replace_vinfo.resize(vinfo.size());
	std::copy(vinfo.begin() + global_search_num, vinfo.end(), replace_vinfo.begin());
	std::copy(vinfo.begin(), vinfo.begin() + global_search_num, replace_vinfo.begin()+global_search_num);
	vinfo.clear();
	vinfo = replace_vinfo;
	replace_vinfo.clear();
	global_search_num = 0;*/
	// 1029 ************

	vector<int> vinfo_order;
	vinfo_order.clear();
	for (int i = 0; i < vinfo.size(); i++)
	{
		vinfo_order.push_back(vinfo[i].idx_i);
	}

	vector<int> posi_position;
	for (int i = 0; i < posi_cone_status.size(); i++)
	{
		if (!posi_cone_status[i])
		{
			posi_position.push_back(i);
		}
	}

	vector<int> unposi_position;
	for (int i = 0; i < unposi_cone_status.size(); i++)
	{
		if (!unposi_cone_status[i])
		{
			unposi_position.push_back(i);
		}
	}

	vector<int> dele_vinfo;
	for (int i = 0; i < posi_position.size(); i++)
	{
		auto iter = find(vinfo_order.begin(), vinfo_order.end(), posi_cone[posi_position[i]]);
		int distannn = distance(vinfo_order.begin(), iter);
		dele_vinfo.push_back(distannn);
		//vinfo.erase(vinfo.begin()+distannn);
	}

	for (int i = 0; i < unposi_position.size(); i++)
	{
		auto iter = find(vinfo_order.begin(), vinfo_order.end(), unposi_cone[unposi_position[i]]);
		int distannn = distance(vinfo_order.begin(), iter);
		dele_vinfo.push_back(distannn);
		//vinfo.erase(vinfo.begin() + distannn);
	}

	vector<VInfo> candi_vinfo;
	for (int i = 0; i < dele_vinfo.size(); i++)
	{
		candi_vinfo.push_back(vinfo[dele_vinfo[i]]);
	}
	vinfo.clear();
	vinfo = candi_vinfo;


	// zhzh 1029
	std::sort(vinfo.begin(), vinfo.end(), com_smaller_gap);
	// zhzh 1029

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

	non_cone.clear();
	for (int i = 0; i < rho_vec.size(); i++)
	{
		if (abs(rho_vec[i]) < 1e-4)
		{
			non_cone.push_back(i);
		}
	}
	/*
	std::ofstream _out3("vinfo_after" + to_string(y) + ".csv");
	for (int i = 0; i < vinfo.size(); i++)
	{
		_out3 << vinfo[i].idx_i << std::endl;
	}
	_out3.close();

	std::ofstream _out4("posi_" + to_string(y) + ".csv");
	for (int i = 0; i < posi_cone.size(); i++)
	{
		_out4<< posi_cone[i]<<","<< posi_cone_status[i] << std::endl;
	}
	_out4.close();

	std::ofstream _out5("unposi_" + to_string(y) + ".csv");
	for (int i = 0; i < unposi_cone.size(); i++)
	{
		_out5 << unposi_cone[i] << "," << unposi_cone_status[i] << std::endl;
	}
	_out5.close();*/



	posi_cone.clear();
	posi_cone_status.clear();
	unposi_cone.clear();
	unposi_cone_status.clear();
}

void IntCones::calc_cones_position(vector<double> socp_res, vector<double> round_res)
{

	clock_t startTime, endTime;
	startTime = clock();
	// 1.load txt
	
	/*
	std::string str1 = "C://Users//USTC-GCL//Desktop//armadillo//018//round_10.txt";
	std::string str2 = "C://Users//USTC-GCL//Desktop//armadillo//018//res_0.18_10.txt";
	vector<double> vec_rho = load_rho_double(str1);
	vector<double> after_reweighed_rho = load_rho_double(str2);*/

	vector<double> vec_rho = round_res;
	vector<double> after_reweighed_rho = socp_res;


	//cout << "end!!!" << endl;

	// 2.init vinfo and non_cone
	init_vinfo_noncone(vec_rho, after_reweighed_rho);

	// 3.init data and calc init distortion ( after around ) 
	rho_vec = Eigen::Map<VectorXd>(vec_rho.data(), mesh.n_vertices());
	init_equatuion_data();

	// 4.loop {global, local(break), merge}_n
	for (int y = 0; y < loop_mix_num; y++)
	{
		int global_iter = 1;
		cout << endl;
		cout << "�� " << y + 1 << " �ε�����ʼ *************************************" << endl;
		for (int j = 0; j < vinfo.size(); j++)
		{
			if (global_iter <= 2)
			{
				rho_vec[vinfo[j].idx_i] = 0;
				std::cout << "�ƶ�ǰ: " << vinfo[j].idx_i << "  **  ";
				non_cone.push_back(vinfo[j].idx_i);
				rj = vinfo[j].after_round;

				global_search();

				std::sort(compare.begin(), compare.end(), cmp_bigger);
				rho_vec[compare[0].first] = rj;

				vector<int>::iterator it = find(non_cone.begin(), non_cone.end(), compare[0].first);
				non_cone.erase(it);
				non_cone.push_back(vinfo[j].idx_i);

				vinfo[j].idx_i = compare[0].first;//����cone��λ��

				std::cout << " ȫ������(����һ��һ��)  ***  �ƶ���(��СŤ������ֵ): " << compare[0].first << " " << "��ӦŤ��ֵ: " << compare[0].second << "  roundֵ: " << " " << rj << endl;

				distortion = compare[0].second;

				compare.clear();
				non_cone.pop_back();

				global_iter++;
				global_search_num++;
			}
			else
			{
				if (distortion > ratio_sigma * sigma_bound) //ratio_sigma*0.2
				{
					rho_vec[vinfo[j].idx_i] = 0;
					std::cout << "�ƶ�ǰ: " << vinfo[j].idx_i << "  **  ";
					non_cone.push_back(vinfo[j].idx_i);
					rj = vinfo[j].after_round;

					global_search();

					std::sort(compare.begin(), compare.end(), cmp_bigger);
					rho_vec[compare[0].first] = rj;

					vector<int>::iterator it = find(non_cone.begin(), non_cone.end(), compare[0].first);
					non_cone.erase(it);
					non_cone.push_back(vinfo[j].idx_i);

					vinfo[j].idx_i = compare[0].first;//����cone��λ��

					std::cout << " ȫ������  ***  �ƶ���(��СŤ������ֵ): " << compare[0].first << " " << "��ӦŤ��ֵ: " << compare[0].second << "  roundֵ: " << " " << rj << endl;

					distortion = compare[0].second;

					compare.clear();
					non_cone.pop_back();
					global_search_num++;
				}
				else
				{
					set<int> n_rings_set;
					vector<int> n_rings_vec;
					n_rings_set.clear();
					n_rings_vec.clear();
					auto v = mesh.vertex_handle(vinfo[j].idx_i);
					n_rings_set.insert(vinfo[j].idx_i);
					// ȥ��cone��
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
					std::cout << "�ƶ�ǰ: " << vinfo[j].idx_i;
					n_rings_vec.push_back(vinfo[j].idx_i);
					rj = vinfo[j].after_round;

					local_search(n_rings_vec);

					std::sort(compare.begin(), compare.end(), cmp_bigger);
					rho_vec[compare[0].first] = rj;
					vinfo[j].idx_i = compare[0].first;//����cone��λ��

					std::cout << " ** ��������  ***  �ƶ���(��СŤ������ֵ): " << compare[0].first << " " << " ��ӦŤ��ֵ: " << compare[0].second << "  roundֵ: " << " " << rj << endl;

					compare.clear();
					n_rings_set.clear();
					n_rings_vec.clear();
				}
			}
		}

		if (y == loop_mix_num - 1)
		{
			double stop_dist = calc_mesh_distortion(rho_vec);
			cout << endl;
			cout << "��ѭ�������Ť��: " << stop_dist << endl;

			search_merge_info.push_back(stop_dist);
			std::ofstream _out5("final_cone_new.txt");
			int cones_number = 0;
			for (int i = 0; i < rho_vec.size(); i++)
			{
				_out5 << rho_vec[i] << std::endl;
				if (rho_vec[i] !=0)
				{
					cones_number++;
				}
			}
			_out5.close();
			search_merge_info.push_back(cones_number);
			break;
		}

		// merge and undata vinfo
		merge_and_updata_vinfo(search_rings_num);

		// calc y_th distortion
		double final_distortion = calc_mesh_distortion(rho_vec);

		cout << "�� " << y + 1 << " �ε�������, ���յ�Ť��: " << final_distortion << " **********************************************************" << endl;
		cout << endl;

		/*
		std::ofstream _out1("rho_mix_merge_points_new_" + to_string(y) + ".txt");
		for (int i = 0; i < rho_vec.size(); i++)
		{
			_out1 << rho_vec[i] << std::endl;
		}
		_out1.close();*/
	}


	endTime = clock();
	std::cout << endl;
	std::cout << "������ �����ϲ� ��ʱ��: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

	search_merge_info.push_back((double)(endTime - startTime) / CLOCKS_PER_SEC);


}

double IntCones::compute_dist(vector<double> a)
{
	VectorXd b = Eigen::Map<VectorXd>(a.data(), mesh.n_vertices());
	cout <<"a: "<< a.size() << endl;
	cout << "b: " << b.size() << endl;
	cout << "v: " << mesh.n_vertices() << endl;
	
	
	Init();
	E.setOnes(mesh.n_vertices());

	//std::cout << "gauss_bonet: " << rho_vec.sum() << endl;
	solver.compute(PT * m_L * P);
	cout << endl;
	double c =calc_mesh_distortion(b);
	cout << "zz"<< endl;
	cout << c << endl;
	return c;
}


void IntCones::ori_mesh_local_search(vector<VInfo> vinfo,int loop_num)
{
	clock_t start, end;
	start = clock();
	Init();
	
	E.setOnes(mesh.n_vertices());
	//cout << "1" << endl;
	rho_vec.setZero(mesh.n_vertices());
	/*
	for (int i = 0; i < rho_vec.size(); i++)
	{
		cout << rho_vec[i] << endl;
	}*/
	//cout << rho_vec.size() << endl;
	//cout << "2" << endl;
	//cout << vinfo.size() << endl;
	//cout << "3" << endl;
	for (int i=0;i<vinfo.size();i++)
	{
		rho_vec[vinfo[i].idx_i] = vinfo[i].after_round;
		//cout << vinfo[i].idx_i << " " << vinfo[i].after_round << endl;
	}
	//cout << "4" << endl;
	//std::cout << "gauss_bonet: " << rho_vec.sum() << endl;
	solver.compute(PT * m_L * P);
	//vector<pair<int, double>> compare;

	// 3.
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
	//Eigen::VectorXd U1 = U - ((U.transpose() * vec_A.asDiagonal() * E) / (E.transpose() * vec_A.asDiagonal() * E)) * E;
	distortion = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
	std::cout << "ͶӰ��ĳ�ʼŤ���� " << distortion << endl;
	//search_merge_info.push_back(distortion);
	ori_mesh_local_info.clear();
	ori_mesh_local_info.push_back(distortion);
	ori_mesh_local_info.push_back(vinfo.size());
	cout << endl;
	double final_local_search_distortion = 0;

	//cout << "2" << endl;

	for (int m = 0; m < loop_num; m++)
	{
		cout << endl;
		//cout << "�� " << m + 1 << " ��ԭʼ����ľֲ�����" << endl;
		for (int j = 0; j < vinfo.size(); j++)
		{
			set<int> n_rings_set;
			vector<int> n_rings_vec;
			n_rings_set.clear();
			n_rings_vec.clear();
			auto v = mesh.vertex_handle(vinfo[j].idx_i);
			n_rings_set.insert(vinfo[j].idx_i);
			// ȥ��cone��
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
			std::cout << "�ƶ�ǰ: " << vinfo[j].idx_i;
			n_rings_vec.push_back(vinfo[j].idx_i);
			rj = vinfo[j].after_round;

			local_search(n_rings_vec);

			std::sort(compare.begin(), compare.end(), cmp_bigger);
			rho_vec[compare[0].first] = rj;
			vinfo[j].idx_i = compare[0].first;//����cone��λ��

			std::cout << " ** ͶӰ��ԭʼ���� ** ��������  ***  �ƶ���(��СŤ������ֵ): " << compare[0].first << " " << " ��ӦŤ��ֵ: " << compare[0].second << "  roundֵ: " << " " << rj << endl;
			final_local_search_distortion = compare[0].second;
			compare.clear();
			n_rings_set.clear();
			n_rings_vec.clear();
		}
	}
	end = clock();
	ori_mesh_local_info.push_back((double)(end - start) / CLOCKS_PER_SEC );
	ori_mesh_local_info.push_back(final_local_search_distortion);

	std::ofstream _out2("final_rho.csv");
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		_out2 << rho_vec[i] << std::endl;
	}
	_out2.close();




}


#if 0
void IntCones::ori_mesh_local_search_kai(vector<KaiCones::VInfo> vinfo, int loop_num)
{
	clock_t start, end;
	start = clock();
	Init();

	E.setOnes(mesh.n_vertices());
	//cout << "1" << endl;
	rho_vec.setZero(mesh.n_vertices());
	/*
	for (int i = 0; i < rho_vec.size(); i++)
	{
		cout << rho_vec[i] << endl;
	}*/
	//cout << rho_vec.size() << endl;
	//cout << "2" << endl;
	//cout << vinfo.size() << endl;
	//cout << "3" << endl;
	for (int i = 0; i < vinfo.size(); i++)
	{
		rho_vec[vinfo[i].idx_i] = vinfo[i].after_round;
		//cout << vinfo[i].idx_i << " " << vinfo[i].after_round << endl;
	}
	//cout << "4" << endl;
	//std::cout << "gauss_bonet: " << rho_vec.sum() << endl;
	solver.compute(PT * m_L * P);
	//vector<pair<int, double>> compare;

	// 3.
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
	//Eigen::VectorXd U1 = U - ((U.transpose() * vec_A.asDiagonal() * E) / (E.transpose() * vec_A.asDiagonal() * E)) * E;
	distortion = sqrt(U1.transpose() * vec_A.asDiagonal() * U1);
	std::cout << "ͶӰ��ĳ�ʼŤ���� " << distortion << endl;
	//search_merge_info.push_back(distortion);
	ori_mesh_local_info.clear();
	ori_mesh_local_info.push_back(distortion);
	ori_mesh_local_info.push_back(vinfo.size());
	cout << endl;
	double final_local_search_distortion = 0;

	//cout << "2" << endl;

	for (int m = 0; m < loop_num; m++)
	{
		cout << endl;
		cout << "�� " << m + 1 << " ��ԭʼ����ľֲ�����" << endl;
		for (int j = 0; j < vinfo.size(); j++)
		{
			set<int> n_rings_set;
			vector<int> n_rings_vec;
			n_rings_set.clear();
			n_rings_vec.clear();
			auto v = mesh.vertex_handle(vinfo[j].idx_i);
			n_rings_set.insert(vinfo[j].idx_i);
			// ȥ��cone��
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
			std::cout << "�ƶ�ǰ: " << vinfo[j].idx_i;
			n_rings_vec.push_back(vinfo[j].idx_i);
			rj = vinfo[j].after_round;

			local_search(n_rings_vec);

			std::sort(compare.begin(), compare.end(), cmp_bigger);
			rho_vec[compare[0].first] = rj;
			vinfo[j].idx_i = compare[0].first;//����cone��λ��

			std::cout << " ** ͶӰ��ԭʼ���� ** ��������  ***  �ƶ���(��СŤ������ֵ): " << compare[0].first << " " << " ��ӦŤ��ֵ: " << compare[0].second << "  roundֵ: " << " " << rj << endl;
			final_local_search_distortion = compare[0].second;
			compare.clear();
			n_rings_set.clear();
			n_rings_vec.clear();
		}
	}
	end = clock();
	ori_mesh_local_info.push_back((double)(end - start) / CLOCKS_PER_SEC);
	ori_mesh_local_info.push_back(final_local_search_distortion);

	std::ofstream _out2("final_rho.csv");
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		_out2 << rho_vec[i] << std::endl;
	}
	_out2.close();




}
#endif
