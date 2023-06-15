#include "Prune.h"

Prune::Prune(Mesh& mesh):mesh(mesh)
{

}

void Prune::vert_area()
{
	A_vec.setConstant(mesh.n_vertices(), 0);
	double fA;
	for (const FH& f_h : mesh.faces())
	{
		fA = mesh.calc_face_area(f_h) / 3;
		for (const VH& fv_h : mesh.fv_range(f_h))
		{
			A_vec[fv_h.idx()] += fA;
		}
	}

	A_vec = A_vec / A_vec.sum();

}
void Prune::face_angle()
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
void Prune::vert_gauss()
{
	K_vec.setConstant(mesh.n_vertices(), 2 * M_PI);

	for (const VH& v_h : mesh.vertices())
	{
		if (mesh.is_boundary(v_h)) {
			K_vec[v_h.idx()] = M_PI;
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

		K_vec[v_ids[0]] -= f_angle[3 * f_h.idx() + 0];
		K_vec[v_ids[1]] -= f_angle[3 * f_h.idx() + 1];
		K_vec[v_ids[2]] -= f_angle[3 * f_h.idx() + 2];
	}

	K_0.resize(mesh.n_vertices());
	for (int j = 0; j < mesh.n_vertices(); j++)
	{
		K_0[j] = K_vec[j];
	}

}

void Prune::calc_laplace()
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

	L_matrix.resize(mesh.n_vertices(), mesh.n_vertices());
	L_matrix.setFromTriplets(trip.begin(), trip.end());

}

void Prune::generate_mesh_info()
{
	vert_area();
	face_angle();
	vert_gauss();
	calc_laplace();

	E.setOnes(mesh.n_vertices());

	P.resize(mesh.n_vertices(), mesh.n_vertices() - 1);
	std::vector<T> element;
	for (int i = 0; i < mesh.n_vertices() - 1; i++)
	{
		element.emplace_back(i, i, 1);
	}

	P.setFromTriplets(element.begin(), element.end());
	PT = P.transpose();

	SpMat L_multi_A_sqrt_inv = L_matrix * A_vec.cwiseSqrt().asDiagonal().inverse();
	int q = 0;
	int nv = mesh.n_vertices();
	mat_rows.clear();
	mat_cols.clear();
	mat_vals.clear();
	for (int k = 0; k < L_matrix.outerSize(); k++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(L_matrix, k); it; ++it)
		{
			double vals = it.value();
			int row = it.row();
			int cols = it.col();
			mat_rows.push_back(row);
			mat_cols.push_back(cols + 1);
			mat_vals.push_back(L_multi_A_sqrt_inv.coeffRef(row, cols));
			q++;
		}
	}
	L_element_num = q;
	mat_rows.resize(L_element_num + nv);
	mat_cols.resize(L_element_num + nv);
	mat_vals.resize(L_element_num + nv);
}

void Prune::socp(int vertex_num,double sigma,int loop_num_socp)
{
	
	generate_mesh_info();
	int n = vertex_num;
	W_vec.resize(n, 1);
	double threshold = 1;
	VectorXd rho_last;
	rho_last.setZero(n);
	int k = 0;
	while (k < loop_num_socp && threshold>1e-6)
	{
		Model::t M = new Model("cqo1"); auto _M = finally([&]() { M->dispose(); });
		Variable::t x = M->variable("x", 2 * n + 2, Domain::unbounded());
		M->constraint(x->index(0), Domain::inRange(0, sigma));
		M->constraint(x->index(n + 1), Domain::greaterThan(0));

		for (int i = 1; i < n + 1; i++)
		{
			M->constraint(x->index(n + 1 + i), Domain::inRange(-sqrt(W_vec[i - 1]), sqrt(W_vec[i - 1])));
		}
		for (int q = 0; q < n; q++)
		{
			mat_rows[L_element_num + q] = q;
			mat_cols[L_element_num + q] = q + n + 2;
			mat_vals[L_element_num + q] = -1 * 0.5 * M_PI * pow(sqrt(W_vec[q]), -1);
		}
		auto rows = monty::new_array_ptr<int>(mat_rows);
		auto cols = monty::new_array_ptr<int>(mat_cols);
		auto values = monty::new_array_ptr<double>(mat_vals);
		auto A = mosek::fusion::Matrix::sparse(n, 2 * n + 2, rows, cols, values);
		auto B = monty::new_array_ptr<double>(K_0);

		M->constraint(Expr::add(Expr::mul(A, x), B), Domain::equalsTo(0.0));
		Variable::t z1 = Var::vstack(x->index(0), x->slice(1, n + 1));
		Variable::t z2 = Var::vstack(x->index(n + 1), x->slice(n + 2, 2 * n + 2));
		Constraint::t qc1 = M->constraint("qc1", z1, Domain::inQCone());
		Constraint::t qc2 = M->constraint("qc2", z2, Domain::inQCone());
		
		// start opt
		M->objective("obj", ObjectiveSense::Minimize, x->index(n + 1));
		M->solve();

		// out result
		ndarray<double, 1> S2 = *((x->slice(n + 2, 2 * n + 2))->level());
		ndarray<double, 1> S1 = *((x->slice(1, n + 1))->level());
		ndarray<double, 1> S3 = *((x->slice(0, 1))->level());
		
		// update w and rho
		rho.clear();
		for (int i = 0; i < S2.size(); i++)
		{
			rho.push_back(S2[i] / sqrt(W_vec[i]));
		}

		for (int i = 0; i < W_vec.size(); i++)
		{
			double W_val = std::pow(std::pow(rho[i], 2) + 1e-10, -2);
			W_vec[i] = W_val / 1e10;
		}

		int nnn = 0;
		for (int i = 0; i < rho.size(); i++)
		{
			if (abs(rho[i]) > 1e-8)
			{
				nnn++;
			}

		}

		std::cout << "Iters:  " << k + 1 << "  rho norm2£º" << vec_norm_2(rho) << "  cones numbers£º " << nnn << "  distortion£º " << sqrt(mosek_norm_2(S1)) << "   Bound: " << S3[0] << endl;//" "<< accumulate(rho.begin(),rho.end(),0.0) <<
		
		VectorXd rho_ = Eigen::Map<VectorXd>(rho.data(), n);
		threshold = (rho_ - rho_last).norm() / rho_.norm();
		rho_last = rho_;
		test_socp_num++;
		socp_rho = rho;
		k++;
	}
	std::cout << "************ IRLS End ************" << endl;
}

void Prune::interval_round(int n, int loop_num,double delta, double left, double left_sigma, double mid, double mid_sigma, double right, double right_sigma)
{
	clock_t   start, end;
	start = clock();
	int vertex_num = n;
	std::cout << endl;
	std::cout << "************ Interval Round Start ************" << endl;
	
	int eplison = 1;
	double delta_dis = 10;
	double distortion_0 = 0;

	int q = 0;
	int nv = n;
	mat_rows.clear();
	mat_cols.clear();
	mat_vals.clear();
	SpMat L_multi_A_sqrt_inv = L_matrix * A_vec.cwiseSqrt().asDiagonal().inverse();
	for (int k = 0; k < L_matrix.outerSize(); k++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(L_matrix, k); it; ++it)
		{
			double vals = it.value();
			int row = it.row();
			int cols = it.col();
			mat_rows.push_back(row);
			mat_cols.push_back(cols + 1);
			mat_vals.push_back(L_multi_A_sqrt_inv.coeffRef(row, cols));
			q++;
		}
	}
	L_element_num = q;
	mat_rows.resize(L_element_num + nv);
	mat_cols.resize(L_element_num + nv);
	mat_vals.resize(L_element_num + nv);

	for (int q = 0; q < n; q++)
	{
		mat_rows[L_element_num + q] = q;
		mat_cols[L_element_num + q] = q + n + 1;// notice idx
		mat_vals[L_element_num + q] = -1 * 0.5 * M_PI;
	}

	while ((eplison <= loop_num) && (abs(delta_dis) > delta))
	{
		Model::t M = new Model("cqo1"); auto _M = finally([&]() { M->dispose(); });
		Variable::t x = M->variable("x", 2 * n + 1, Domain::unbounded());
		M->constraint(x->index(0), Domain::greaterThan(0));

		for (int i = 0; i < mesh.n_vertices(); i++)
		{
			if (abs(rho[i]) < mid + mid_sigma * eplison)
			{
				M->constraint(x->index(n + 1 + i), Domain::equalsTo(0));
			}
			else if (rho[i] >= right + right_sigma * eplison)
			{
				M->constraint(x->index(n + 1 + i), Domain::equalsTo(1));
			}
			else if (rho[i] <= left + left_sigma * eplison)
			{
				M->constraint(x->index(n + 1 + i), Domain::equalsTo(-1));
			}
			else
			{
				M->constraint(x->index(n + 1 + i), Domain::inRange(-1, 1));
			}
		}

		auto rows = monty::new_array_ptr<int>(mat_rows);
		auto cols = monty::new_array_ptr<int>(mat_cols);
		auto values = monty::new_array_ptr<double>(mat_vals);
		auto A = mosek::fusion::Matrix::sparse(n, 2 * n + 1, rows, cols, values);
		auto B = monty::new_array_ptr<double>(K_0);

		M->constraint(Expr::add(Expr::mul(A, x), B), Domain::equalsTo(0.0));

		Variable::t z1 = Var::vstack(x->index(0), x->slice(1, n + 1));
		Constraint::t qc1 = M->constraint("qc1", z1, Domain::inQCone());

		M->objective("obj", ObjectiveSense::Minimize, x->index(0));
		M->solve();

		ndarray<double, 1> S1 = *((x->slice(1, n + 1))->level());
		ndarray<double, 1> rho_array = *((x->slice(n + 1, 2 * n + 1))->level());

		for (int i = 0; i < n; i++)
		{
			rho[i] = rho_array(i);
		}

		double dis_1 = sqrt(mosek_norm_2(S1));
		delta_dis = dis_1 - distortion_0;
		distortion_0 = dis_1;
		eplison += 1;
		std::cout << "Iters: " << eplison-1 << " distortion:  " << sqrt(mosek_norm_2(S1)) << endl;
	}

	interval_rho.resize(mesh.n_vertices());
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		interval_rho[i] = rho[i];
	}

	int posi_iteger_cone_numbers = 0;
	int unposi_iteger_cone_numbers = 0;
	int non_cone_numbers = 0;

	non_iteger_cone_numbers = 0;
	iteger_cone_numbers = 0;
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		if (rho[i] == 1)
		{
			posi_iteger_cone_numbers++;
		}
		else if (rho[i] == -1)
		{
			unposi_iteger_cone_numbers++;
		}
		else if (rho[i] == 0)
		{
			non_cone_numbers++;
		}
		else
		{
			non_iteger_cone_numbers++;
		}
	}

	iteger_cone_numbers = posi_iteger_cone_numbers + unposi_iteger_cone_numbers;
	std::cout << "************ Interval Round End ************" << endl;
}

void Prune::single_round(int n, int rest_cones)
{
	std::cout << endl;
	std::cout << "************ Single Round Start ************" << endl;
	int loop_num = 0;
	vector<pair<int, double>> cone_pair;
	while (rest_cones > 10)
	{
		for (int y = 0; y < mesh.n_vertices(); y++)
		{
			if (abs(rho[y]) < 1 && abs(rho[y]) > 0)
			{
				cone_pair.push_back(make_pair(y, abs(rho[y] - std::round(rho[y]))));
			}
		}
		sort(cone_pair.begin(), cone_pair.end(), cmp_bigger);

		rho[cone_pair[0].first] = std::round(rho[cone_pair[0].first]);

		Model::t M = new Model("cqo1"); auto _M = finally([&]() { M->dispose(); });
		Variable::t x = M->variable("x", 2 * n + 1, Domain::unbounded());
		M->constraint(x->index(0), Domain::greaterThan(0));

		for (int i = 0; i < mesh.n_vertices(); i++)
		{
			if (rho[i] == 0)
			{
				M->constraint(x->index(n + 1 + i), Domain::equalsTo(0));
			}
			else if (rho[i] == 1)
			{
				M->constraint(x->index(n + 1 + i), Domain::equalsTo(1));
			}
			else if (rho[i] == -1)
			{
				M->constraint(x->index(n + 1 + i), Domain::equalsTo(-1));
			}
			else
			{
				M->constraint(x->index(n + 1 + i), Domain::inRange(-1, 1));
			}
		}

		auto rows = monty::new_array_ptr<int>(mat_rows);
		auto cols = monty::new_array_ptr<int>(mat_cols);
		auto values = monty::new_array_ptr<double>(mat_vals);
		auto A = mosek::fusion::Matrix::sparse(n, 2 * n + 1, rows, cols, values);
		auto B = monty::new_array_ptr<double>(K_0);

		M->constraint(Expr::add(Expr::mul(A, x), B), Domain::equalsTo(0.0));

		Variable::t z1 = Var::vstack(x->index(0), x->slice(1, n + 1));
		Constraint::t qc1 = M->constraint("qc1", z1, Domain::inQCone());

		M->objective("obj", ObjectiveSense::Minimize, x->index(0));
		M->solve();

		ndarray<double, 1> S1 = *((x->slice(1, n + 1))->level());
		ndarray<double, 1> rho_array = *((x->slice(n + 1, 2 * n + 1))->level());

		for (int i = 0; i < n; i++)
		{
			rho[i] = rho_array(i);
		}
		loop_num++;
		double dis_1 = sqrt(mosek_norm_2(S1));
		cone_pair.clear();
		rest_cones--;
		std::cout << "rest " << rest_cones << " noninteger  distortion:  " << sqrt(mosek_norm_2(S1)) << endl;
	}
	non_iteger_cone_numbers = rest_cones;
	std::cout << "************ Single Round End ************" << endl;
	//cout << rest_cones << endl;
}

void Prune::enum_round(int enum_num)
{
	std::cout << endl;
	std::cout << "************ Enum Round Start ************ "<<endl; 
	std::cout<<"noninteger number : " << enum_num << endl;
	int posi_num = 0;
	int nega_num = 0;
	vector<int> non_iterge_cone_vec(1);
	non_iterge_cone_vec.clear();
	VectorXd ei_rho;
	ei_rho.resize(mesh.n_vertices());
	ei_rho.setZero();
	for (int t = 0; t < mesh.n_vertices(); t++)
	{
		if (rho[t] == 1)
		{
			ei_rho[t] = 1;
			posi_num++;
		}
		else if (rho[t] == -1)
		{
			ei_rho[t] = -1;
			nega_num++;
		}
		else if (rho[t] == 0)
		{

		}
		else
		{
			non_iterge_cone_vec.push_back(t);
		}
	}
	
	int diff_ = 8 - (posi_num - nega_num);
	CommonFunction cfct;
	enum_result = cfct.enumeration(enum_num, diff_);
	std::cout << "the number of enumerations£º "<<enum_result.size() << endl;
	vector<pair<int, double>> cmp_dist_num;

	solver.compute(PT* L_matrix *P);
	for (int t = 0; t < enum_result.size(); t++)
	{
		for (int t1 = 0; t1 < enum_result[t].size(); t1++)
		{
			ei_rho[non_iterge_cone_vec[t1]] = enum_result[t][t1];
		}
		double distortion_2 = slove_equation(ei_rho);
		cmp_dist_num.push_back(make_pair(t, distortion_2));
	}
	sort(cmp_dist_num.begin(), cmp_dist_num.end(), cmp_bigger);

	int increase_num = 0;
	for (int t2 = 0; t2 < enum_result[cmp_dist_num[0].first].size(); t2++)
	{
		rho[non_iterge_cone_vec[t2]] = enum_result[cmp_dist_num[0].first][t2];
		if (enum_result[cmp_dist_num[0].first][t2] != 0)
		{
			increase_num++;
		}
	}
	std::cout << "enum round distortion: " << cmp_dist_num[0].second << endl;
	cmp_dist_num.clear();

	std::cout << "************ Enum Round End ************ " << endl;
}

void Prune::round()
{

	clock_t   start, end_interval,end_round;
	start = clock();
	interval_round(mesh.n_vertices(), 10, 0, -0.9, 0, 0.1, 0, 0.9, 0);
	end_interval= clock();
	interval_round_time = (double)(end_interval - start) / CLOCKS_PER_SEC;
	std::cout << "Interval Round Time:  " << interval_round_time << "s" << endl;


	if (non_iteger_cone_numbers == 0) 
	{
		std::cout << "Round Sucess!" << endl;
		round_rho = rho;
	}
	else if (non_iteger_cone_numbers>10)
	{
		std::cout << "noninteger number£º " << non_iteger_cone_numbers << endl;
		single_round(mesh.n_vertices(), non_iteger_cone_numbers);
		enum_round(non_iteger_cone_numbers);
		round_rho = rho;
		std::cout << "Round Sucess!" << endl;
	}
	else
	{
		std::cout << "noninteger number£º " << non_iteger_cone_numbers << endl;
		enum_round(non_iteger_cone_numbers);
		round_rho = rho;
		std::cout << "Round Sucess!" << endl;
	}

	end_round= clock();
	total_round_time = (double)(end_round - start) / CLOCKS_PER_SEC;
	std::cout << "Total Round Time: " << total_round_time << "s" << endl;
}

double Prune::slove_equation(VectorXd b_right)
{
	int nv = mesh.n_vertices();
	auto B = 0.5 * M_PI * b_right - K_vec;
	VectorXd U = solver.solve(PT * B);
	Eigen::VectorXd U2;
	U2.resize(mesh.n_vertices());
	for (int k = 0; k < U.size(); k++)
	{
		U2[k] = U[k];
	}
	U2[mesh.n_vertices() - 1] = 0;

	double xishu1 = U2.transpose() * A_vec.asDiagonal() * E;
	double xishi2 = E.transpose() * A_vec.asDiagonal() * E;
	Eigen::VectorXd U1 = U2 - ((double)(xishu1 / xishi2)) * E;
	U_test = U1;
	double  distortion = sqrt(U1.transpose() * A_vec.asDiagonal() * U1);
	return distortion;
}

double Prune::mosek_norm_2(ndarray<double, 1> rho2)
{
	double n0rm_2 = 0.0;
	for (int i = 0; i < rho2.size(); i++)
	{
		n0rm_2 += std::pow(rho2[i], 2);
	}
	return n0rm_2;
}

double Prune::vec_norm_2(vector<double> rho3)
{
	double n0rm_2 = 0.0;
	for (int i = 0; i < rho3.size(); i++)
	{
		n0rm_2 += std::pow(rho3[i], 2);
	}
	return n0rm_2;
}

void Prune::prune_clear()
{
	K_0.clear();

	mat_rows.clear();
	mat_cols.clear();
	mat_vals.clear();
	L_element_num=0;
	rho.clear();

	socp_rho.clear();
	round_rho.clear();

	W_vec.clear();
	interval_rho.clear();
	non_iteger_cone_numbers=0;
	iteger_cone_numbers=0;

	enum_result.clear();
	test_socp_num = 0;

	interval_round_time = 0.0;
	total_round_time = 0.0;
}


