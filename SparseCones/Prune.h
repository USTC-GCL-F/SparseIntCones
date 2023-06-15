#pragma once
#include <iostream>
#include <fstream>
#include <fusion.h>
#include <monty.h>
#include <eigen/Eigen>
#include "CommonFunction.h"
#include "Mesh/MeshDefinition.h"
#include "Mesh/MeshData.h"

using namespace std;
using namespace Eigen;
using namespace mosek::fusion;
using namespace monty;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class Prune
{
public:
	Prune(Mesh& mesh);
	void vert_area();
	void face_angle();
	void vert_gauss();
	void calc_laplace();
	void generate_mesh_info();
	void socp(int vertex_num, double sigma, int loop_num_socp);
	void interval_round(int n, int loop_num, double delta, double left, double left_sigma, double mid, double mid_sigma, double right, double right_sigma);
	void single_round(int n, int non_iteger_cone_number);
	void enum_round(int enum_num);
	void round();
	double slove_equation(VectorXd b_right);

	double mosek_norm_2(ndarray<double, 1> rho2);
	double vec_norm_2(vector<double> rho3);

	void prune_clear();

	static bool cmp_bigger(const pair<int, double>& a, const pair<int, double>& b) {
		return a.second < b.second;
	}

public:
	Mesh& mesh;
	Eigen::VectorXd f_angle;
	SpMat L_matrix;
	SpMat P;
	SpMat PT;
	VectorXd A_vec;
	VectorXd K_vec;
	vector<double> K_0;

	std::vector<int> mat_rows;
	std::vector<int> mat_cols;
	std::vector<double> mat_vals;
	int L_element_num;
	vector<double> rho;

	vector<double> socp_rho;
	vector<double> round_rho;

	vector<double> W_vec;
	vector<double> interval_rho;
	int non_iteger_cone_numbers;
	int iteger_cone_numbers;

	VectorXd E;
	vector<vector<int>> enum_result;
	Eigen::SimplicialLDLT<SpMat> solver;
	int test_socp_num = 0;
	VectorXd U_test;

	double interval_round_time = 0.0;
	double total_round_time = 0.0;
};

