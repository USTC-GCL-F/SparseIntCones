#pragma once
#include "Mesh/MeshDefinition.h"
#include <eigen/Eigen>
#include <numeric>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string>
#include "ShortestPath/Auxiliary.h"
#include "MeshCache.h"

using namespace Eigen;
using namespace std;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class Opt
{
public:

	// 1.
	Opt(Mesh& mesh);

	// 2. compute basic data
	void Init();
	void VertArea();
	void FaceAngle();
	void VertGauss();
	void MatrixL();

	// 3. 
	void calc_cones_position(vector<double> socp_res, vector<double> round_res, double p_time);
	double calc_mesh_distortion(Eigen::VectorXd rho_vec);
	void init_vinfo_noncone(vector<double>& a, vector<double>& b);
	void init_equatuion_data();
	void global_search();
	void local_search(vector<int>& a);

public:

	struct VInfo {
		int idx_i;
		double after_reweighed_rho;
		double gap;
		double after_round;
		bool is_cone = false;
	};

	struct MoveInfo {
		int before_idx_i;
		int after_idx_i;
		double distortion;

	};

	struct DeleteInfo {
		int posi_idx_i;
		int unposi_idx_i;
		double delta_dist;
		double before_dele_dist;
		double after_dele_dist;
	};

	static bool cmp_smaller(const pair<int, double>& a, const pair<int, double>& b) {
		return a.second > b.second;
	}

	static bool cmp_bigger(const pair<int, double>& a, const pair<int, double>& b) {
		return a.second < b.second;
	}

	static bool com_smaller_delta_dist(const DeleteInfo& a, const DeleteInfo& b) {
		return a.delta_dist > b.delta_dist;
	}

	static bool com_smaller_gap(const VInfo& a, const VInfo& b) {
		return a.gap > b.gap;
	}

	static bool com_bigger_gap(const VInfo& a, const VInfo& b) {
		return a.gap < b.gap;
	}

	static bool com_bigger_idx(const VInfo& a, const VInfo& b) {
		return a.idx_i < b.idx_i;
	}

	void ori_mesh_local_search(vector<VInfo> vinfo, int loop_num);
	double local_search2(VectorXd merge_rho,vector<VInfo> vv_info);
	void merge_cones_pair(int rings_num);

public:

	Mesh& mesh;

	Eigen::VectorXd f_angle;
	Eigen::VectorXd vec_A; 
	Eigen::VectorXd K_ori;
	SpMat m_L;
	Eigen::SimplicialLDLT<SpMat> solver;
	Eigen::VectorXd rho_vec;
	double rj = 0.0;
	SpMat P;
	SpMat PT;
	Eigen::VectorXd U;
	Eigen::VectorXd B;
	Eigen::VectorXd E;
	vector<VInfo> vinfo;
	
	vector<int> non_cone;
	vector<pair<int, double>> compare;

	double distortion;
	int search_rings_num = 10;
	int loop_mix_num = 3;	  
	double ratio_sigma = 1.2;
	double delta_sigma = 0.1;
	double sigma_bound = 0.2;				
	double merge_bound = 1;

	vector<double> search_merge_info;
	vector<double> ori_mesh_local_info;
	int global_search_num = 0;
	double distortion_last = 0;
	Eigen::VectorXd final_U;
	int global_search_num_each_loop = 0;

	double ori_mesh_distortion = 0.0;
	VectorXd judge_merge_rho_vec;
	vector<vector<VInfo>> vec_vec_vinfo;
};


