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

class IntCones
{
public:

	// 1.
	IntCones(Mesh& mesh);

	// 2.
	void Init();
	void VertArea();
	void FaceAngle();
	void VertGauss();
	void MatrixL();

	// 3.
	void calc_cones_position(vector<double> socp_res, vector<double> round_res);
	double calc_mesh_distortion(Eigen::VectorXd rho_vec);
	vector<int> load_rho_int(string str);
	vector<double> load_rho_double(string str);
	void init_vinfo_noncone(vector<double>& a, vector<double>& b);
	void init_equatuion_data();
	void global_search();
	void local_search(vector<int>& a);
	void merge_and_updata_vinfo(int rings_num);

	double compute_dist(vector<double> a);
	//void ori_mesh_local_search_kai(vector<KaiCones::VInfo> vinfo, int loop_num);
	//void projection_and_search_merge(const Mesh& mesh,const Mesh& mesh_ori);

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

	void ori_mesh_local_search(vector<VInfo> vinfo, int loop_num);

public:

	Mesh& mesh;
	//Mesh& mesh_ori;

	Eigen::VectorXd f_angle;
	Eigen::VectorXd vec_A; // Areas of vertices
	Eigen::VectorXd K_ori; // Curvatures of vertices
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
	int search_rings_num = 7;// �ϲ���Χ				5K=7
	int loop_mix_num = 3;	  //						5K=3
	double ratio_sigma = 1.2;// ȫ�������Ľ�			5K=1.2
	double delta_sigma = 0.1;// �ϲ��� delta Ť��		5K=0.1
	double sigma_bound = 0.2;// Ť����					5K=0.2
	double merge_bound = 1.25;
	
	vector<double> search_merge_info;
	vector<double> ori_mesh_local_info;
	int global_search_num=0;
	VectorXd final_U;
};
//extern vector<IntCones::VInfo> vinfo;;

