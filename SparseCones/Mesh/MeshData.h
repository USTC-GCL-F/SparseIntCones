#pragma once
#include "MeshDefinition.h"
#include <eigen/eigen>

using namespace Eigen;
using namespace std;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class MeshData
{
public:
	MeshData(Mesh& mesh);
	void Init();
	void VertArea();
	void FaceAngle();
	void VertGauss();
	void MatrixL(); 
	void MatrixP();
public:
	Mesh& mesh;
	Eigen::VectorXd f_angle;
	SpMat L_matrix;
	SpMat P;
	SpMat PT;
	VectorXd A_vec;
	VectorXd K_vec;
	vector<double> socp_info;
	
};