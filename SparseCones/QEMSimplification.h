#pragma once
#include <eigen/Eigen>
#include <queue>
#include "General_Struct.h"
class QEMSimplification
{
public:
	QEMSimplification(const Mesh & m);
	~QEMSimplification(void);
	void Simplify(int tar_num_v, double err_threshhold);
	int OriginalIndex(int currentid) const;
	const Mesh & GetMesh(void) const;
private:
	bool ChooseCollapasedEdge(double err_threshhold);
	void InitialEdgeCost(void);
	void UpdateEdgeCost(int i);
	Mesh mesh;
	std::vector<std::vector<int>> index_collapse_v_to_v; // 存储恢复时需要用到的拓扑信息
	std::vector<Mesh::Point> o_p;
	OpenMesh::VPropHandleT<int> vindex; // initial index of each vertex
	OpenMesh::FPropHandleT<double> area_f;            // area of each face
	OpenMesh::FPropHandleT<Eigen::Matrix4d> QuaMat_f; // ?? matrix of each face
	OpenMesh::VPropHandleT<Eigen::Matrix4d> QuaMat_v; // ?? matrix of each vertex
	OpenMesh::EPropHandleT<double> QuaMat_e;
	OpenMesh::EPropHandleT<OpenMesh::Vec3d> new_point; // new position after collapse
	//OpenMesh::EPropHandleT<double> edge_l;
	OpenMesh::EPropHandleT<bool> is_cal_e;
	OpenMesh::HPropHandleT<bool> is_cal_he;
	std::priority_queue<EdgeQEM> qpq;
};

