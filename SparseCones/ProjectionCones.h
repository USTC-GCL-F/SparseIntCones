#pragma once
// kdtree lib
#include "IntCones.h"
#include "KdTree/KdTree.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include "Opt.h"

namespace GG = GCLF::Geometry;
using GCLF::Geometry::KdTree;

class ProjectionCones {

public:
	ProjectionCones();
	~ProjectionCones();
	

public:
	//Mesh& mesh;
	//Mesh& mesh_ori;

	struct VInfo {
		int idx_i;
		double after_reweighed_rho;
		double gap;
		double after_round;
		bool is_cone = false;
	};

	struct SimplifyInfo {
		int idx_i;
		//double after_reweighed_rho;
		//double gap;
		//double after_round;
		int ori_idx_i;
		GG::Vec3d pi;
	};

	void projection_onto_ori_mesh(const Mesh& mesh, const Mesh& mesh_ori, vector<Opt::VInfo> vinfo_projection);
	//void projection_onto_ori_mesh_kai(const Mesh& mesh, const Mesh& mesh_ori, vector<KaiCones::VInfo> vinfo_projection);
	vector<int> after_projection;
	//vector<IntCones::VInfo> vinfo_projrction;
	//vector<SimplifyInfo> sinfo_vec;
};


