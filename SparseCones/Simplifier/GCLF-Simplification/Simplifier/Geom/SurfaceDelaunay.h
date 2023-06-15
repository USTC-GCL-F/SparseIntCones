#pragma once
#include "Mesh/SurfaceMeshDefinition.h"
#include "GeometryMath/GeometryMath.h"
#include "Simplifier/Topo/TopoOperations.h"
#include "Simplifier/Topo/EdgeFlipper.h"
#include "Simplifier/Topo/EdgeSplitter.h"
#include "Simplifier/Geom/GeometryCheck.h"
#include "Simplifier/SpaceSearch/FaceTree.h"
#include "Simplifier/SpaceSearch/DFaceTree.h"
#include "Simplifier/SpaceSearch/VertexTree.h"
#include "Simplifier/SpaceSearch/FaceGrid.h"
#include "Exact/Predicates.h"
#include <set>
#include <queue>

namespace GCLF
{
using namespace Geometry;
namespace Simplification
{
// TODO: extrac these functions to a higher level.
using namespace Geometry;
using namespace SMesh;
using Geometry::FaceTree;
using Geometry::DFaceTree;
using Geometry::VertexTree;

using OpenMesh::VertexHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::FaceHandle;

typedef const Vec3d& CR_Vec3d;
extern double infinite_fp;

/*******************/
/*    Delaunay     */
/*******************/

class SurfaceDelaunay
{
public:
  SMeshT* mesh;

  SurfaceDelaunay() = default;

  void initialize(SMeshT* _mesh);
  void generate();

  void check_delaunay();
public:
  SMeshT original_mesh;
	double rho_v, rho_e;
	OpenMesh::EPropHandleT<int> original_edge;

  struct EdgePriority
	{
		int e_id;
		int state;
		double priority;

		EdgePriority(int e, double p, int s) { e_id = e; priority = p; state = s; }
		bool operator> (const EdgePriority& right) const { return priority > right.priority; }
	};
	std::priority_queue<EdgePriority, std::vector<EdgePriority>, std::greater<EdgePriority>> NLD_edges;
	std::vector<int> queue_state;

  static bool test_local_delaunay_points_cot(const Vec3d& a, const Vec3d& b, const Vec3d& c);
  static bool test_local_delaunay_points_cot(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d);
  static double calc_local_delaunay_points_cot(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d);

  bool test_local_delaunay_edge_cot(EdgeHandle eh)const;
  double calc_local_delaunay_edge_cot(EdgeHandle eh)const;

  bool satisfy_Delaunay(VertexHandle vh, const Vec3d& new_point)const;
  bool satisfy_Delaunay(const std::vector<HalfedgeHandle>& halfedges, const Vec3d& new_point)const;

  Vec3d find_split_pos_make_edge_Delaunay(EdgeHandle eh)const;
  double find_pos_nearest_drop_feet(HalfedgeHandle hh)const;
  double find_drop_feet_farthest(HalfedgeHandle cur_hh)const;
  void find_split_pos_make_edge_Delaunay(const std::vector<double>& edge_len, const std::vector<double>& face_area, double& lower_bound, double& upper_bound)const;
  double find_split_pos_make_surround_edge_Delaunay(HalfedgeHandle hh, HalfedgeHandle sur_hh)const;
  double find_drop_feet(HalfedgeHandle hh)const;

  inline EdgeFlipper new_edge_flipper() { return EdgeFlipper(&original_mesh, mesh, nullptr, nullptr, nullptr); }
  inline EdgeSplitter new_edge_splitter() { return EdgeSplitter(&original_mesh, mesh, nullptr, nullptr); }
};

}// namespace Simplification
}// namespace GCLF