#pragma once
#include <queue>
#include <set>
#include <vector>

#include "Mesh/SurfaceMeshDefinition.h"
#include "GeometryMath/GeometryMath.h"
#include "Simplifier/Geom/GeometryCheck.h"
#include "Simplifier/Geom/SurfaceDelaunay.h"
#include "Simplifier/Topo/EdgeCollapser.h"
#include "Simplifier/Topo/EdgeFlipper.h"
#include "Simplifier/Topo/VertexRelocater.h"
#include "Simplifier/Topo/EdgeSplitter.h"
#include "GCLF_Simplification_Config.h"

namespace GCLF
{
namespace Simplification
{
namespace DE
{
using namespace Utils;
using namespace Geometry;
using namespace SMesh;

class DEDM
{
public:
  // input
  SMeshT* om;
  SMeshT* rm;
  // tree
  VertexTree* vt;
  DFaceTree* ot;
  DFaceTree* rt;
  FaceGrid* og;
  // parameters
  ParamDEDM* param;
  // -- for error bound
  double original_diagonal_length;
  double max_distance_error_value;
  Vec3d p_epsilon;
  // -- for DE
	int DE_iter_num = 100;
	double convergence_rate = 1.0e-4;
	int max_consecutive_iter = 5;

  DEDM(SMeshT* original, SMeshT* remeshing,
    VertexTree* vertex_tree, DFaceTree* original_tree, DFaceTree* remeshing_tree,
    FaceGrid* original_grid, ParamDEDM* _param):
    om(original), rm(remeshing), vt(vertex_tree), ot(original_tree), rt(remeshing_tree), og(original_grid), param(_param)
  {}

  void run();
private:
  // Delaunay refinement algorithm
  SurfaceDelaunay surf_delaunay;

  // update remeshing mesh after make it delaunay.
  void update_rm();

  // relocate DE
  BoundingBox calc_local_epsilon_bbox(SMeshT* mesh, const std::vector<FaceHandle>& faces);
  std::vector<Vec3d> generate_samples_by_links(SMeshT* mesh, const std::vector<FaceHandle>& one_ring_faces, size_t sample_num);
  std::vector<Vec3d> generate_random_samples(SMeshT* mesh, FaceHandle fh, size_t sample_num);

  void relocate_Hausdorff();
  bool predict_relocate_pos(VertexHandle vh, VertexRelocater& vertex_relocater, double max_error, Vec3d& point, std::map<Vec3d, LinkVector>& local_in_links);

  // collapse DE
  struct EdgePriority
	{
		double cost;
		int e_id;
		int state;
    EdgePriority(double _cost, int _e_id, int _state):cost(_cost), e_id(_e_id), state(_state) {}
		bool operator >(const EdgePriority& rhs) const { return cost > rhs.cost; }
	};
	std::priority_queue<EdgePriority, std::vector<EdgePriority>, std::greater<EdgePriority>> queue_cost;
	std::vector<int> queue_state;
  std::vector<Vec3d> queue_new_points;
  std::vector<std::map<Vec3d, LinkVector>> queue_new_local_in_links;

  void initialize_simplification();
  void update_cost(EdgeHandle eh);
  void update_cost(VertexHandle vh);
  bool constrained_DM_simplification(EdgeHandle eh, EdgeCollapser& edge_collapser, double max_error,
    double& cost, Vec3d& point, std::map<Vec3d, LinkVector>& local_in_links);
  void simplification();

  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, rt, og); }
  inline EdgeFlipper new_edge_flipper() { return EdgeFlipper(om, rm, ot, rt, og); }
  inline VertexRelocater new_vertex_relocater() { return VertexRelocater(om, rm, ot, rt, vt, og); }
  inline EdgeSplitter new_edge_splitter() { return EdgeSplitter(om, rm, ot, rt); }
};

}// namespace DE
}// namespace Simplification
}// namespace GCLF