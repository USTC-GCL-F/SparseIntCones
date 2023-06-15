#pragma once
#include <queue>
#include <set>
#include <vector>

#include "Mesh/SurfaceMeshDefinition.h"
#include "Mesh/MeshCurvature.h"
#include "GeometryMath/GeometryMath.h"
#include "Simplifier/Geom/GeometryCheck.h"
#include "Simplifier/Topo/EdgeCollapser.h"
#include "Simplifier/Topo/EdgeFlipper.h"
#include "Simplifier/Topo/VertexRelocater.h"
#include "Simplifier/Topo/EdgeSplitter.h"
#include "GCLF_Simplification_Config.h"

namespace GCLF
{
namespace Simplification
{
using namespace Utils;
using namespace Geometry;
using namespace SMesh;
namespace DE
{

class DERemesher
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
  DEParamRemesher* param;
  // constraints
  double original_diagonal_length;
  double max_distance_error_value;
  double end_angle;

  DERemesher(SMeshT* original, SMeshT* remeshing,
    VertexTree* vertex_tree, DFaceTree* original_tree, DFaceTree* remeshing_tree,
    FaceGrid* original_grid,
    DEParamRemesher* _param):
    om(original), rm(remeshing), vt(vertex_tree), ot(original_tree), rt(remeshing_tree), og(original_grid), param(_param)
  {}

  void remesh();
private:
  /* target length */
  void calc_target_length_on_origin();
  void enlarge_target_length_on_origin();
  void find_target_length_on_remeshing();

  /* minimal angle */
  double calc_ave_min_angle(SMeshT* mesh);
  double calc_small_angle_of_tri(SMeshT* mesh, FaceHandle fh);

  /* project and target length*/
  std::pair<Vec3d, double> project_and_target_length(const Vec3d& p);

  size_t collapse_short_edges();

  size_t Laplacian_smooth();

  size_t equalize_valences();

  size_t split_long_edges();

  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, rt, og); }
  inline EdgeFlipper new_edge_flipper() { return EdgeFlipper(om, rm, ot, rt, og); }
  inline VertexRelocater new_vertex_relocater() { return VertexRelocater(om, rm, ot, rt, vt, og); }
  inline EdgeSplitter new_edge_splitter() { return EdgeSplitter(om, rm, ot, rt); }
};
}// namespace DE
}// namespace Simplification
}// namespace GCLF