#pragma once

#include "Mesh/SurfaceMeshDefinition.h"
#include "Simplifier/Topo/TopoOperations.h"
#include "Simplifier/Geom/GeometryCheck.h"

using OpenMesh::VertexHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::FaceHandle;
namespace GCLF
{
namespace Simplification
{
using namespace SMesh;
using namespace Utils;

class EdgeCollapser
{
public:
  EdgeCollapser(
    SMeshT* _origin, SMeshT* _remeshing,
    DFaceTree* _origin_tree, DFaceTree* _remeshing_tree, FaceGrid* _origin_grid)
    :om(_origin), rm(_remeshing),
    ot(_origin_tree), rt(_remeshing_tree), og(_origin_grid),
    initialized(false),
    f_update_links(false),
    f_update_target_length(false),
    f_update_normals(false),
    f_check_wrinkle(false),
    f_check_selfinter(false),
    f_check_inter(false)
  {}

  bool init(EdgeHandle _e);
  void clear();

  void set_flags(
    bool _f_update_links, bool _f_update_target_length, bool _f_update_normals,
    bool _f_update_tree, bool _f_update_one_ring_faces,
    bool _f_check_wrinkle, bool _f_check_selfinter, bool _f_check_inter
  );
  void set_local_mesh_in_links(std::map<Vec3d, LinkVector>&& local_in_links);

  bool try_collapse_edge(const Vec3d& new_point, const ExactPoint* new_ep);
  bool try_collapse_edge(EdgeHandle e, const Vec3d& new_point, const ExactPoint* new_ep);

  // some particular strategies
  bool try_collapse_avoid_long_short_edge(EdgeHandle e);
  bool try_collapse_avoid_long_short_edge_avoid_out_of_bound(EdgeHandle e, double bound);
  bool try_collapse_almost_degenerate_edge(EdgeHandle e);

  double local_Hausdorff_before_collapsing()const;
  double local_Hausdorff_after_collapsing(
    SMeshT* local_rm, const Vec3d& new_point, double& threshold,
    std::map<Vec3d, LinkVector>& local_in_links)const;

  void predict_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const;
  void predict_tangential_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const;
  void predict_weighted_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const;
  void predict_tangential_weighted_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const;

  const std::vector<HalfedgeHandle>& get_halfedges()const { ASSERT(initialized, "not initialized."); return halfedges; }
  const VertexHandle get_collapsed_center() const { ASSERT(initialized, "not initialized."); return center_vh; }
public:
  // input
  SMeshT* om;
  SMeshT* rm;

  DFaceTree* ot;
  DFaceTree* rt;
  FaceGrid* og;
  // flags
  bool f_update_links;
  bool f_update_target_length;
  bool f_update_normals;
  bool f_update_tree;
  bool f_update_one_ring_faces;
  bool f_check_wrinkle;
  bool f_check_selfinter;
  bool f_check_inter;

  void predict_faces_after_collapse();
  VertexHandle find_closer_end_point()const;

  /* safe for parallel */

  bool initialized;
  /// collapse edge
  HalfedgeHandle collapse_he;
  HalfedgeHandle collapse_he_opp;
  /// halfedges(with respect to faces) around collapse edge.
  std::vector<HalfedgeHandle> halfedges;
  /// one ring faces around collapse halfedge 
  std::set<FaceHandle> one_ring_faces;

  bool collapse_would_cause_wrinkle(const Vec3d& new_point)const;
  bool collapse_would_cause_intersection(const Vec3d& new_point, const ExactPoint* new_ep)const;
  bool collapse_would_cause_degenerate(const Vec3d& new_point, const ExactPoint* new_ep)const;

  /* only serial */

  /// in links on result mesh
  bool local_mesh_in_links_set;
  std::map<Vec3d, LinkVector> local_mesh_in_links;

  /// in links on affected faces.
  LinkVector faces_in_links;
  /// new center vertex
  VertexHandle center_vh;
  /// target length
  double new_target_length;

  void backup_links();
  void generate_links();

  void collapse_edge(HalfedgeHandle he, const Vec3d& new_point, const ExactPoint* new_ep);

  void update_target_len();
  void update_length_and_area();
  void update_normals();
  void update_remeshing_tree_deleted();
  void update_remeshing_tree_updated();
  void update_one_ring_faces();
};
}// namespace Simplification
}// namespace GCLF