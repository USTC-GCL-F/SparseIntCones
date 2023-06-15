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

class VertexRelocater
{
public:
  VertexRelocater(
    SMeshT* _origin, SMeshT* _remeshing,
    DFaceTree* _origin_tree, DFaceTree* _remeshing_tree,
    VertexTree* _vertex_tree,
    FaceGrid* _origin_grid)
    :om(_origin), rm(_remeshing),
    ot(_origin_tree), rt(_remeshing_tree),
    vt(_vertex_tree), og(_origin_grid),
    initialized(false),
    f_update_links(false),
    f_update_target_length(false),
    f_update_normals(false),
    f_check_wrinkle(false),
    f_check_intersection(false),
    f_check_self_intersection(false)
  {}

  bool init(VertexHandle _v);
  void clear();

  void set_flags(
    bool _f_update_links, bool _f_update_target_length, bool _f_update_normals, bool _f_update_tree,
    bool _f_check_wrinkle, bool _f_check_intersection, bool _f_check_self_intersection
  );
  void set_new_target_length(double ntl) { new_target_length = ntl; }
  void set_local_mesh_in_links(std::map<Vec3d, LinkVector>&& local_in_links);

  const std::vector<HalfedgeHandle>& get_halfedges()
  {
    ASSERT(initialized, "relocater not initialized");
    return halfedges;
  }

  void relocate(const Vec3d& new_point);
  bool try_relocate_vertex(const Vec3d& new_point);
  bool try_relocate_vertex(VertexHandle _v, const Vec3d& new_point);

  // some particular strategies
  bool try_reloate_vertex_avoid_out_of_bound(const Vec3d& new_point, double bound);

  double local_Hausdorff_before_relocating()const;
  double local_Hausdorff_after_relocating(
    SMeshT* local_rm, const Vec3d& new_point, double threshold,
    std::map<Vec3d, LinkVector>& local_in_links)const;

  Vec3d find_smooth_target()const;
  Vec3d find_tangential_smooth_target()const;
  Vec3d find_weighted_smooth_target()const;
  Vec3d find_weighted_tangential_smooth_target()const;
  Vec3d find_smooth_target_with_target_length()const;
public:
  // input
  SMeshT* om;
  SMeshT* rm;

  VertexTree* vt;
  DFaceTree* ot;
  DFaceTree* rt;
  FaceGrid* og;
  // flags
  bool f_update_links;
  bool f_update_target_length;
  bool f_update_normals;
  bool f_update_tree;
  bool f_check_wrinkle;
  bool f_check_intersection;
  bool f_check_self_intersection;

  double new_target_length;

  void find_faces_affected();

  /* safe for parallel */

  bool initialized;
  /// relocate vertex
  VertexHandle relocate_vh;
  /// halfedges(with respect to faces) opposite to vertex.
  std::vector<HalfedgeHandle> halfedges;
  /// one ring faces around vertex.
  std::set<FaceHandle> one_ring_faces;

  bool relocate_would_cause_wrinkle(const Vec3d& new_point)const;
  bool relocate_would_cause_intersection(const Vec3d& new_point)const;
  bool relocate_would_cause_degenerate(const Vec3d& new_point)const;

  /* only serial */

  /// in links on affected faces.
  LinkVector faces_in_links;

  bool local_mesh_in_links_set;
  std::map<Vec3d, LinkVector> local_mesh_in_links;

  void backup_links();
  void generate_links();

  void update_target_len();
  void update_length_and_area();
  void update_remeshing_tree();
  void update_normals();
};
}// namespace Simplification
}// namespace GCLF