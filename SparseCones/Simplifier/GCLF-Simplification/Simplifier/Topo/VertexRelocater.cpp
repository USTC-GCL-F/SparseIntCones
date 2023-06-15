#include "VertexRelocater.h"

namespace GCLF
{
namespace Simplification
{

bool VertexRelocater::init(VertexHandle _v)
{
  clear();
  relocate_vh = _v;

  find_faces_affected();

  initialized = true;

  return true;
}

void VertexRelocater::clear()
{
  halfedges.clear();
  one_ring_faces.clear();
  faces_in_links.clear();
  local_mesh_in_links.clear();
  initialized = false;
  local_mesh_in_links_set = false;
}

void VertexRelocater::set_flags(
  bool _f_update_links, bool _f_update_target_length, bool _f_update_normals, bool _f_update_tree,
  bool _f_check_wrinkle, bool _f_check_intersection, bool _f_check_self_intersection)
{
  f_update_links = _f_update_links;
  f_update_target_length = _f_update_target_length;
  f_update_normals = _f_update_normals;
  f_update_tree = _f_update_tree;
  f_check_wrinkle = _f_check_wrinkle;
  f_check_intersection = _f_check_intersection;
  f_check_self_intersection = _f_check_self_intersection;
}

void VertexRelocater::set_local_mesh_in_links(std::map<Vec3d, LinkVector>&& local_in_links)
{
  ASSERT(!local_mesh_in_links_set, "duplicate local in links.");
  local_mesh_in_links = std::move(local_in_links);
  local_mesh_in_links_set = true;
}

void VertexRelocater::relocate(const Vec3d& new_point)
{
  // pre-update
  if (f_update_links)
    backup_links();

  // real relocating
  rm->set_point(relocate_vh, new_point);

  // post-update
  update_length_and_area();
  if (f_update_target_length)
    update_target_len();
  if (f_update_normals)
    update_normals();
  if (f_update_tree)
    update_remeshing_tree();
  if (f_update_links)
    generate_links();
}

bool VertexRelocater::try_relocate_vertex(const Vec3d& new_point)
{
  ASSERT(initialized, "unintialized vertex relocater.");
  // check constraints
  if (relocate_would_cause_degenerate(new_point))
    return false;
  if (relocate_would_cause_intersection(new_point))
    return false;
  if (f_check_wrinkle && relocate_would_cause_wrinkle(new_point))
    return false;

  relocate(new_point);

  return true;
}

/// @brief A simple wrap of try_relocate_vertex.
bool VertexRelocater::try_relocate_vertex(VertexHandle _v, const Vec3d& new_point)
{
  if (init(_v))
    return try_relocate_vertex(new_point);
  else
    return false;
}

bool VertexRelocater::try_reloate_vertex_avoid_out_of_bound(const Vec3d& new_point, double bound)
{
  ASSERT(initialized, "unintialized vertex relocater.");

  OpenMesh::VertexHandle local_center_v;
  std::shared_ptr<SMeshT> local_mesh = construct_local_mesh(rm, get_halfedges(), new_point, local_center_v, /*is_disk*/!rm->is_boundary(relocate_vh));
  std::map<Vec3d, LinkVector> local_in_links;

  // check new point will violate error bounds
  if (local_Hausdorff_after_relocating(local_mesh.get(), new_point, bound, local_in_links) != DBL_MAX)
  {
    set_local_mesh_in_links(std::move(local_in_links));
    relocate(new_point);
    return true;
  }
  else return false;
}

double VertexRelocater::local_Hausdorff_before_relocating()const
{
  ASSERT(initialized, "relocater not initialized.");

  return rm->data(relocate_vh).out_surround_error;
}

/// @brief Assume relocater is initialized, calculate local Hausdorff distance.
/// @param [in] local_rm local remeshing mesh, neighbor faces of vertex.
/// @param [in] new_point new point of relocating vertex.
/// @param [in] threshold hausdorff distance threshold.
/// @param [in] local_in_links In links for local mesh, reused to update rm after relocating.
/// @return false if violate constraints.
double VertexRelocater::local_Hausdorff_after_relocating(
  SMeshT* local_rm, const Vec3d& new_point, double threshold,
  std::map<Vec3d, LinkVector>& local_in_links)const
{
  ASSERT(initialized, "relocater not initialized.");

  if (relocate_would_cause_degenerate(new_point))
    return DBL_MAX;

  if (f_check_wrinkle && relocate_would_cause_wrinkle(new_point))
    return DBL_MAX;

  if (relocate_would_cause_intersection(new_point))
    return DBL_MAX;

  double hd = 0.0;

  // 0. reserve memory for store new in_links
  size_t n_links = 0;
  for (FaceHandle f : one_ring_faces)
      n_links += rm->data(f).face_in_links.size();
  n_links /= local_rm->n_faces();
  for (FaceHandle f : local_rm->faces())
    local_rm->data(f).face_in_links.reserve(n_links);
  // 1. construct local tree, calculate "in" hausdorff distance.
  FaceTree local_rt(*local_rm);
  for (FaceHandle f : one_ring_faces)
  {
    if (!local_rt.pre_hint_set && !rm->data(f).face_in_links.empty())
      local_rt.set_hint(rm->data(f).face_in_links[0].first);
    for (Link& link : rm->data(f).face_in_links)
    {
      std::pair<Vec3d, FaceHandle> cp = local_rt.closest_point_and_face_handle(link.first);
      double cd = (cp.first - link.first).norm();
      if (cd > threshold)
        return DBL_MAX;
      if (cd > hd)
        hd = cd;
      local_rm->data(cp.second).face_in_links.emplace_back(link.first, cp.first);
    }
  }
  // 2. calculate "out" hausdorff distance.
#ifdef USE_TREE_SEARCH
  ot->set_hint(new_point);
  hd = std::max(hd, calc_out_error(local_rm, om, ot, threshold));
  ot->unset_hint();
#else
  hd = std::max(hd, calc_out_error(local_rm, om, og, threshold));
#endif
  // 3. store face_in_links into map.
  if (hd != DBL_MAX)
  {
    // it's possible we finally adopt the new_point and relocate the mesh.
    local_in_links.clear();
    for (auto voh : local_rm->voh_range(VertexHandle(0)))
    {
      if (!local_rm->is_boundary(voh))
      {
        local_in_links[local_rm->point(local_rm->to_vertex_handle(voh))] =
          std::move(local_rm->data(local_rm->face_handle(voh)).face_in_links);
      }
    }
  }
  return hd;
}

Vec3d VertexRelocater::find_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  Vec3d q(0.0, 0.0, 0.0);
  for (VertexHandle vv : rm->vv_range(relocate_vh))
  {
    q += rm->point(vv);
  }
  q /= rm->valence(relocate_vh);
  return q;
}

/// @brief Uniform Laplacian smoothing.
Vec3d VertexRelocater::find_tangential_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  Vec3d q(0.0, 0.0, 0.0);
  for (VertexHandle vv : rm->vv_range(relocate_vh))
  {
    q += rm->point(vv);
  }
  q /= rm->valence(relocate_vh);
  return q + (rm->normal(relocate_vh) | (rm->point(relocate_vh) - q)) * rm->normal(relocate_vh);
}

Vec3d VertexRelocater::find_weighted_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  const Vec3d& vertex_normal = rm->normal(relocate_vh);
  Vec3d q(0., 0., 0.);
  double weight = 0.0;

  for (FaceHandle f : rm->vf_range(relocate_vh))
  {
    Vec3d barycenter(0., 0., 0.);
    for (VertexHandle fv : rm->fv_range(f))
      barycenter += rm->point(fv);
    barycenter /= 3.0;
    q += barycenter * rm->data(f).face_area;
    weight += rm->data(f).face_area;
  }
  if (weight != 0.0)
  {
    q /= weight;
    return q;
  }
  else return find_smooth_target();
}

Vec3d VertexRelocater::find_weighted_tangential_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  const Vec3d& vertex_normal = rm->normal(relocate_vh);
  Vec3d q(0., 0., 0.);
  double weight = 0.0;

  for (FaceHandle f : rm->vf_range(relocate_vh))
  {
    Vec3d barycenter(0., 0., 0.);
    for (VertexHandle fv : rm->fv_range(f))
      barycenter += rm->point(fv);
    barycenter /= 3.0;
    q += barycenter * rm->data(f).face_area;
    weight += rm->data(f).face_area;
  }
  if (weight != 0.0)
  {
    q /= weight;
    return q + (vertex_normal | (rm->point(relocate_vh) - q)) * vertex_normal;
  }
  else return find_tangential_smooth_target();
}

/// @brief tangential relaxation in "Isotropic Surface Remeshing without Large
/// and Small Angles".
Vec3d VertexRelocater::find_smooth_target_with_target_length()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  Vec3d q(0.0, 0.0, 0.0);
  double w = 0.0;
  for (FaceHandle fh : rm->vf_range(relocate_vh))
  {
    double t_i = rm->data(fh).face_area;
    double L_b_i = 0.0;
    Vec3d b_i(0.0, 0.0, 0.0);
    for (VertexHandle vv : rm->fv_range(fh))
    {
      b_i += rm->point(vv);
      L_b_i += rm->data(vv).target_length;
    }
    b_i /= 3.0;
    L_b_i /= 3.0;
    w += (t_i * L_b_i);
    q += (t_i * L_b_i) * b_i;
  }
  if (w != 0.0)
  {
    q /= w;
    return q + (rm->normal(relocate_vh) | (rm->point(relocate_vh) - q)) * rm->normal(relocate_vh);
  }
  else return find_tangential_smooth_target();
}

void VertexRelocater::update_normals()
{
  for (FaceHandle vf : rm->vf_range(relocate_vh))
    rm->update_normal(vf);
  rm->update_normal(relocate_vh);
  for (VertexHandle vv : rm->vv_range(relocate_vh))
    rm->update_normal(vv);
}

void VertexRelocater::find_faces_affected()
{
  halfedges.reserve(rm->valence(relocate_vh));
  for (HalfedgeHandle voh : rm->voh_range(relocate_vh))
  {
    if (rm->is_boundary(voh))
      continue;
    halfedges.push_back(rm->next_halfedge_handle(voh));
    one_ring_faces.insert(rm->face_handle(voh));
  }
}

void VertexRelocater::backup_links()
{
  if (!local_mesh_in_links_set)
  {
    faces_in_links = backup_local_in_links(
      rm, std::vector<FaceHandle>(one_ring_faces.begin(), one_ring_faces.end()));
  }
  // else we use local_mesh_in_links
}

void VertexRelocater::generate_links()
{
  if (local_mesh_in_links_set)
  {
    for (OpenMesh::HalfedgeHandle voh : rm->voh_range(relocate_vh))
    {
      if (!rm->is_boundary(voh))
      {
        rm->data(rm->face_handle(voh)).face_in_links =
          std::move(local_mesh_in_links[rm->point(rm->to_vertex_handle(voh))]);
      }
    }
  }
  else
  {
    FaceTree local_face_tree(*rm, std::vector<FaceHandle>(one_ring_faces.begin(), one_ring_faces.end()));

    for (Link& link : faces_in_links)
    {
      auto cp = local_face_tree.closest_point_and_face_handle(link.first);
      link.second = cp.first;
      rm->data(cp.second).face_in_links.push_back(link);
    }
  }
}

bool VertexRelocater::relocate_would_cause_wrinkle(const Vec3d& new_point)const
{
  return check_wrinkle(rm, halfedges, new_point);
}

bool VertexRelocater::relocate_would_cause_intersection(const Vec3d& new_point)const
{
  return check_intersection(
    rm, halfedges, new_point, nullptr,
    ot, og, rt, one_ring_faces,
    f_check_self_intersection, f_check_intersection);
}

bool VertexRelocater::relocate_would_cause_degenerate(const Vec3d& new_point)const
{
  return check_degenerate(rm, halfedges, new_point, nullptr);
}

void VertexRelocater::update_target_len()
{
  // update target length
  rm->data(relocate_vh).target_length = new_target_length;
  for (EdgeHandle ve : rm->ve_range(relocate_vh))
  {
    rm->data(ve).target_length = std::min(new_target_length, rm->data(ve).target_length);
  }
}

void VertexRelocater::update_length_and_area()
{
  // update edge length
  for (EdgeHandle ve : rm->ve_range(relocate_vh))
    rm->data(ve).edge_length = rm->calc_edge_length(ve);
  // calculate face area after updating all edges' length.
  for (FaceHandle f : one_ring_faces)
    rm->data(f).face_area = calc_face_area(rm, f);
}

void VertexRelocater::update_remeshing_tree()
{
  for (FaceHandle f : one_ring_faces)
    rt->update(rm, f);
}

}// namespace Simplification
}// namespace GCLF