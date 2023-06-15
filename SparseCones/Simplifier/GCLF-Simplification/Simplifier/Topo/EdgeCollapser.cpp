#include "EdgeCollapser.h"
#include <set>

namespace GCLF
{
namespace Simplification
{

/// @brief initialize essential data, prepare for collapsing.
/// @return false if collapse can't be done.
bool EdgeCollapser::init(EdgeHandle _e)
{
  clear();
  collapse_he = rm->halfedge_handle(_e, 0);
  collapse_he_opp = rm->halfedge_handle(_e, 1);

  // can't collapse
  if (!rm->is_collapse_ok(collapse_he))
    return false;

  predict_faces_after_collapse();
  initialized = true;
  return true;
}

/// @brief clear all data, set collapser to uninitialized.
void EdgeCollapser::clear()
{
  halfedges.clear();
  one_ring_faces.clear();
  faces_in_links.clear();
  local_mesh_in_links.clear();
  initialized = false;
  local_mesh_in_links_set = false;
}

void EdgeCollapser::set_flags(
  bool _f_update_links, bool _f_update_target_length, bool _f_update_normals,
  bool _f_update_tree, bool _f_update_one_ring_faces,
  bool _f_check_wrinkle, bool _f_check_selfinter, bool _f_check_inter
)
{
  f_update_links = _f_update_links;
  f_update_target_length = _f_update_target_length;
  f_update_normals = _f_update_normals;
  f_update_tree = _f_update_tree;
  f_update_one_ring_faces = _f_update_one_ring_faces;
  f_check_wrinkle = _f_check_wrinkle;
  f_check_selfinter = _f_check_selfinter;
  f_check_inter = _f_check_inter;
}

void EdgeCollapser::set_local_mesh_in_links(std::map<Vec3d, LinkVector>&& local_in_links)
{
  ASSERT(!local_mesh_in_links_set, "duplicate local in links.");
  local_mesh_in_links = std::move(local_in_links);
  local_mesh_in_links_set = true;
}

/// @brief A general purpose function.
/// Try to collapse edge. It only check intersection constraint.
/// Other constraints should be checked before.
/// @return true if succeed to collapse.
bool EdgeCollapser::try_collapse_edge(const Vec3d& new_point, const ExactPoint* new_ep)
{
  ASSERT(initialized, "uninitialized edge collapser");

  // check and backup before collapsing
  if (f_check_wrinkle && collapse_would_cause_wrinkle(new_point))
    return false;
  if (collapse_would_cause_degenerate(new_point, new_ep))
    return false;
  if (collapse_would_cause_intersection(new_point, new_ep))
    return false;

  collapse_edge(collapse_he, new_point, new_ep);

  return true;
}

/// @brief A general purpose function.
/// Try to collapse edge. It check topo and intersection constraint.
/// Other constraints should be checked before.
bool EdgeCollapser::try_collapse_edge(EdgeHandle e, const Vec3d& new_point, const ExactPoint* new_ep)
{
  if (init(e))
    return try_collapse_edge(new_point, new_ep);
  else
    return false;
}

/// @brief (Used in collapse_short_edges)Check collpase would cause long edge.
bool EdgeCollapser::try_collapse_avoid_long_short_edge(EdgeHandle e)
{
  if (rm->data(e).edge_length > (0.8 * rm->data(e).target_length))
    return false; // e is not short edge

  if (!init(e))
    return false;

  // check constrained vertex
  VertexHandle from_v = rm->to_vertex_handle(collapse_he_opp);
  VertexHandle to_v = rm->to_vertex_handle(collapse_he);

  // find new point closer to original mesh than the other.
  VertexHandle closer_vh = find_closer_end_point();
  const Vec3d& new_p = rm->point(closer_vh);
  const ExactPoint* new_ep = rm->data(closer_vh).ep.get();
  new_target_length = rm->data(closer_vh).target_length;

  // check new point will cause long edge
  if (check_long_edge(rm, halfedges, new_p, new_target_length))
    return false;
  return try_collapse_edge(new_p, new_ep);
}

/// @brief (Used in collapse_short_edges in DE)Check collpase would cause long edge.
bool EdgeCollapser::try_collapse_avoid_long_short_edge_avoid_out_of_bound(EdgeHandle e, double bound)
{
  if (rm->data(e).edge_length > (0.8 * rm->data(e).target_length))
    return false; // e is not short edge

  if (!init(e))
    return false;

  VertexHandle from_v = rm->to_vertex_handle(collapse_he_opp);
  VertexHandle to_v = rm->to_vertex_handle(collapse_he);

  bool is_disk = !rm->is_boundary(e) && !rm->is_boundary(from_v) && !rm->is_boundary(to_v);
  OpenMesh::VertexHandle local_center_v;
  std::shared_ptr<SMeshT> local_mesh = construct_local_mesh(rm, get_halfedges(), rm->point(from_v), local_center_v, is_disk);
  std::map<Vec3d, LinkVector> local_in_links;

  // check new point will violate error bounds
  new_target_length = rm->data(from_v).target_length;
  if (!check_long_edge(rm, halfedges, rm->point(from_v), new_target_length)
    && local_Hausdorff_after_collapsing(local_mesh.get(), rm->point(from_v), bound, local_in_links) != DBL_MAX)
  {
    set_local_mesh_in_links(std::move(local_in_links));
    collapse_edge(collapse_he, rm->point(from_v), nullptr);
    return true;
  }
  // else we reset properties and try another point

  // check another point will violate error bounds
  local_mesh = construct_local_mesh(rm, get_halfedges(), rm->point(to_v), local_center_v, is_disk);
  new_target_length = rm->data(to_v).target_length;
  if (!check_long_edge(rm, halfedges, rm->point(to_v), new_target_length)
    && local_Hausdorff_after_collapsing(local_mesh.get(), rm->point(to_v), bound, local_in_links) != DBL_MAX)
  {
    set_local_mesh_in_links(std::move(local_in_links));
    collapse_edge(collapse_he, rm->point(to_v), nullptr);
    return true;
  }

  return false;
}

/// @brief (Used in eliminate_degenerations)collapse very short edge.
bool EdgeCollapser::try_collapse_almost_degenerate_edge(EdgeHandle e)
{
  if (!init(e))
    return false;

  // find new point closer to original mesh than the other.
  VertexHandle closer_vh = find_closer_end_point();
  if (!try_collapse_edge(rm->point(closer_vh), rm->data(closer_vh).ep.get()))
  {
    if (rm->to_vertex_handle(collapse_he) == closer_vh)
      return try_collapse_edge(
        rm->point(rm->from_vertex_handle(collapse_he)),
        rm->data(rm->from_vertex_handle(collapse_he)).ep.get());
    else
      return try_collapse_edge(
        rm->point(rm->to_vertex_handle(collapse_he)),
        rm->data(rm->to_vertex_handle(collapse_he)).ep.get());
  }
  return true;
}

void EdgeCollapser::predict_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const
{
  ASSERT(initialized, "edge collapser not initialized.");

  vertex_normal = Vec3d(0., 0., 0.);
  target = Vec3d(0., 0., 0.);
  for (HalfedgeHandle h : halfedges)
  {
    const Vec3d& from_p = rm->point(rm->from_vertex_handle(h));
    const Vec3d& to_p = rm->point(rm->to_vertex_handle(h));

    Vec3d face_normal = (from_p - new_point).cross(to_p - new_point);

    vertex_normal += face_normal.normalize();

    Vec3d barycenter = (from_p + to_p + new_point) / 3.0;
    target += barycenter;
  }
  vertex_normal.normalize();
  if (vertex_normal == Vec3d(0., 0., 0.))  vertex_normal = Vec3d(0., 0., 1.0);
}

void EdgeCollapser::predict_tangential_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const
{
  ASSERT(initialized, "edge collapser not initialized.");

  vertex_normal = Vec3d(0., 0., 0.);
  target = Vec3d(0., 0., 0.);
  for (HalfedgeHandle h : halfedges)
  {
    const Vec3d& from_p = rm->point(rm->from_vertex_handle(h));
    const Vec3d& to_p = rm->point(rm->to_vertex_handle(h));

    Vec3d face_normal = (from_p - new_point).cross(to_p - new_point);

    vertex_normal += face_normal.normalize();

    Vec3d barycenter = (from_p + to_p + new_point) / 3.0;
    target += barycenter;
  }
  vertex_normal.normalize();
  target = target + (vertex_normal | (new_point - target)) * vertex_normal;
  if (vertex_normal == Vec3d(0., 0., 0.)) vertex_normal = Vec3d(0., 0., 1.0);
}

void EdgeCollapser::predict_weighted_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const
{
  ASSERT(initialized, "edge collapser not initialized.");

  vertex_normal = Vec3d(0., 0., 0.);
  target = Vec3d(0., 0., 0.);
  double weight = 0.0;
  for (HalfedgeHandle h : halfedges)
  {
    const Vec3d& from_p = rm->point(rm->from_vertex_handle(h));
    const Vec3d& to_p = rm->point(rm->to_vertex_handle(h));

    Vec3d face_normal = (from_p - new_point).cross(to_p - new_point);
    double face_area = 0.5 * face_normal.length();
    weight += face_area;

    vertex_normal += face_normal.normalize();

    Vec3d barycenter = (from_p + to_p + new_point) / 3.0;
    target += barycenter * face_area;
  }
  vertex_normal.normalize();
  if (weight != 0.0)
    target /= weight;
  else predict_smooth_target(new_point, vertex_normal, target);
}

void EdgeCollapser::predict_tangential_weighted_smooth_target(const Vec3d& new_point, Vec3d& vertex_normal, Vec3d& target)const
{
  ASSERT(initialized, "edge collapser not initialized.");

  vertex_normal = Vec3d(0., 0., 0.);
  target = Vec3d(0., 0., 0.);
  double weight = 0.0;
  for (HalfedgeHandle h : halfedges)
  {
    const Vec3d& from_p = rm->point(rm->from_vertex_handle(h));
    const Vec3d& to_p = rm->point(rm->to_vertex_handle(h));

    Vec3d face_normal = (from_p - new_point).cross(to_p - new_point);
    double face_area = 0.5 * face_normal.length();
    weight += face_area;

    vertex_normal += face_normal.normalize();

    Vec3d barycenter = (from_p + to_p + new_point) / 3.0;
    target += barycenter * face_area;
  }
  vertex_normal.normalize();
  if (weight != 0.0)
  {
    target /= weight;
    target = target + (vertex_normal | (new_point - target)) * vertex_normal;
  }
  else predict_tangential_smooth_target(new_point, vertex_normal, target);
}

/// @brief calculate local Hausdorff distance before collapsing.
/// Assume in_error and out_error are correctly set.
double EdgeCollapser::local_Hausdorff_before_collapsing()const
{
  ASSERT(initialized, "edge collapser not initialized.");

  return std::max(rm->data(rm->to_vertex_handle(collapse_he)).out_surround_error,
    rm->data(rm->from_vertex_handle(collapse_he)).out_surround_error);
}

/// @brief Assume relocater is initialized, calculate local Hausdorff distance.
/// @param [in] new_point new point of relocating vertex.
/// @return DBL_MAX if violate constraints.
double EdgeCollapser::local_Hausdorff_after_collapsing(
  SMeshT* local_rm, const Vec3d& new_point, double& threshold,
  std::map<Vec3d, LinkVector>& local_in_links)const
{
  ASSERT(initialized, "collapser not initialized.");

  if (f_check_wrinkle && collapse_would_cause_wrinkle(new_point))
    return DBL_MAX;
  if (collapse_would_cause_degenerate(new_point, nullptr))
    return DBL_MAX;
  if (collapse_would_cause_intersection(new_point, nullptr))
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
    // it's possible we finally adopt the new_point and collapse the mesh.
    local_in_links.clear();
    for (auto vih : local_rm->vih_range(VertexHandle(0)))
    {
      if (!local_rm->is_boundary(vih))
      {
        local_in_links[local_rm->point(local_rm->from_vertex_handle(vih))] =
          std::move(local_rm->data(local_rm->face_handle(vih)).face_in_links);
      }
    }
  }
  return hd;
}

/// @brief predict new faces after collapse. new faces are represented by halfedges.
/// @return true if collapse is ok.
/// @see predict_faces_after_collapse.pdf in folder "figs"
void EdgeCollapser::predict_faces_after_collapse()
{
  int n_valance = rm->valence(rm->from_vertex_handle(collapse_he)) +
    rm->valence(rm->to_vertex_handle(collapse_he));
  halfedges.clear();
  halfedges.reserve(n_valance);

  HalfedgeHandle hh_next = rm->next_halfedge_handle(collapse_he);
  HalfedgeHandle hh_oppo_next = rm->next_halfedge_handle(collapse_he_opp);

  HalfedgeHandle moving_hh = rm->opposite_halfedge_handle(rm->prev_halfedge_handle(collapse_he));
  if (rm->is_boundary(moving_hh) && moving_hh != hh_oppo_next)
    moving_hh = rm->opposite_halfedge_handle(rm->prev_halfedge_handle(moving_hh));
  while (moving_hh != hh_oppo_next)
  {
    moving_hh = rm->next_halfedge_handle(moving_hh);
    if (!rm->is_boundary(moving_hh))
      halfedges.push_back(moving_hh);
    moving_hh = rm->opposite_halfedge_handle(rm->next_halfedge_handle(moving_hh));
    if (rm->is_boundary(moving_hh) && moving_hh != hh_oppo_next)
      moving_hh = rm->opposite_halfedge_handle(rm->prev_halfedge_handle(moving_hh));
  }

  moving_hh = rm->opposite_halfedge_handle(rm->prev_halfedge_handle(collapse_he_opp));
  if (rm->is_boundary(moving_hh) && moving_hh != hh_next)
    moving_hh = rm->opposite_halfedge_handle(rm->prev_halfedge_handle(moving_hh));
  while (moving_hh != hh_next)
  {
    moving_hh = rm->next_halfedge_handle(moving_hh);
    if (!rm->is_boundary(moving_hh))
      halfedges.push_back(moving_hh);
    moving_hh = rm->opposite_halfedge_handle(rm->next_halfedge_handle(moving_hh));
    if (rm->is_boundary(moving_hh) && moving_hh != hh_next)
      moving_hh = rm->opposite_halfedge_handle(rm->prev_halfedge_handle(moving_hh));
  }

  one_ring_faces.clear();
  for (HalfedgeHandle hh : halfedges)
    one_ring_faces.insert(rm->face_handle(hh));
  if (!rm->is_boundary(collapse_he))
    one_ring_faces.insert(rm->face_handle(collapse_he));
  if (!rm->is_boundary(collapse_he_opp))
    one_ring_faces.insert(rm->face_handle(collapse_he_opp));
}

void EdgeCollapser::backup_links()
{
  if (!local_mesh_in_links_set)
  {
    faces_in_links = backup_local_in_links(
      rm, std::vector<FaceHandle>(one_ring_faces.begin(), one_ring_faces.end()));
  }
}

void EdgeCollapser::generate_links()
{
  if (local_mesh_in_links_set)
  {
    // the iterator is different from what we used during storing, it's not a bug.
    for (OpenMesh::HalfedgeHandle voh : rm->voh_range(center_vh))
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
    std::vector<FaceHandle> faces; faces.reserve(halfedges.size());
    for (HalfedgeHandle heh : halfedges)
      faces.push_back(rm->face_handle(heh));

    FaceTree local_face_tree(*rm, faces);
    for (Link& link : faces_in_links)
    {
      auto cp = local_face_tree.closest_point_and_face_handle(link.first);
      link.second = cp.first;
      rm->data(cp.second).face_in_links.push_back(link);
    }
  }
}

VertexHandle EdgeCollapser::find_closer_end_point()const
{
  if (rm->data(rm->to_vertex_handle(collapse_he)).out_error <
    rm->data(rm->from_vertex_handle(collapse_he)).out_error)
    return rm->to_vertex_handle(collapse_he);
  else
    return rm->from_vertex_handle(collapse_he);
}

bool EdgeCollapser::collapse_would_cause_wrinkle(const Vec3d& new_point)const
{
  return check_wrinkle(rm, halfedges, new_point);
}

bool EdgeCollapser::collapse_would_cause_intersection(const Vec3d& new_point, const ExactPoint* new_ep)const
{
  return check_intersection(rm, halfedges, new_point, new_ep, ot, og, rt, one_ring_faces, f_check_selfinter, f_check_inter);
}

bool EdgeCollapser::collapse_would_cause_degenerate(const Vec3d& new_point, const ExactPoint* new_ep)const
{
  return check_degenerate(rm, halfedges, new_point, new_ep);
}

void EdgeCollapser::collapse_edge(HalfedgeHandle he, const Vec3d& new_point, const ExactPoint* new_ep)
{
  if (f_update_links)
    backup_links();

  // pre-update
  if (f_update_one_ring_faces)
    update_one_ring_faces();
  if (f_update_tree)
    update_remeshing_tree_deleted();

  // do real collapse
  center_vh = rm->to_vertex_handle(he);
  rm->set_point(center_vh, new_point);
  if (new_ep)
    rm->data(center_vh).ep = std::make_unique<ExactPoint>(*new_ep);
  else
    rm->data(center_vh).ep = nullptr;
  rm->collapse(he);

  // post-update 
  update_length_and_area();
  if (f_update_target_length)
    update_target_len();
  if (f_update_tree)
    update_remeshing_tree_updated();
  if (f_update_normals)
    update_normals();
  if (f_update_links)
    generate_links();
}

void EdgeCollapser::update_target_len()
{
  rm->data(center_vh).target_length = new_target_length;
  for (HalfedgeHandle voh : rm->voh_range(center_vh))
  {
    EdgeHandle e = rm->edge_handle(voh);
    rm->data(e).target_length = std::min(
      new_target_length,
      rm->data(rm->to_vertex_handle(voh)).target_length
    );
  }
}

void EdgeCollapser::update_length_and_area()
{
  // update edge length
  for (HalfedgeHandle hh : halfedges)
  {
    EdgeHandle next_eh = rm->edge_handle(rm->next_halfedge_handle(hh));
    rm->data(next_eh).edge_length = rm->calc_edge_length(next_eh);
  }
  // calculate face area after updating all edges' length.
  for (HalfedgeHandle hh : halfedges)
  {
    FaceHandle fh = rm->face_handle(hh);
    rm->data(fh).face_area = calc_face_area(rm, fh);
  }
}

void EdgeCollapser::update_normals()
{
  for (FaceHandle vf : rm->vf_range(center_vh))
    rm->update_normal(vf);
  rm->update_normal(center_vh);
  for (VertexHandle vv : rm->vv_range(center_vh))
    rm->update_normal(vv);
}

void EdgeCollapser::update_remeshing_tree_deleted()
{
  if (!rm->is_boundary(collapse_he))
    rt->remove(rm->face_handle(collapse_he));
  if (!rm->is_boundary(collapse_he_opp))
    rt->remove(rm->face_handle(collapse_he_opp));
}

void EdgeCollapser::update_remeshing_tree_updated()
{
  for (HalfedgeHandle heh : halfedges)
    rt->update(rm, rm->face_handle(heh));
}

void EdgeCollapser::update_one_ring_faces()
{
  FaceHandle deleted_fh = rm->face_handle(collapse_he);
  FaceHandle deleted_fh_opp = rm->face_handle(collapse_he_opp);
  // after collapsing, some faces will have more one ring faces.
  for (FaceHandle from_fh : rm->vf_range(rm->to_vertex_handle(collapse_he_opp)))
  {
    for (FaceHandle to_fh : rm->vf_range(rm->to_vertex_handle(collapse_he)))
    {
      rm->data(from_fh).one_ring_faces.insert(to_fh);
      rm->data(to_fh).one_ring_faces.insert(from_fh);
    }
  }
  // after collapsing, remove two collapsed faces from other faces' one ring.
  for (FaceHandle fh : one_ring_faces)
  {
    rm->data(fh).one_ring_faces.erase(deleted_fh);
    rm->data(fh).one_ring_faces.erase(deleted_fh_opp);
  }
  if (!rm->is_boundary(collapse_he))
  {
    for (FaceHandle vfh : rm->vf_range(rm->opposite_vh(collapse_he)))
      rm->data(vfh).one_ring_faces.erase(deleted_fh);
  }
  if (!rm->is_boundary(collapse_he_opp))
  {
    for (FaceHandle vfh : rm->vf_range(rm->opposite_vh(collapse_he_opp)))
      rm->data(vfh).one_ring_faces.erase(deleted_fh_opp);
  }
}
}// namespace Simplification
}// namespace GCLF