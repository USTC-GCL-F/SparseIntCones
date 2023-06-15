#include "EdgeSplitter.h"

namespace GCLF
{
namespace Simplification
{

void EdgeSplitter::set_flags(
  bool _f_update_links, bool _f_update_target_length, bool _f_update_normals,
  bool _f_update_tree, bool _f_update_one_ring_faces)
{
  f_update_links = _f_update_links;
  f_update_target_length = _f_update_target_length;
  f_update_normals = _f_update_normals;
  f_update_tree = _f_update_tree;
  f_update_one_ring_faces = _f_update_one_ring_faces;
}

/// @brief split edge at midpoint, 
VertexHandle EdgeSplitter::split(EdgeHandle e)
{
  return split(e, rm->calc_edge_midpoint(e));
}

VertexHandle EdgeSplitter::split(EdgeHandle e, const Vec3d& _split_point)
{
  // notations come from "split" source code of OpenMesh.
  h0 = rm->halfedge_handle(e, 0);

  // pre-update 
  if (f_update_target_length)
  {
    VertexHandle v_from = rm->from_vertex_handle(h0);
    VertexHandle v_to = rm->to_vertex_handle(h0);
    // calculate new target length
    new_target_length = 0.5 * (
      rm->data(v_from).target_length + rm->data(v_to).target_length);
  }

  // split at midpoint
  split_point = _split_point;
  split_v = rm->split(e, split_point);
  set_handles_after_split();

  // post-update
  update_length_and_area();
  if (f_update_target_length)
    update_target_len();
  if (f_update_normals)
    update_normals();
  if (f_update_tree)
    update_tree();
  if (f_update_one_ring_faces)
    update_one_ring_faces();
  if (f_update_links)
    split_in_links();
  return split_v;
}

void EdgeSplitter::set_handles_after_split()
{
  if (!rm->is_boundary(h0))
  {
    f0 = rm->face_handle(h0);
    t0 = rm->opposite_halfedge_handle(rm->prev_halfedge_handle(h0));
    f1 = rm->face_handle(t0);
  }

  o0 = rm->opposite_halfedge_handle(h0);
  if (!rm->is_boundary(o0))
  {
    f3 = rm->face_handle(o0);
    t2 = rm->opposite_halfedge_handle(rm->next_halfedge_handle(o0));
    f2 = rm->face_handle(t2);
  }
}

/// @note notations come from source code of OpenMesh
void EdgeSplitter::update_tree()
{
  // update dynamic aabb tree
  if (!rm->is_boundary(h0))
  {
    // update f0
    rt->update(rm, f0);
    // insert f1
    rt->insert(rm, f1);
  }
  if (!rm->is_boundary(o0))
  {
    // update f3
    rt->update(rm, f3);
    // inert f2
    rt->insert(rm, f2);
  }
}

void EdgeSplitter::update_target_len()
{
  rm->data(split_v).target_length = new_target_length;
  for (HalfedgeHandle voh : rm->voh_range(split_v))
  {
    rm->data(rm->edge_handle(voh)).target_length =
      std::min(
        rm->data(rm->to_vertex_handle(voh)).target_length,
        new_target_length
      );
  }
}

void EdgeSplitter::update_length_and_area()
{
  for (EdgeHandle ve : rm->ve_range(split_v))
    rm->data(ve).edge_length = rm->calc_edge_length(ve);

  if (!rm->is_boundary(h0))
  {
    rm->data(f0).face_area = calc_face_area(rm, f0);
    rm->data(f1).face_area = calc_face_area(rm, f1);
  }
  if (!rm->is_boundary(o0))
  {
    rm->data(f2).face_area = calc_face_area(rm, f2);
    rm->data(f3).face_area = calc_face_area(rm, f3);
  }
}

void EdgeSplitter::update_one_ring_faces()
{
  // erase f0 f3 from some faces.
  if (!rm->is_boundary(h0))
  {
    for (FaceHandle vfh : rm->vf_range(rm->to_vertex_handle(o0)))
      rm->data(vfh).one_ring_faces.erase(f0);
  }
  if (!rm->is_boundary(o0))
  {
    for (FaceHandle vfh : rm->vf_range(rm->to_vertex_handle(o0)))
      rm->data(vfh).one_ring_faces.erase(f3);
  }
  // insert f0 f1 f2 f3 to some faces.
  if (!rm->is_boundary(h0))
  {
    for (FaceHandle vfh : rm->vf_range(rm->to_vertex_handle(o0)))
    {
      rm->data(vfh).one_ring_faces.insert(f1);
    }
    for (FaceHandle vfh : rm->vf_range(rm->opposite_vh(h0)))
    {
      rm->data(vfh).one_ring_faces.insert(f0);
      rm->data(vfh).one_ring_faces.insert(f1);
    }
  }
  if (!rm->is_boundary(o0))
  {
    for (FaceHandle vfh : rm->vf_range(rm->to_vertex_handle(o0)))
    {
      rm->data(vfh).one_ring_faces.insert(f2);
    }
    for (FaceHandle vfh : rm->vf_range(rm->opposite_vh(o0)))
    {
      rm->data(vfh).one_ring_faces.insert(f2);
      rm->data(vfh).one_ring_faces.insert(f3);
    }
  }

  if (!rm->is_boundary(h0))
  {
    init_one_ring_faces(rm, f0);
    init_one_ring_faces(rm, f1);
  }
  if (!rm->is_boundary(o0))
  {
    init_one_ring_faces(rm, f2);
    init_one_ring_faces(rm, f3);
  }
}

void EdgeSplitter::update_normals()
{
  if (!rm->is_boundary(h0))
    rm->update_normal(f1);
  if (!rm->is_boundary(o0))
    rm->update_normal(f2);
  rm->update_normal(split_v);
}

/// @brief Split face_in_links on original f0 on new f0 and f1.
/// Split face_in_links on original f3 on new f3 and f2.
void EdgeSplitter::split_in_links()
{
  auto& new_p = rm->point(rm->to_vertex_handle(o0));
  if (!rm->is_boundary(h0) && !rm->data(rm->face_handle(h0)).face_in_links.empty())
  {
    // split links on original f0 to new f0 and f1.
    // use the orientation of link.second with respect to middle line.
    // (1) find middle line
    auto& right_p = rm->point(rm->to_vertex_handle(h0));
    auto& middle_p = rm->point(rm->to_vertex_handle(t0));

    int right_ori = coplanar_orient2d(new_p, middle_p, right_p);
    // (2) split links
    auto& f0_links = rm->data(f0).face_in_links;
    auto& f1_links = rm->data(f1).face_in_links;  f1_links.reserve(f0_links.size() / 2);
    std::vector<uint8_t> is_deleted(f0_links.size(), false);
    for (size_t i = 0;i < f0_links.size();i++)
    {
      if (coplanar_orient2d(new_p, middle_p, f0_links[i].second) != right_ori)
      {
        is_deleted[i] = true;
        f1_links.push_back(f0_links[i]);
      }
    }
    // (3) really delete deleted links.
    size_t front = 0, back = f0_links.size() - 1;
    while (true)
    {
      while (!is_deleted[front] && front < back)front++;
      while (is_deleted[back] && front < back)back--;
      if (front >= back)break;

      f0_links[front] = f0_links[back];
      is_deleted[front] = false;
      is_deleted[back] = true;
    }
    f0_links.resize(is_deleted[front] ? front : front + 1);
  }
  if (!rm->is_boundary(o0) && !rm->data(rm->face_handle(o0)).face_in_links.empty())
  {
    // split links on original f3 to new f3 and f2.
    // use the orientation of link.second with respect to middle line.
    // (1) find middle line
    auto& right_p = rm->point(rm->to_vertex_handle(h0));
    auto& middle_p = rm->point(rm->from_vertex_handle(t2));

    int right_ori = coplanar_orient2d(new_p, middle_p, right_p);
    // (2) split links
    auto& f3_links = rm->data(f3).face_in_links;
    auto& f2_links = rm->data(f2).face_in_links;  f2_links.reserve(f3_links.size() / 2);
    std::vector<uint8_t> is_deleted(f3_links.size(), false);
    for (size_t i = 0;i < f3_links.size();i++)
    {
      if (coplanar_orient2d(new_p, middle_p, f3_links[i].second) != right_ori)
      {
        is_deleted[i] = true;
        f2_links.push_back(f3_links[i]);
      }
    }
    // (3) really delete deleted links.
    size_t front = 0, back = f3_links.size() - 1;
    while (true)
    {
      while (!is_deleted[front] && front < back)front++;
      while (is_deleted[back] && front < back)back--;
      if (front >= back)break;

      f3_links[front] = f3_links[back];
      is_deleted[front] = false;
      is_deleted[back] = true;
    }
    f3_links.resize(is_deleted[front] ? front : front + 1);
  }
}
}// namespace Simplification
}// namespace GCLF