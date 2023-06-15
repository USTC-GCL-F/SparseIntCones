#include "SurfaceDelaunay.h"

namespace GCLF
{
namespace Simplification
{
// handles
#define EH(hh) (mesh->edge_handle(hh))  
#define NEXT(hh) (mesh->next_halfedge_handle(hh))
#define PREV(hh) (mesh->prev_halfedge_handle(hh))
#define OPPH(hh) (mesh->opposite_halfedge_handle(hh))
#define TO(hh) (mesh->to_vertex_handle(hh))
#define FROM(hh) (mesh->from_vertex_handle(hh))
#define OPPV(hh) (mesh->opposite_vh(hh))
#define OPPH_OPPV(hh) (mesh->opposite_he_opposite_vh(hh))
#define FH(hh) (mesh->face_handle(hh))

// data
#define PNT(vh) (mesh->point(vh))
#define LEN(eh) (mesh->data(eh).edge_length)
#define AREA(fh) (mesh->data(fh).face_area)

void SurfaceDelaunay::initialize(SMeshT* _mesh)
{
  mesh = _mesh;
  original_mesh = *mesh;

  double l_min = DBL_MAX, theta_min = DBL_MAX;
  double l_temp, theta_temp;
  for (HalfedgeHandle hh : mesh->halfedges())
  {
    if (!mesh->is_boundary(hh))
    {
      l_temp = mesh->calc_edge_length(hh);
      if (l_temp < l_min)
      {
        l_min = l_temp;
      }
      theta_temp = mesh->calc_sector_angle(hh);
      if (theta_temp < theta_min)
      {
        theta_min = theta_temp;
      }
    }
  }
  rho_v = std::min((l_min * std::sin(theta_min)) / (0.5 + std::sin(theta_min)), 0.5 * l_min);
  rho_e = 2 * rho_v * sin(theta_min);

  Logger::user_logger->info("Initialize surface Delaunay. rho_v {}  rho_e {}.", rho_v, rho_e);

  mesh->add_property(original_edge, "original_edge");
  for (int i = 0; i < mesh->n_edges(); ++i)
  {
    mesh->property(original_edge, mesh->edge_handle(i)) = i;
  }
}

void SurfaceDelaunay::generate()
{
  // Initialize priority queue.
  queue_state.resize(mesh->n_edges(), 0);
  for (EdgeHandle eh : mesh->edges())
  {
    double priority = calc_local_delaunay_edge_cot(eh);
    if (priority < 0)
      NLD_edges.emplace(eh.idx(), priority, queue_state[eh.idx()]);
  }
  Logger::user_logger->info("Got {} Not-Local-Delaunay edges.", NLD_edges.size());

  // Enter loop.
  // Only use flip and split to remesh.
  auto edge_flipper = new_edge_flipper();
  edge_flipper.set_flags(
    false, false, false,
    false, false,
    false, false, false);

  auto edge_splitter = new_edge_splitter();
  edge_splitter.set_flags(
    false, false, false,
    false, false);

  std::vector<EdgeHandle> surrounding_edges; surrounding_edges.reserve(4);
  while (!NLD_edges.empty())
  {
    // get an edge to process
    EdgePriority tmp = NLD_edges.top();
    NLD_edges.pop();
    int e_process_id = tmp.e_id;
    if (tmp.state != queue_state[e_process_id])
      continue;

    EdgeHandle eh = mesh->edge_handle(e_process_id);
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0), hh_oppo = mesh->halfedge_handle(eh, 1);
    if (mesh->is_boundary(hh))
      std::swap(hh, hh_oppo);

    // get adjacent halfedges
    surrounding_edges.clear();
    surrounding_edges.push_back(mesh->edge_handle(mesh->next_halfedge_handle(hh)));
    surrounding_edges.push_back(mesh->edge_handle(mesh->prev_halfedge_handle(hh)));
    if (!mesh->is_boundary(eh))
    {
      surrounding_edges.push_back(mesh->edge_handle(mesh->next_halfedge_handle(hh_oppo)));
      surrounding_edges.push_back(mesh->edge_handle(mesh->prev_halfedge_handle(hh_oppo)));
    }

    // check if we can flip this edge
    bool flip_done = false;
    if (mesh->property(original_edge, eh) == -1)
    {
      // an edge formed by split
      flip_done = edge_flipper.try_flip_edge(eh);
    }
    else if (edge_flipper.init(eh) && adjacent_triangles_coplanar_filtered(mesh, eh) == 0)
    {
      double priority = calc_local_delaunay_points_cot(
        PNT(edge_flipper.vb1), PNT(edge_flipper.va1),
        PNT(edge_flipper.vb0), PNT(edge_flipper.va0));
      flip_done = priority >= 0 && edge_flipper.try_flip_edge();
    }

    // if we flip the edge, update edges around it.
    if (flip_done)
    {
      for (EdgeHandle sur_eh : surrounding_edges)
      {
        queue_state[sur_eh.idx()]++;
        double priority = calc_local_delaunay_edge_cot(sur_eh);
        if (priority < 0)
          NLD_edges.emplace(sur_eh.idx(), priority, queue_state[sur_eh.idx()]);
      }
      continue;
    }
    // else we split this edge to make it Delaunay
    int original_idx = mesh->property(original_edge, eh);
    // prepare to update edges
    std::set<VertexHandle> opp_vhs;
    if (mesh->is_boundary(eh))
    {
      queue_state.insert(queue_state.end(), 2, 0);
      opp_vhs.insert(OPPV(hh));
    }
    else
    {
      queue_state.insert(queue_state.end(), 3, 0);
      opp_vhs.insert(OPPV(hh));
      opp_vhs.insert(OPPH_OPPV(hh));
    }

    edge_splitter.split(eh, find_split_pos_make_edge_Delaunay(eh));

    // update edges incident to split vertex
    for (HalfedgeHandle voh : mesh->voh_range(edge_splitter.split_v))
    {
      EdgeHandle eh_voh = EH(voh);
      if (opp_vhs.count(TO(voh)))
        mesh->property(original_edge, eh_voh) = -1;
      else
      {
        queue_state[eh_voh.idx()]++;
        mesh->property(original_edge, eh_voh) = original_idx;
        double priority = calc_local_delaunay_edge_cot(eh_voh);
        if (priority < 0)
          NLD_edges.emplace(eh_voh.idx(), priority, queue_state[eh_voh.idx()]);
      }
    }

    // update edges adjacent to split edge
    for (EdgeHandle sur_eh : surrounding_edges)
    {
      queue_state[sur_eh.idx()]++;
      double priority = calc_local_delaunay_edge_cot(sur_eh);
      if (priority < 0)
        NLD_edges.push(EdgePriority(sur_eh.idx(), priority, queue_state[sur_eh.idx()]));
    }
  }
  // check
  check_delaunay();
  mesh->remove_property(original_edge);
}

void SurfaceDelaunay::check_delaunay()
{
  size_t not_local_delaunay = 0;
  for (EdgeHandle eh : mesh->edges())
    not_local_delaunay += !test_local_delaunay_edge_cot(eh);
  Logger::user_logger->info("We have {} Not-Local-Delaunay edges", not_local_delaunay);
}

bool SurfaceDelaunay::test_local_delaunay_points_cot(const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
  double a_len = (a - b).norm(), b_len = (b - c).norm(), c_len = (a - c).norm();
  double cotAlpha1 = (b_len * b_len + c_len * c_len - a_len * a_len)
    / sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));
  return cotAlpha1 >= 0;
}

bool SurfaceDelaunay::test_local_delaunay_points_cot(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d)
{
  double a_len, b_len, c_len;
  double cotAlpha1, cotAlpha2;
  // we use a approximate check for speed.
  bool is_coplanar = orient3d_filtered(a.x(), a.y(), a.z(), b.x(), b.y(), b.z(), c.x(), c.y(), c.z(), d.x(), d.y(), d.z()) == 0;

  if (is_coplanar)
  {
    a_len = (a - b).norm(); b_len = (b - c).norm(); c_len = (a - c).norm();
    cotAlpha1 = (b_len * b_len + c_len * c_len - a_len * a_len) /
      sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));
    b_len = (a - d).norm(); c_len = (b - d).norm();
    cotAlpha2 = (b_len * b_len + c_len * c_len - a_len * a_len) /
      sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));

    if (_fpclass(cotAlpha1) == _FPCLASS_PINF || _fpclass(cotAlpha2) == _FPCLASS_PINF) return true;
    if (cotAlpha1 != cotAlpha1 || cotAlpha2 != cotAlpha2) return true;
    if (cotAlpha1 + cotAlpha2 >= 0)
      return true;

    a_len = (c - d).norm(); b_len = (d - a).norm(); c_len = (c - a).norm();
    cotAlpha1 = (b_len * b_len + c_len * c_len - a_len * a_len) /
      sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));
    b_len = (c - b).norm(); c_len = (d - b).norm();
    cotAlpha2 = (b_len * b_len + c_len * c_len - a_len * a_len) /
      sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));

    if (_fpclass(cotAlpha1) == _FPCLASS_PINF || _fpclass(cotAlpha2) == _FPCLASS_PINF) return false;
    if (cotAlpha1 != cotAlpha1 || cotAlpha2 != cotAlpha2) return false;
    return cotAlpha1 + cotAlpha2 < 0;
  }
  else
  {
    a_len = (a - b).norm(); b_len = (b - c).norm(); c_len = (a - c).norm();
    cotAlpha1 = (b_len * b_len + c_len * c_len - a_len * a_len) /
      sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));
    b_len = (a - d).norm(), c_len = (b - d).norm();
    cotAlpha2 = (b_len * b_len + c_len * c_len - a_len * a_len) /
      sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));

    if (_fpclass(cotAlpha1) == _FPCLASS_PINF || _fpclass(cotAlpha2) == _FPCLASS_PINF) return true;
    if (cotAlpha1 != cotAlpha1 || cotAlpha2 != cotAlpha2) return true;
    return cotAlpha1 + cotAlpha2 >= 0;
  }
}

double SurfaceDelaunay::calc_local_delaunay_points_cot(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d)
{
  double a_len, b_len, c_len;
  double cotAlpha1, cotAlpha2;

  a_len = (a - b).norm(); b_len = (b - c).norm(); c_len = (a - c).norm();
  cotAlpha1 = (b_len * b_len + c_len * c_len - a_len * a_len) /
    sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));

  b_len = (a - d).norm(), c_len = (b - d).norm();
  cotAlpha2 = (b_len * b_len + c_len * c_len - a_len * a_len) /
    sqrt((b_len + c_len + a_len) * (b_len + c_len - a_len) * (a_len + b_len - c_len) * (a_len - b_len + c_len));
  return cotAlpha1 + cotAlpha2;
}

/// @brief After relocating vertex "vh" to new_point, check if the one ring satisfy Delaunay triangulation property.
/// @param vh center vertex
/// @param new_point put center vertex to this new position
/// @return true if satify
bool SurfaceDelaunay::satisfy_Delaunay(VertexHandle vh, const Vec3d& new_point)const
{
  for (HalfedgeHandle voh : mesh->voh_range(vh))
  {
    if (mesh->is_boundary(voh))
    {
      const Vec3d& p_b = PNT(TO(voh)),
        & p_c = PNT(OPPH_OPPV(voh));
      if (!test_local_delaunay_points_cot(new_point, p_b, p_c))
        return false;
    }
    else
    {
      const Vec3d& p_b = PNT(TO(voh)), & p_c = PNT(OPPV(voh));
      if (mesh->is_boundary(OPPH(voh)))
      {
        if (!test_local_delaunay_points_cot(new_point, p_b, p_c))
          return false;
      }
      else
      {
        const Vec3d& p_d = PNT(OPPH_OPPV(voh));
        if (!test_local_delaunay_points_cot(new_point, p_b, p_c, p_d))
          return false;
      }
      HalfedgeHandle next_voh = NEXT(voh);
      if (mesh->is_boundary(OPPH(next_voh)))
      {
        if (!test_local_delaunay_points_cot(p_b, p_c, new_point))
          return false;
      }
      else
      {
        const Vec3d& p_d = PNT(OPPH_OPPV(next_voh));
        if (!test_local_delaunay_points_cot(p_b, p_c, new_point, p_d))
          return false;
      }
    }
  }
  return true;
}

/// @brief After collapsing edge "hh" to new_point, check if the one ring satisfy Delaunay triangulation property.
/// @param halfedges one ring halfedges after collapsing
/// @param new_point put center point to this new position
/// @return true if satisfy
bool SurfaceDelaunay::satisfy_Delaunay(const std::vector<HalfedgeHandle>& halfedges, const Vec3d& new_point)const
{
  auto prev_i = [&](size_t i) { return i == 0 ? halfedges.size() - 1 : i - 1; };
  auto next_i = [&](size_t i) { return i == halfedges.size() - 1 ? 0 : i + 1; };

  for (size_t i = 0;i < halfedges.size();i++)
  {
    HalfedgeHandle hi = halfedges[i];
    HalfedgeHandle hj = halfedges[prev_i(i)];
    HalfedgeHandle hk = halfedges[next_i(i)];
    if (FROM(hi) != TO(hj))
    {
      if (!test_local_delaunay_points_cot(new_point, PNT(FROM(hi)), PNT(TO(hi))))
        return false;
    }

    if (mesh->is_boundary(OPPH(hi)))
    {
      if (!test_local_delaunay_points_cot(PNT(FROM(hi)), PNT(TO(hi)), new_point))
        return false;
    }
    else
    {
      if (!test_local_delaunay_points_cot(PNT(FROM(hi)), PNT(TO(hi)), new_point, PNT(OPPH_OPPV(hi))))
        return false;
    }

    if (TO(hi) == FROM(hk))
    {
      if (!test_local_delaunay_points_cot(new_point, PNT(TO(hi)), PNT(FROM(hi)), PNT(TO(hk))))
        return false;
    }
    else
    {
      if (!test_local_delaunay_points_cot(new_point, PNT(TO(hi)), PNT(FROM(hi))))
        return false;
    }
  }
  return true;
}

bool SurfaceDelaunay::test_local_delaunay_edge_cot(EdgeHandle eh)const
{
  if (mesh->is_boundary(eh))
  {
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
    if (mesh->is_boundary(hh))
      hh = OPPH(hh);
    double a = LEN(eh);
    double b = LEN(EH(NEXT(hh)));
    double c = LEN(EH(PREV(hh)));
    double cotAlpha1 = (b * b + c * c - a * a) /
      sqrt((b + c + a) * (b + c - a) * (a + b - c) * (a - b + c));
    return cotAlpha1 >= 0;
  }
  else
  {
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0), hh_oppo = OPPH(hh);
    double a = LEN(eh);
    double b = LEN(EH(NEXT(hh)));
    double c = LEN(EH(PREV(hh)));
    double cotAlpha1 = (b * b + c * c - a * a) /
      sqrt((b + c + a) * (b + c - a) * (a + b - c) * (a - b + c));

    b = LEN(EH(NEXT(hh_oppo)));
    c = LEN(EH(PREV(hh_oppo)));
    double cotAlpha2 = (b * b + c * c - a * a) /
      sqrt((b + c + a) * (b + c - a) * (a + b - c) * (a - b + c));

    if (_fpclass(cotAlpha1) == _FPCLASS_PINF || _fpclass(cotAlpha2) == _FPCLASS_PINF) return true;
    if (cotAlpha1 != cotAlpha1 || cotAlpha2 != cotAlpha2) return true;
    return cotAlpha1 + cotAlpha2 >= 0;
  }
}

double SurfaceDelaunay::calc_local_delaunay_edge_cot(EdgeHandle eh)const
{
  if (mesh->is_boundary(eh))
  {
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
    if (mesh->is_boundary(hh))
      hh = OPPH(hh);
    double a = LEN(eh);
    double b = LEN(EH(NEXT(hh)));
    double c = LEN(EH(PREV(hh)));
    double cotAlpha1 = (b * b + c * c - a * a) /
      sqrt((b + c + a) * (b + c - a) * (a + b - c) * (a - b + c));
    return cotAlpha1;
  }
  else
  {
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0), hh_oppo = OPPH(hh);
    double a = LEN(eh);
    double b = LEN(EH(NEXT(hh)));
    double c = LEN(EH(PREV(hh)));
    double cotAlpha1 = (b * b + c * c - a * a) /
      sqrt((b + c + a) * (b + c - a) * (a + b - c) * (a - b + c));
    b = LEN(EH(NEXT(hh_oppo)));
    c = LEN(EH(PREV(hh_oppo)));
    double cotAlpha2 = (b * b + c * c - a * a) /
      sqrt((b + c + a) * (b + c - a) * (a + b - c) * (a - b + c));
    return cotAlpha1 + cotAlpha2;
  }
}

double SurfaceDelaunay::find_drop_feet(HalfedgeHandle hh)const
{
  if (mesh->is_boundary(hh))
  {
    return DBL_MAX;
  }
  else
  {
    double len0 = LEN(EH(hh));
    double len1 = LEN(EH(NEXT(hh)));
    double len2 = LEN(EH(PREV(hh)));
    return (len0 * len0 + len2 * len2 - len1 * len1) / (2 * len0);
  }
}

/// @brief see name
/// @param mesh remeshed mesh.
/// @param original_mesh original mesh, it is the mesh when we record original_edge property.
/// @param hh halfedge
/// @return position to split
double SurfaceDelaunay::find_pos_nearest_drop_feet(HalfedgeHandle hh)const
{
  double pos = find_drop_feet(hh);
  if (mesh->is_boundary(hh) || pos<0 || pos>LEN(EH(hh)))
    return pos;

  EdgeHandle ori_eh = original_mesh.edge_handle(mesh->property(original_edge, EH(hh)));
  HalfedgeHandle ori_hh = original_mesh.halfedge_handle(ori_eh, 0);
  const Vec3d& ori_p_from = original_mesh.point(original_mesh.from_vertex_handle(ori_hh));
  const Vec3d& ori_p_to = original_mesh.point(original_mesh.to_vertex_handle(ori_hh));
  const Vec3d& a = PNT(FROM(hh));
  const Vec3d& b = PNT(TO(hh));
  double len1 = (ori_p_from - a).norm(), len2 = (ori_p_to - a).norm(), len_left = len1;
  if (len2 < len1)
    len_left = len2;
  double kL = (len_left + pos - rho_v) / rho_e;
  kL = kL < 0 ? 0 : kL;

  if (kL < INT_MAX) kL = (int)kL;

  double kU = kL + 1;
  double pos1 = kL * rho_e + rho_v - len_left, pos2 = kU * rho_e + rho_v - len_left;
  if (pos1 < rho_v) { pos1 = rho_v; pos2 = rho_v + rho_e; }

  if ((LEN(EH(hh)) - pos2) / LEN(EH(hh)) < 1e-5)
    pos = pos1;
  else
  {
    Vec3d d0 = b - a;
    d0.normalize();
    Vec3d p1 = a + d0 * pos1, p2 = a + d0 * pos2;
    Vec3d q = PNT(OPPV(hh));

    if ((p1 - q).length() < (p2 - q).length()) pos = pos1;
    else pos = pos2;
  }
  return pos;
}

double SurfaceDelaunay::find_drop_feet_farthest(HalfedgeHandle cur_hh)const
{
  HalfedgeHandle cur_hh_oppo = OPPH(cur_hh);
  double edge_len = LEN(EH(cur_hh));
  double pos1 = find_pos_nearest_drop_feet(cur_hh);
  double pos2 = find_pos_nearest_drop_feet(cur_hh_oppo);
  double toBound1 = std::min(pos1 / edge_len, 1 - pos1 / edge_len);
  double toBound2 = std::min(pos2 / edge_len, 1 - pos2 / edge_len);
  if (toBound1 < 0 && toBound2 < 0)
  {
    return -1.0;
  }
  if (toBound1 > toBound2)
  {
    return pos1;
  }
  else
  {
    return edge_len - pos2;
  }
}

/// @brief see name
/// @param [in] edge_len [0] hh [1] hh_next [2] hh_prev [3] hh_opp_next [4] hh_opp_prev.
///                        if hh_opp is boundary, only contains the first three elements.
/// @param [in] face_area [0] hh_face [1] hh_opp_face.
///                        if hh_opp is boundary, only contains the first element.
void SurfaceDelaunay::find_split_pos_make_edge_Delaunay(const std::vector<double>& edge_len, const std::vector<double>& face_area, double& lower_bound, double& upper_bound)const
{
  auto& len = edge_len;   // short name
  auto& area = face_area; // short name
  double cosAngle1 = (len[0] * len[0] + len[2] * len[2] - len[1] * len[1]) / (2 * len[0] * len[2]);
  double cosAngle2 = (len[0] * len[0] + len[1] * len[1] - len[2] * len[2]) / (2 * len[0] * len[1]);
  double h1 = 2.0 * area[0] / len[1], h2 = 2.0 * area[0] / len[2];
  if (len.size() == 5)
  {
    double cosAngle3 = (len[0] * len[0] + len[4] * len[4] - len[3] * len[3]) / (2 * len[0] * len[4]);
    double cosAngle4 = (len[0] * len[0] + len[3] * len[3] - len[4] * len[4]) / (2 * len[0] * len[3]);
    double h3 = 2.0 * area[1] / len[3], h4 = 2.0 * area[1] / len[4];
    upper_bound = (len[2] * h3 + len[3] * h2) / (h3 * cosAngle1 + h2 * cosAngle4);
    lower_bound = (len[1] * h4 + len[4] * h1) / (h4 * cosAngle2 + h1 * cosAngle3);
    lower_bound = len[0] - lower_bound;
  }
  else
  {
    upper_bound = len[2] / cosAngle1;
    lower_bound = len[1] / cosAngle2;
    lower_bound = len[0] - lower_bound;
  }
}

double SurfaceDelaunay::find_split_pos_make_surround_edge_Delaunay(HalfedgeHandle hh, HalfedgeHandle sur_hh)const
{
  HalfedgeHandle sur_hh_oppo = OPPH(sur_hh);
  if (mesh->is_boundary(sur_hh_oppo))
  {
    return find_drop_feet(hh);
  }
  double len0 = LEN(EH(sur_hh_oppo)), len1 = LEN(EH(NEXT(sur_hh_oppo))), len2 = LEN(EH(PREV(sur_hh_oppo)));
  double cosAngle1 = (len1 * len1 + len2 * len2 - len0 * len0) / (2.0 * len1 * len2);
  double angle1 = std::acos(cosAngle1);

  HalfedgeHandle next_hh = NEXT(hh), prev_hh = PREV(hh);
  len0 = LEN(EH(hh)), len1 = LEN(EH(next_hh)), len2 = LEN(EH(prev_hh));
  if (sur_hh == next_hh)
  {
    double cosAngle2 = (len0 * len0 + len1 * len1 - len2 * len2) / (2.0 * len0 * len1);
    double angle2 = std::acos(cosAngle2);
    double angle3 = angle1 - angle2;
    return len0 - sin(angle3) * LEN(EH(sur_hh_oppo)) / sin(angle1);
  }
  else
  {
    double cosAngle2 = (len0 * len0 + len2 * len2 - len1 * len1) / (2.0 * len0 * len2);
    double angle2 = std::acos(cosAngle2);
    double angle3 = angle1 - angle2;
    return sin(angle3) * LEN(EH(sur_hh_oppo)) / sin(angle1);
  }
}

/// @brief see name
Vec3d SurfaceDelaunay::find_split_pos_make_edge_Delaunay(EdgeHandle eh)const
{
  HalfedgeHandle hh = mesh->halfedge_handle(eh, 0), hh_oppo = mesh->halfedge_handle(eh, 1);
  if (mesh->is_boundary(hh))
    std::swap(hh, hh_oppo);

  // get adjacent halfedges, length and area.
  std::vector<HalfedgeHandle> all_halfedges; all_halfedges.reserve(5);
  std::vector<double> edge_len; edge_len.reserve(5);
  std::vector<double> face_area; face_area.reserve(2);

  all_halfedges.push_back(hh);
  edge_len.push_back(LEN(EH(all_halfedges.back())));
  all_halfedges.push_back(NEXT(hh));
  edge_len.push_back(LEN(EH(all_halfedges.back())));
  all_halfedges.push_back(PREV(hh));
  edge_len.push_back(LEN(EH(all_halfedges.back())));
  face_area.push_back(AREA(FH(hh)));
  if (!mesh->is_boundary(eh))
  {
    all_halfedges.push_back(NEXT(hh_oppo));
    edge_len.push_back(LEN(EH(all_halfedges.back())));
    all_halfedges.push_back(PREV(hh_oppo));
    edge_len.push_back(LEN(EH(all_halfedges.back())));
    face_area.push_back(AREA(FH(hh_oppo)));
  }

  std::vector<std::pair<double, double>> bounds;
  double lower_bound, upper_bound;
  find_split_pos_make_edge_Delaunay(edge_len, face_area, lower_bound, upper_bound);
  double this_edge_lowerBound = lower_bound, this_edge_upperBound = upper_bound;
  bounds.push_back(std::make_pair(lower_bound, upper_bound));
  upper_bound = find_split_pos_make_surround_edge_Delaunay(all_halfedges[0], all_halfedges[1]);
  lower_bound = find_split_pos_make_surround_edge_Delaunay(all_halfedges[0], all_halfedges[2]);
  bounds.push_back(std::make_pair(0.0, upper_bound));
  bounds.push_back(std::make_pair(lower_bound, edge_len[0]));
  if (all_halfedges.size() == 5)
  {
    lower_bound = find_split_pos_make_surround_edge_Delaunay(OPPH(all_halfedges[0]), all_halfedges[3]);
    upper_bound = find_split_pos_make_surround_edge_Delaunay(OPPH(all_halfedges[0]), all_halfedges[4]);
    bounds.push_back(std::make_pair(0.0, edge_len[0] - upper_bound));
    bounds.push_back(std::make_pair(edge_len[0] - lower_bound, edge_len[0]));
  }
  for (int i = 0; i < bounds.size(); ++i)
  {
    if (bounds[i].first < 0) bounds[i].first = 0;
    if (bounds[i].first > edge_len[0]) bounds[i].first = edge_len[0];
    if (bounds[i].second < 0) bounds[i].second = 0;
    if (bounds[i].second > edge_len[0]) bounds[i].second = edge_len[0];
  }
  std::set<double> end_points;
  for (int i = 1; i < bounds.size(); ++i)
  {
    end_points.insert(bounds[i].first); end_points.insert(bounds[i].second);
  }
  lower_bound = bounds[0].first;
  upper_bound = bounds[0].second;

  int max_Delaunay_num = 0;
  double opt_LB = lower_bound, opt_UB = upper_bound;
  for (double temp : end_points)
  {
    if (temp < bounds[0].first) continue;
    if (temp > bounds[0].second) break;
    double mid = (lower_bound + temp) / 2.0;
    int cur_Delaunay_num = 0;
    for (int i = 1; i < bounds.size(); ++i)
      if (mid > bounds[i].first && mid < bounds[i].second)
        ++cur_Delaunay_num;
    if (cur_Delaunay_num > max_Delaunay_num)
    {
      max_Delaunay_num = cur_Delaunay_num;
      opt_LB = lower_bound; opt_UB = temp;
    }
    lower_bound = temp;
  }
  double mid = (lower_bound + upper_bound) / 2.0;
  int cur_Delaunay_num = 0;
  for (int i = 1; i < bounds.size(); ++i)
    if (mid > bounds[i].first && mid < bounds[i].second)
      ++cur_Delaunay_num;
  if (cur_Delaunay_num > max_Delaunay_num)
  {
    max_Delaunay_num = cur_Delaunay_num;
    opt_LB = lower_bound; opt_UB = upper_bound;
  }

  const Vec3d& a = PNT(FROM(hh)), & b = PNT(TO(hh));
  double pos = find_drop_feet_farthest(all_halfedges[0]);
  Vec3d new_point;
  if (pos < opt_LB || pos > opt_UB)
    new_point = a + (0.5 * (opt_LB + opt_UB)) * ((b - a).normalized());
  else
    new_point = a + pos * ((b - a).normalized());
  if (isnan(new_point[0]) || isnan(new_point[1]) || isnan(new_point[2]))
    new_point = a + (0.5 * (this_edge_lowerBound + this_edge_upperBound)) * (b - a).normalized();

  return new_point;
}

#undef EH
}// namespace Simplification
}// namespace GCLF