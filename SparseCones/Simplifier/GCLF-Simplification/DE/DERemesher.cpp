#include "DERemesher.h"

#define CHECK_SELF_INTER true

namespace GCLF
{
namespace Simplification
{
namespace DE
{

void DERemesher::remesh()
{
  calc_target_length_on_origin();
  find_target_length_on_remeshing();
  double ave_min_angle = calc_ave_min_angle(rm);
  end_angle = ave_min_angle > 45.0 ? 45.0 : std::floor(ave_min_angle);

  size_t outer_loop_iter = 1;
  size_t last_step_vn = rm->n_vertices();
  SMeshT small_step_mesh;
  while (true)
  {
    Logger::user_logger->info("remeshing OUTER loop iter {}.", outer_loop_iter);
    for (size_t inner_loop_iter = 1; inner_loop_iter < 6; ++inner_loop_iter)
    {
      Logger::user_logger->info("remeshing inner loop iter {}.", inner_loop_iter);
      if (ave_min_angle > end_angle)
        small_step_mesh = *rm;
      size_t n_split = split_long_edges();
      size_t n_collapse = collapse_short_edges();
      // end condition 1
      if (n_split == 0 && n_collapse == 0)
        break;
      equalize_valences();
      Laplacian_smooth();
      //std::string name = "output/remesh_" + std::to_string(outer_loop_iter) + "_" + std::to_string(inner_loop_iter) + ".obj";
      //OpenMesh::IO::write_mesh(*rm, name);
      // end condition 2
      ave_min_angle = calc_ave_min_angle(rm);
      if (outer_loop_iter > 1 && rm->n_vertices() <= param->target_vn)
        break;
    }
    rt->rebuild();
    Logger::user_logger->info("average minimal angle {}, end angle {}.", ave_min_angle, end_angle);
    Logger::user_logger->info("vertices number {}, last step vertices number {}", rm->n_vertices(), last_step_vn);
    // end condition
    if (ave_min_angle < end_angle ||
      rm->n_vertices() <= param->target_vn||
      rm->n_vertices() >= last_step_vn)
    {
      // check if a new outer loop is needed
      if (outer_loop_iter == 1)
      {
        if (ave_min_angle < end_angle)
        {
          if ((end_angle - ave_min_angle) < 1.0)
            end_angle -= 1.0;
          else
          {
            small_step_mesh = *rm;
            break;
          }
        }
      }
      else break;
    }
    last_step_vn = rm->n_vertices();
    enlarge_target_length_on_origin();
    outer_loop_iter++;
  }
  *rm = small_step_mesh;

  rm->garbage_collection();
  init_one_ring_faces(rm);
  rt->collect_garbage();
  rt->rebuild();
}

void DERemesher::calc_target_length_on_origin()
{
  // principal curvature
  PrincipalCurvature mesh_curvature(om);
  mesh_curvature.compute_principal_curvature();
  auto& K1 = mesh_curvature.principal_K1;
  auto& K2 = mesh_curvature.principal_K2;
  auto& dir1 = mesh_curvature.principal_dir1;
  auto& dir2 = mesh_curvature.principal_dir2;

  size_t V_N = om->n_vertices();
  std::vector<double> curvature; curvature.reserve(V_N);
  for (size_t i = 0; i < V_N; i++)
    curvature.push_back(std::max(abs(K1[i]), abs(K2[i])));

  // calculate target length on original mesh
  double epsilon = original_diagonal_length * 0.001;
  for (VertexHandle vh : om->vertices())
  {
    double length_temp = sqrt(6 * epsilon / curvature[vh.idx()] - 3 * epsilon * epsilon);
    if (length_temp < epsilon || isnan(length_temp))
      length_temp = epsilon;
    if (length_temp > 40 * epsilon)
      length_temp = 40 * epsilon;
    om->data(vh).target_length = length_temp;
  }
  for (EdgeHandle eh : om->edges())
  {
    HalfedgeHandle hh = om->halfedge_handle(eh, 0);
    VertexHandle v_from = om->from_vertex_handle(hh), v_to = om->to_vertex_handle(hh);
    om->data(eh).target_length = std::min(om->data(v_from).target_length,
      om->data(v_to).target_length);
  }

  // smooth target length on original mesh
  size_t max_smooth_iter = 10;
  for (size_t iter = 0;iter < max_smooth_iter;iter++)
  {
    for (VertexHandle vh : om->vertices())
    {
      double len_range = 0.0;
      for (VertexHandle vv : om->vv_range(vh))
        len_range += om->data(vv).target_length;
      om->data(vh).target_length = len_range / om->valence(vh);
    }
  }
}

void DERemesher::enlarge_target_length_on_origin()
{
  for (VertexHandle vh : om->vertices())
  {
    om->data(vh).target_length *= param->enlargeTargetLengthRatio;
  }
  for (EdgeHandle eh : om->edges())
  {
    HalfedgeHandle hh = om->halfedge_handle(eh, 0);
    VertexHandle v_from = om->from_vertex_handle(hh), v_to = om->to_vertex_handle(hh);
    om->data(eh).target_length = std::min(om->data(v_from).target_length,
      om->data(v_to).target_length);
  }
}

void DERemesher::find_target_length_on_remeshing()
{
  for (VertexHandle vh : rm->vertices())
  {
    rm->data(vh).target_length = om->data(vh).target_length;
  }
  for (EdgeHandle eh : rm->edges())
  {
    HalfedgeHandle hh = rm->halfedge_handle(eh, 0);
    VertexHandle v_from = rm->from_vertex_handle(hh), v_to = rm->to_vertex_handle(hh);
    rm->data(eh).target_length = std::min(rm->data(v_from).target_length,
      rm->data(v_to).target_length);
  }
}

double DERemesher::calc_ave_min_angle(SMeshT* mesh)
{
  double sum_min_angle = 0.0;
  for (OpenMesh::FaceHandle fh : mesh->faces())
  {
    sum_min_angle += calc_small_angle_of_tri(mesh, fh);
  }
  return (sum_min_angle / mesh->n_faces()) * 180.0 / M_PI;
}

double DERemesher::calc_small_angle_of_tri(SMeshT* mesh, OpenMesh::FaceHandle fh)
{
#define HE(p) p->first
#define COS(p) p->second
  auto sectors_angle = face_sector_angle(rm, fh);
  // find min and max angle
  // angle_a > angle_b === cos(angle_a) < cos(angle_b) when angle in [0, pi]
  auto max_angle = sectors_angle.begin(), min_angle = sectors_angle.begin();
  for (auto sa = sectors_angle.begin() + 1;sa != sectors_angle.end();sa++)
  {
    if (COS(sa) < COS(max_angle))max_angle = sa;
    else if (COS(sa) > COS(min_angle))min_angle = sa;
  }
  if (std::isnan(COS(min_angle)))
    return 0.0;
  else
    return std::acos(std::max(std::min(COS(min_angle), 1.0), -1.0));
#undef HE
#undef COS
}

std::pair<Vec3d, double> DERemesher::project_and_target_length(const Vec3d& p)
{
  auto cpp = og->closest_point(p); // closest point and primitive
  const Vec3d& cp = cpp.first.first;     // closest point
  OpenMesh::FaceHandle cf(cpp.second);   // closest face handle
  OpenMesh::HalfedgeHandle hh = om->halfedge_handle(cf);
  OpenMesh::VertexHandle v_from = om->from_vertex_handle(hh), v_to = om->to_vertex_handle(hh), v_oppo = om->opposite_vh(hh);
  double new_tar = 0.0;
  double face_area = om->data(cf).face_area;
  if (face_area < GeomThr::area_de_thr)
    new_tar = (om->data(v_oppo).target_length + om->data(v_from).target_length
      + om->data(v_to).target_length) / 3.0;
  else
  {
    const Vec3d& p_from = om->point(v_from), & p_to = om->point(v_to), & p_oppo = om->point(v_oppo);
    double oppo_area = triangle_area(cp, p_from, p_to);
    double oppo_weight = oppo_area / face_area;
    double from_area = triangle_area(cp, p_to, p_oppo);
    double from_weight = from_area / face_area;
    double to_area = triangle_area(cp, p_oppo, p_from);
    double to_weight = to_area / face_area;
    new_tar = om->data(v_oppo).target_length * oppo_weight + om->data(v_from).target_length * from_weight
      + om->data(v_to).target_length * to_weight;
  }
  return std::make_pair(cp, new_tar);
}

/// @brief collapse edges shorter than target length.
size_t DERemesher::collapse_short_edges()
{
  auto edge_collapser = new_edge_collapser();
  edge_collapser.set_flags(true, true, true, true, true,
    true, CHECK_SELF_INTER, /*check_inter*/false);

  size_t collapsed_num = 0;
  size_t m_num = rm->n_vertices();
  for (size_t iter = 0; iter < param->collapseIter;iter++)
  {
    for (OpenMesh::EdgeHandle eh : rm->edges())
    {
      if (rm->status(eh).deleted())
        continue;

      if (edge_collapser.try_collapse_avoid_long_short_edge_avoid_out_of_bound(eh, max_distance_error_value))
      {
          collapsed_num++;
          if ((m_num - collapsed_num) == 5000)
          {
              goto flag1;
          }
      }
    }
  }
  flag1:

  Logger::user_logger->info("collapsed [{}] edges.", collapsed_num);
  rm->garbage_collection();
  rt->collect_garbage();
  init_one_ring_faces(rm);
  Logger::user_logger->info("[{}] vertices and [{}] faces remained.", rm->n_vertices(), rm->n_faces());
  return collapsed_num;
}

size_t DERemesher::Laplacian_smooth()
{
  auto vertex_relocater = new_vertex_relocater();
  vertex_relocater.set_flags(
    true, true, true, true,
    true, false, CHECK_SELF_INTER);

  size_t smooth_num = 0;
  for (size_t it = 0;it < param->smoothIter;it++)
  {
    size_t smooth_num_each_iter = 0;
    for (VertexHandle v : rm->vertices())
    {
      if (rm->status(v).deleted())
        continue;
      vertex_relocater.init(v);
      Vec3d new_point = vertex_relocater.find_smooth_target_with_target_length();
      auto pp_nt = project_and_target_length(new_point);
      const Vec3d& project_point = pp_nt.first;
      if ((project_point - rm->point(v)).norm() < GeomThr::relocate_length_thr)
        continue; // too short relocate length

      vertex_relocater.set_new_target_length(pp_nt.second);
      if (vertex_relocater.try_reloate_vertex_avoid_out_of_bound(project_point, max_distance_error_value))
      {
        smooth_num_each_iter++;
      }
    }
    if (smooth_num_each_iter == 0)
      break;
    smooth_num += smooth_num_each_iter;
  }
  Logger::user_logger->info("smoothed [{}] vertices.", smooth_num);
  return smooth_num;
}

size_t DERemesher::equalize_valences()
{
  auto edge_flipper = new_edge_flipper();
  edge_flipper.set_flags(
    true, true, true, true, true,
    true, false, CHECK_SELF_INTER);

  size_t flipped_num = 0;
  for (size_t it = 0;it < param->equalizeValenceIter;it++)
  {
    for (EdgeHandle e : rm->edges())
    {
      if (edge_flipper.try_flip_edge_decrease_valence_avoid_out_of_bound(e, max_distance_error_value))
      {
        flipped_num++;
      }
    }
  }

  Logger::user_logger->info("flipped [{}] edges.", flipped_num);
  return flipped_num;
}

size_t DERemesher::split_long_edges()
{
  auto edge_splitter = new_edge_splitter();
  edge_splitter.set_flags(true, true, true, true, true);

  size_t split_num = 0;
  for (size_t it = 0;it < param->splitIter;it++)
  {
    for (EdgeHandle e : rm->edges())
    {
      double high_edge_len = 4.0 / 3.0 * rm->data(e).target_length;
      if (rm->data(e).edge_length >= high_edge_len)
      {
        // split edge
        edge_splitter.split(e);
        split_num++;
      }
    }
    // end condition
    if (split_num == 0)
      break;
  }

  Logger::user_logger->info("splitted [{}] edges.", split_num);
  init_one_ring_faces(rm);
  return split_num;
}
}// namespace DE
}// namespace Simplification
}// namespace GCLF

#undef CHECK_SELF_INTER