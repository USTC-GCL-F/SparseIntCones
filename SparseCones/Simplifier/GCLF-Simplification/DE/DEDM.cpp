#include "DEDM.h"
#include <random>

#define CHECK_SELF_INTER true

namespace GCLF
{
namespace Simplification
{
namespace DE
{

void DEDM::run()
{
  surf_delaunay.initialize(rm);
  surf_delaunay.generate();

  update_rm();

  initialize_simplification();
  simplification();

  surf_delaunay.check_delaunay();
}

void DEDM::update_rm()
{
  *rt = DFaceTree(*rm);
  if (!rm->has_face_normals())
    rm->request_face_normals();
  if (!rm->has_vertex_normals())
    rm->request_vertex_normals();
  rm->update_normals();
  pre_calculate_edge_length(rm);
  pre_calculate_face_area(rm);
  init_one_ring_faces(rm);
  // initialize hausdorff distance.
  clear_links(rm);
  generate_out_links(om, rm, rt);
  calc_face_in_error(rm, std::vector<FaceHandle>(rm->faces_begin(), rm->faces_end()));
#ifdef USE_TREE_SEARCH
  calc_out_error(rm, om, ot, infinite_fp);
#else
  calc_out_error(rm, om, og, infinite_fp);
#endif
}

BoundingBox DEDM::calc_local_epsilon_bbox(SMeshT* mesh, const std::vector<FaceHandle>& faces)
{
  BoundingBox bbox;
  // check if links exist on these faces.
  bool have_links = false;
  for (FaceHandle fh : faces)
  {
    if (!mesh->data(fh).face_in_links.empty())
    {
      have_links = true;
      break;
    }
  }

  if (have_links)
  {
    // we contruct box from links
    for (FaceHandle fh : faces)
      for (const Link& link : mesh->data(fh).face_in_links)
        bbox += link.first;
  }
  else
  {
    // we construct box from faces
    for (FaceHandle fh : faces)
      for (VertexHandle vh : mesh->fv_range(fh))
        bbox += mesh->point(vh);
  }
  bbox.enlarge(p_epsilon);
  return bbox;
}

std::vector<Vec3d> DEDM::generate_samples_by_links(SMeshT* mesh, const std::vector<FaceHandle>& one_ring_faces, size_t sample_num)
{
  std::vector<Vec3d> samples; samples.reserve(sample_num);
  size_t faces_num = one_ring_faces.size();
  size_t samples_per_face = sample_num / faces_num;
  for (FaceHandle fh : one_ring_faces)
  {
    size_t sample_cnt = 0;
    for (const Link& link : mesh->data(fh).face_in_links)
    {
      if (sample_cnt >= samples_per_face)
        break;

      samples.push_back(link.first);
      sample_cnt++;
    }
  }
  if (samples.size() < sample_num)
  {
    size_t remain_samples_num = sample_num - samples.size();
    std::vector<Vec3d> remain_samples = generate_random_samples(mesh, one_ring_faces[0], remain_samples_num);
    samples.insert(samples.end(), remain_samples.begin(), remain_samples.end());
  }
  return samples;
}

std::vector<Vec3d> DEDM::generate_random_samples(SMeshT* mesh, FaceHandle fh, size_t sample_num)
{
  std::vector<Vec3d> samples; samples.reserve(sample_num);

  HalfedgeHandle he = mesh->halfedge_handle(fh);
  Vec3d a = mesh->point(mesh->to_vertex_handle(he));
  Vec3d b = mesh->point(mesh->opposite_vh(he));
  Vec3d c = mesh->point(mesh->from_vertex_handle(he));

  if (sample_num == 1)
  {
    samples.push_back(mesh->calc_face_centroid(fh));
  }
  else if (sample_num == 2)
  {
    double bc = (b - c).norm();
    double ca = (c - a).norm();
    double ab = (a - b).norm();
    Vec3d d;
    int longest_edge_index;
    if (bc > ca)
    {
      longest_edge_index = bc > ab ? 0 : 2;
    }
    else
    {
      longest_edge_index = ab < ca ? 1 : 2;
    }
    switch (longest_edge_index)
    {
    case 0:
      d = (b + c) / 2.0;
      samples.push_back((a + b + d) / 3.0);
      samples.push_back((a + d + c) / 3.0);
      break;
    case 1:
      d = (a + c) / 2.0;
      samples.push_back((b + c + d) / 3.0);
      samples.push_back((b + d + a) / 3.0);
      break;
    case 2:
      d = (a + b) / 2.0;
      samples.push_back((c + a + d) / 3.0);
      samples.push_back((c + d + b) / 3.0);
      break;
    default:
      break;
    }
  }
  else
  {
    Vec3d ab = b - a, ac = c - a;
    while (samples.size() < sample_num)
    {
      double  u = 0.0, v = 0.0;
      u = (double)rand() / (double)RAND_MAX;
      v = (double)rand() / (double)RAND_MAX;
      if (u + v > 1.0)
      {
        u = 1.0 - u, v = 1.0 - v;
      }
      samples.push_back(a + u * ab + v * ac);
    }
  }
  return samples;
}

void DEDM::relocate_Hausdorff()
{
  generate_out_links(om, rm, rt);
  calc_face_in_error(rm, std::vector<FaceHandle>(rm->faces_begin(), rm->faces_end()));
#ifdef USE_TREE_SEARCH
  calc_out_error(rm, om, ot, infinite_fp);
#else
  calc_out_error(rm, om, og, infinite_fp);
#endif

  auto vertex_relocater = new_vertex_relocater();
  vertex_relocater.set_flags(
    true, true, true, true,
    true, false, CHECK_SELF_INTER);

  bool try_relocating_any_vertex = false;

  for (size_t iter = 0;iter < 3;iter++)
  {
    try_relocating_any_vertex = rm->n_vertices() <= 3000;
    try_relocating_any_vertex = true;
    for (VertexHandle vh : rm->vertices())
    {
      Logger::dev_logger->trace("try relocate vertex {}", vh.idx());
      // check if we need to relocate this vertex
      double max_local_error = rm->data(vh).out_surround_error;
      if (!try_relocating_any_vertex && max_local_error < max_distance_error_value)
        continue;

      Vec3d new_point;
      std::map<Vec3d, LinkVector> new_local_in_links;
      if (predict_relocate_pos(vh, vertex_relocater, max_local_error, new_point, new_local_in_links))
      {
        // relocate and update links
        vertex_relocater.set_local_mesh_in_links(std::move(new_local_in_links));
        vertex_relocater.relocate(new_point);

        // update error
        std::vector<FaceHandle> one_ring_faces = std::vector<FaceHandle>(
          vertex_relocater.one_ring_faces.begin(),
          vertex_relocater.one_ring_faces.end());
        calc_face_in_error(rm, one_ring_faces);
        calc_vertex_out_error(rm, om, og, vh);
        calc_out_error(rm, om, og, one_ring_faces, infinite_fp);
        calc_out_surround_error(rm, vh);
        for (VertexHandle vv : rm->vv_range(vh))
          calc_out_surround_error(rm, vv);
      }
    }
  }
}


bool DEDM::predict_relocate_pos(VertexHandle vh, VertexRelocater& vertex_relocater, double max_error,
  Vec3d& point, std::map<Vec3d, LinkVector>& local_in_links)
{
  if (!vertex_relocater.init(vh))
    return false;

  VertexHandle local_center_vh;
  std::shared_ptr<SMeshT> local_mesh = construct_local_mesh(rm, vertex_relocater.get_halfedges(),
    rm->point(vh), local_center_vh, !rm->is_boundary(vh));

  std::vector<FaceHandle> one_ring_faces = std::vector<FaceHandle>(
    vertex_relocater.one_ring_faces.begin(),
    vertex_relocater.one_ring_faces.end());

  BoundingBox bbox = calc_local_epsilon_bbox(rm, one_ring_faces);

  // random generators
  std::default_random_engine generator1;
  std::uniform_real_distribution<double> dis1(0.0, std::nextafter(1.0, DBL_MAX));
  std::default_random_engine generator2;
  std::uniform_int_distribution<int> dis2(0, (int)param->Np - 1);
  std::default_random_engine generator3;
  std::uniform_int_distribution<int> dis3(0, 2);
  int rand1, rand2, rand3;

  // populations
  std::vector<Vec3d> old_population, mutation_population, crossover_population, new_population;

  // Hausdorff distances and links
  std::vector<double> Hausdorff_distances(param->Np, DBL_MAX);
  std::vector<std::map<Vec3d, LinkVector>> new_local_in_links(param->Np);
  size_t new_min_HD_idx = 0;      // HD: Hausdorff distance
  double new_min_HD = DBL_MAX, old_min_HD = DBL_MAX;

  // local meshes for multi-threads
  std::vector<SMeshT> local_meshes; local_meshes.resize(omp_get_max_threads(), *local_mesh);

  old_population = generate_samples_by_links(rm, one_ring_faces, param->Np);
  new_population = old_population;

  /****** Begin DE *******/

  // initialize
#pragma omp parallel for
  for (int j = 0; j < param->Np; ++j)
  {
    if (surf_delaunay.satisfy_Delaunay(vh, old_population[j]))
    {
      int current_thread_num = omp_get_thread_num();
      ASSERT(current_thread_num < local_meshes.size(), "threads size error.");

      SMeshT* lm = &local_meshes[current_thread_num]; // local mesh
      clear_links(lm);
      lm->point(local_center_vh) = old_population[j];
      pre_calculate_edge_length(lm);
      pre_calculate_face_area(lm);
      Hausdorff_distances[j] = vertex_relocater.local_Hausdorff_after_relocating(
        lm, old_population[j], infinite_fp, new_local_in_links[j]
      );
    }
  }

  new_min_HD_idx = std::distance(Hausdorff_distances.begin(), std::min_element(Hausdorff_distances.begin(), Hausdorff_distances.end()));
  new_min_HD = Hausdorff_distances[new_min_HD_idx];
  old_min_HD = new_min_HD;

  //Logger::dev_logger->trace("Initialize relocate_DE: minimal HD is {}.", new_min_HD);

  // enter loop
  mutation_population.resize(param->Np);
  size_t iter = 0, consecutive_iter = 0;
  while (true)
  {
    // generate mutation population
  #pragma omp parallel for
    for (int j = 0; j < param->Np; ++j)
    {
      rand1 = dis2(generator2);
      while (rand1 == j)
        rand1 = dis2(generator2);
      rand2 = dis2(generator2);
      while (rand2 == j || rand2 == rand1)
        rand2 = dis2(generator2);
      rand3 = dis2(generator2);
      while (rand3 == j || rand3 == rand1 || rand3 == rand2)
        rand3 = dis2(generator2);
      mutation_population[j] = old_population[rand1] + param->mutationScale * (old_population[rand2] - old_population[rand3]);
      mutation_population[j].maximize(bbox.min());
      mutation_population[j].minimize(bbox.max());
    }
    // generate crossover population
    crossover_population = old_population;
  #pragma omp parallel for
    for (int j = 0; j < param->Np; ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        if (dis1(generator1) <= param->crossoverRate || i == dis3(generator3))
          crossover_population[j][i] = mutation_population[j][i];
      }
    }
    // generate new population
  #pragma omp parallel for
    for (int j = 0; j < param->Np; ++j)
    {
      if (surf_delaunay.satisfy_Delaunay(vh, crossover_population[j]))
      {
        int current_thread_num = omp_get_thread_num();
        ASSERT(current_thread_num < local_meshes.size(), "threads size error.");

        SMeshT* lm = &local_meshes[current_thread_num]; // local mesh
        clear_links(lm);
        lm->point(local_center_vh) = crossover_population[j];
        pre_calculate_edge_length(lm);
        pre_calculate_face_area(lm);
        std::map<Vec3d, LinkVector> tmp_local_in_links;

        double tmp_local_Hausdorff = vertex_relocater.local_Hausdorff_after_relocating(
          lm, crossover_population[j], infinite_fp, tmp_local_in_links);
        if (tmp_local_Hausdorff <= Hausdorff_distances[j])
        {
          new_population[j] = crossover_population[j];
          Hausdorff_distances[j] = tmp_local_Hausdorff;
          new_local_in_links[j] = std::move(tmp_local_in_links);
        }
      }
    }
    new_min_HD_idx = distance(Hausdorff_distances.begin(), min_element(Hausdorff_distances.begin(), Hausdorff_distances.end()));
    new_min_HD = Hausdorff_distances[new_min_HD_idx];
    //Logger::dev_logger->trace("Loop iter {}: minimal HD is {}.", iter, new_min_HD);
    // end condition one: check if loop is convergent
    if (old_min_HD != DBL_MAX)
    {
      double relative_change = (old_min_HD - new_min_HD) / old_min_HD;
      if (relative_change <= convergence_rate)
        consecutive_iter++;
      else
        consecutive_iter = 0;
      if (consecutive_iter == max_consecutive_iter)
      {
        //Logger::dev_logger->trace("loop is convergent.");
        //Logger::dev_logger->trace("end minimal HD is {}. relocate? {}", new_min_HD, new_min_HD <= max_error);
        if (new_min_HD <= max_error)
        {
          point = new_population[new_min_HD_idx];
          local_in_links = std::move(new_local_in_links[new_min_HD_idx]);
          return true;
        }
        else
          return false;
      }
    }
    old_min_HD = new_min_HD;
    old_population = new_population;
    ++iter;
    // end condition two: check if loop reaches the maximum iteration number.
    if (iter == DE_iter_num)
    {
      //Logger::dev_logger->trace("maximum iteration number.");
      //Logger::dev_logger->trace("end minimal HD is {}. relocate? {}", new_min_HD, new_min_HD <= max_error);
      if (new_min_HD <= max_error)
      {
        point = new_population[new_min_HD_idx];
        local_in_links = std::move(new_local_in_links[new_min_HD_idx]);
        return true;
      }
      else
        return false;
    }
  }
}

void DEDM::initialize_simplification()
{
  queue_cost = decltype(queue_cost)();
  queue_state.clear(); queue_state.resize(rm->n_edges(), 0);
  queue_new_points.resize(rm->n_edges());
  queue_new_local_in_links.resize(rm->n_edges());

  Logger::user_logger->info("begin intiallizing simplification.");
  for (EdgeHandle eh : rm->edges())
  {
    if (!rm->status(eh).deleted())
      update_cost(eh);
  }
  Logger::user_logger->info("end intiallizing simplification.");
}

void DEDM::update_cost(EdgeHandle eh)
{
  auto edge_collapser = new_edge_collapser();
  edge_collapser.set_flags(true, false, true, true, true,
    true, true, false);

  queue_state[eh.idx()]++;

  Vec3d new_point;
  double cost;
  std::map<Vec3d, LinkVector> new_local_in_links;
  if (constrained_DM_simplification(eh, edge_collapser, max_distance_error_value, cost, new_point, new_local_in_links))
  {
    queue_cost.emplace(cost, eh.idx(), queue_state[eh.idx()]);
    queue_new_points[eh.idx()] = new_point;
    queue_new_local_in_links[eh.idx()] = std::move(new_local_in_links);
  }
}

void DEDM::update_cost(VertexHandle vh)
{
  std::vector<VertexHandle> surrounding_vertices;
  surrounding_vertices.reserve(2 * rm->valence(vh));
  for (HalfedgeHandle voh : rm->voh_range(vh))
  {
    surrounding_vertices.push_back(rm->to_vertex_handle(voh));
    if (!rm->is_boundary(rm->opposite_halfedge_handle(rm->next_halfedge_handle(voh))))
      surrounding_vertices.push_back(rm->opposite_he_opposite_vh(rm->next_halfedge_handle(voh)));
  }

  std::set<EdgeHandle> surrounding_edges;
  for (VertexHandle temp_vh : surrounding_vertices)
    for (EdgeHandle ve : rm->ve_range(temp_vh))
      surrounding_edges.insert(ve);

  for (EdgeHandle eh : surrounding_edges)
  {
    update_cost(eh);
  }
}

bool DEDM::constrained_DM_simplification(
  EdgeHandle eh, EdgeCollapser& edge_collapser, double max_error,
  double& cost, Vec3d& point, std::map<Vec3d, LinkVector>& local_in_links)
{
  if (!edge_collapser.init(eh))
    return false;

  HalfedgeHandle hh = rm->halfedge_handle(eh, 0);
  HalfedgeHandle hh_opp = rm->halfedge_handle(eh, 1);
  if (rm->is_boundary(hh))
    std::swap(hh, hh_opp);
  VertexHandle from_v = rm->from_vertex_handle(hh);
  VertexHandle to_v = rm->to_vertex_handle(hh);

  const std::vector<HalfedgeHandle>& halfedges = edge_collapser.halfedges;

  bool is_disk = !rm->is_boundary(eh) && !rm->is_boundary(from_v) && !rm->is_boundary(to_v);
  VertexHandle local_center_vh;
  std::shared_ptr<SMeshT> local_mesh = construct_local_mesh(rm, halfedges,
    rm->point(to_v), local_center_vh, is_disk);

  std::vector<FaceHandle> one_ring_faces = std::vector<FaceHandle>(
    edge_collapser.one_ring_faces.begin(),
    edge_collapser.one_ring_faces.end());
  BoundingBox bbox = calc_local_epsilon_bbox(rm, one_ring_faces);

  // random generators
  std::default_random_engine generator1;
  std::uniform_real_distribution<double> dis1(0.0, std::nextafter(1.0, DBL_MAX));
  std::default_random_engine generator2;
  std::uniform_int_distribution<int> dis2(0, (int)param->Np - 1);
  std::default_random_engine generator3;
  std::uniform_int_distribution<int> dis3(0, 2);
  int rand1, rand2, rand3;

  // populations
  std::vector<Vec3d> old_population, mutation_population, crossover_population, new_population;

  // Hausdorff distances and links
  std::vector<double> Hausdorff_distances(param->Np, DBL_MAX);
  std::vector<std::map<Vec3d, LinkVector>> new_local_in_links(param->Np);
  size_t new_min_HD_idx = 0;      // HD: Hausdorff distance
  double new_min_HD = DBL_MAX, old_min_HD = DBL_MAX;

  // local meshes for multi-threads
  std::vector<SMeshT> local_meshes; local_meshes.resize(omp_get_max_threads(), *local_mesh);

  old_population = generate_samples_by_links(rm, one_ring_faces, param->Np - 3);
  old_population.push_back(rm->point(from_v));
  old_population.push_back(rm->point(to_v));
  old_population.push_back(rm->calc_edge_midpoint(eh));
  new_population = old_population;

  /****** Begin DE *******/

  // initialize
#pragma omp parallel for
  for (int j = 0; j < param->Np; ++j)
  {
    if (surf_delaunay.satisfy_Delaunay(halfedges, old_population[j]))
    {
      int current_thread_num = omp_get_thread_num();
      ASSERT(current_thread_num < local_meshes.size(), "threads size error.");

      SMeshT* lm = &local_meshes[current_thread_num]; // local mesh
      clear_links(lm);
      lm->point(local_center_vh) = old_population[j];
      pre_calculate_edge_length(lm);
      pre_calculate_face_area(lm);
      Hausdorff_distances[j] = edge_collapser.local_Hausdorff_after_collapsing(
        lm, old_population[j], infinite_fp, new_local_in_links[j]
      );
    }
  }

  new_min_HD_idx = std::distance(Hausdorff_distances.begin(), std::min_element(Hausdorff_distances.begin(), Hausdorff_distances.end()));
  new_min_HD = Hausdorff_distances[new_min_HD_idx];
  old_min_HD = new_min_HD;

  if (new_min_HD <= max_error)
  {
    point = old_population[new_min_HD_idx];
    cost = new_min_HD;
    local_in_links = new_local_in_links[new_min_HD_idx];
    return true;
  }
  old_min_HD = new_min_HD;

  //Logger::dev_logger->trace("Initialize collapse_DE: minimal HD is {}.", new_min_HD);

  // enter loop
  mutation_population.resize(param->Np);
  size_t iter = 0, consecutive_iter = 0;
  while (true)
  {
    // generate mutation population
  #pragma omp parallel for
    for (int j = 0; j < param->Np; ++j)
    {
      rand1 = dis2(generator2);
      while (rand1 == j)
        rand1 = dis2(generator2);
      rand2 = dis2(generator2);
      while (rand2 == j || rand2 == rand1)
        rand2 = dis2(generator2);
      rand3 = dis2(generator2);
      while (rand3 == j || rand3 == rand1 || rand3 == rand2)
        rand3 = dis2(generator2);
      mutation_population[j] = old_population[rand1] + param->mutationScale * (old_population[rand2] - old_population[rand3]);
      mutation_population[j].maximize(bbox.min());
      mutation_population[j].minimize(bbox.max());
    }
    // generate crossover population
    crossover_population = old_population;
  #pragma omp parallel for
    for (int j = 0; j < param->Np; ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        if (dis1(generator1) <= param->crossoverRate || i == dis3(generator3))
          crossover_population[j][i] = mutation_population[j][i];
      }
    }
    // generate new population
  #pragma omp parallel for
    for (int j = 0; j < param->Np; ++j)
    {
      if (surf_delaunay.satisfy_Delaunay(halfedges, crossover_population[j]))
      {
        int current_thread_num = omp_get_thread_num();
        ASSERT(current_thread_num < local_meshes.size(), "threads size error.");

        SMeshT* lm = &local_meshes[current_thread_num]; // local mesh
        clear_links(lm);
        lm->point(local_center_vh) = crossover_population[j];
        pre_calculate_edge_length(lm);
        pre_calculate_face_area(lm);
        std::map<Vec3d, LinkVector> tmp_local_in_links;

        double tmp_local_Hausdorff = edge_collapser.local_Hausdorff_after_collapsing(
          lm, crossover_population[j], infinite_fp, tmp_local_in_links);
        if (tmp_local_Hausdorff <= Hausdorff_distances[j])
        {
          new_population[j] = crossover_population[j];
          Hausdorff_distances[j] = tmp_local_Hausdorff;
          new_local_in_links[j] = std::move(tmp_local_in_links);
        }
      }
    }
    new_min_HD_idx = distance(Hausdorff_distances.begin(), min_element(Hausdorff_distances.begin(), Hausdorff_distances.end()));
    new_min_HD = Hausdorff_distances[new_min_HD_idx];
    //Logger::dev_logger->trace("Loop iter {}: minimal HD is {}.", iter, new_min_HD);
    // end condition one: check if loop is convergent
    if (old_min_HD != DBL_MAX)
    {
      double relative_change = (old_min_HD - new_min_HD) / old_min_HD;
      if (relative_change <= convergence_rate)
        consecutive_iter++;
      else
        consecutive_iter = 0;
      if (consecutive_iter == max_consecutive_iter)
      {
        //Logger::dev_logger->trace("loop is convergent.");
        //Logger::dev_logger->trace("end minimal HD is {}. collapse? {}", new_min_HD, new_min_HD <= max_error);
        if (new_min_HD <= max_error)
        {
          point = new_population[new_min_HD_idx];
          local_in_links = std::move(new_local_in_links[new_min_HD_idx]);
          return true;
        }
        else
          return false;
      }
    }
    old_min_HD = new_min_HD;
    old_population = new_population;
    ++iter;
    // end condition two: check if loop reaches the maximum iteration number.
    if (iter == DE_iter_num)
    {
      //Logger::dev_logger->trace("maximum iteration number.");
      //Logger::dev_logger->trace("end minimal HD is {}. relocate? {}", new_min_HD, new_min_HD <= max_error);
      if (new_min_HD <= max_error)
      {
        point = new_population[new_min_HD_idx];
        local_in_links = std::move(new_local_in_links[new_min_HD_idx]);
        return true;
      }
      else
        return false;
    }
  }
}

void DEDM::simplification()
{
  auto edge_collapser = new_edge_collapser();
  edge_collapser.set_flags(true, false, true, true, true,
    true, true, false);

  size_t total_vertices = rm->n_vertices();
  size_t collapsed_edge = 0;
  while (!queue_cost.empty())
  {
    EdgePriority temp_cost = queue_cost.top();
    queue_cost.pop();

    int e_id = temp_cost.e_id;
    EdgeHandle eh = rm->edge_handle(e_id);
    if (rm->status(eh).deleted() || temp_cost.state != queue_state[e_id])
      continue;

    if (edge_collapser.init(eh))
    {
      collapsed_edge++;
      if (collapsed_edge % 10 == 0)
        Logger::user_logger->info("collapse {} edges, remain {} vertices", collapsed_edge, total_vertices - collapsed_edge);
      if (total_vertices - collapsed_edge <= param->target_vn)
        break;
      edge_collapser.set_local_mesh_in_links(std::move(queue_new_local_in_links[eh.idx()]));
      edge_collapser.collapse_edge(edge_collapser.collapse_he, queue_new_points[eh.idx()], nullptr);
      update_cost(edge_collapser.center_vh);
    }
  }
  rm->garbage_collection();
  rm->update_normals();
  init_one_ring_faces(rm);
  rt->collect_garbage();
  rt->rebuild();
}

}// namespace DE
}// namespace Simplification
}// namespace GCLF

#undef CHECK_SELF_INTER