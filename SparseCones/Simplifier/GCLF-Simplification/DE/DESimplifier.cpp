#include "DESimplifier.h"

namespace GCLF
{
namespace Simplification
{
namespace DE
{

void DESimplifier::calc_max_distance_error_value()
{
  if (param->errorRatio <= 0.1)
    max_distance_error_value = original_diagonal_length * (param->errorRatio * 0.90) * 0.01;
  else
    max_distance_error_value = original_diagonal_length * (param->errorRatio * 0.85) * 0.01;
  double real_error_value = original_diagonal_length * param->errorRatio * 0.01;
  p_epsilon = { real_error_value,real_error_value,real_error_value };
}

void DESimplifier::initialize(const SMeshT& original)
{
  om = std::make_unique<SMeshT>();
  *om = original;
  rm = std::make_unique<SMeshT>();
  *rm = *om;

  pre_calculate_edge_length(om.get());
  pre_calculate_edge_length(rm.get());
  pre_calculate_face_area(om.get());
  pre_calculate_face_area(rm.get());

  if (!om->has_face_normals())
    om->request_face_normals();
  if (!om->has_vertex_normals())
    om->request_vertex_normals();
  om->update_normals();

  if (!rm->has_face_normals())
    rm->request_face_normals();
  if (!rm->has_vertex_normals())
    rm->request_vertex_normals();
  rm->update_normals();

  init_one_ring_faces(rm.get());

  vt = std::make_unique<VertexTree>(*om);
  ot = std::make_unique<DFaceTree>(*om);
  rt = std::make_unique<DFaceTree>(*rm);
  og = std::make_unique<FaceGrid>(*om);

  // initialize hausdorff distance.
  generate_out_links(om.get(), rm.get(), rt.get());
  calc_face_in_error(rm.get(), std::vector<FaceHandle>(rm->faces_begin(), rm->faces_end()));
#ifdef USE_TREE_SEARCH
  calc_out_error(rm.get(), om.get(), ot.get(), infinite_fp);
#else
  calc_out_error(rm.get(), om.get(), og.get(), infinite_fp);
#endif

  BoundingBox bbox;
  for (const auto& vh : om->vertices())
  {
    const auto& point = om->point(vh);
    bbox += point;
  }
  original_diagonal_length = (bbox.max() - bbox.min()).norm();
  calc_max_distance_error_value();

  remesher = std::make_unique<DERemesher>(
    om.get(), rm.get(), vt.get(), ot.get(), rt.get(), og.get(), &param->paramRemesher_DE);
  remesher->original_diagonal_length = original_diagonal_length;
  remesher->max_distance_error_value = max_distance_error_value;

  delaunay_mesher = std::make_unique<DEDM>(
    om.get(), rm.get(), vt.get(), ot.get(), rt.get(), og.get(), &param->paramDM_DE);
  delaunay_mesher->original_diagonal_length = original_diagonal_length;
  delaunay_mesher->max_distance_error_value = max_distance_error_value;
  delaunay_mesher->p_epsilon = p_epsilon;
}

void DESimplifier::simplify()
{
  remesher->remesh();
  //delaunay_mesher->run();
}
}// namespace DE
}// namespace Simplification
}// namespace GCLF