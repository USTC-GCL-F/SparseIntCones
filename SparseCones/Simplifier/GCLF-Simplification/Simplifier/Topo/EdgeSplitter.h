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

class EdgeSplitter
{
public:
  EdgeSplitter(
    SMeshT* _origin, SMeshT* _remeshing,
    DFaceTree* _origin_tree, DFaceTree* _remeshing_tree)
    :om(_origin), rm(_remeshing),
    ot(_origin_tree), rt(_remeshing_tree),
    f_update_links(false),
    f_update_target_length(false),
    f_update_normals(false)
  {}
public:
  void set_flags(
    bool _f_update_links, bool _f_update_target_length, bool _f_update_normals,
    bool _f_update_tree, bool _f_update_one_ring_faces
  );
  VertexHandle split(EdgeHandle e);
  VertexHandle split(EdgeHandle e, const Vec3d& split_point);
public:
  // input
  SMeshT* om;
  SMeshT* rm;

  DFaceTree* ot;
  DFaceTree* rt;
  // flags
  bool f_update_links;
  bool f_update_target_length;
  bool f_update_normals;
  bool f_update_tree;
  bool f_update_one_ring_faces;
  // handles after split, notations come from OpenMesh
  VertexHandle split_v;
  Vec3d split_point;
  HalfedgeHandle h0, o0, t0, t2;
  FaceHandle f0, f1, f2, f3;
  double new_target_length;

  void set_handles_after_split();

  void update_target_len();
  void update_tree();
  void update_length_and_area();
  void update_one_ring_faces();
  void update_normals();

  void split_in_links();
};
}// namespace Simplification
}// namespace GCLF