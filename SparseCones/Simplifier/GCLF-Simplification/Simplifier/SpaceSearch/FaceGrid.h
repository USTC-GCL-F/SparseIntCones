#pragma once

#include "Grid/StaticUniformGrid.h"
#include "Exact/ExactTriangle.h"
#include "Exact/TriTriIntersect.h"
#include "Mesh/SurfaceMeshDefinition.h"

namespace GCLF
{
using namespace Geometry;
namespace Simplification
{
using namespace SMesh;

namespace Geometry
{

class FaceGridKernel
{
public:
  typedef ExactTriangle<HeavyTriangle> Primitive;
  typedef CalcBoxForExactTriangle CalcBox;
};

class FaceGrid : public StaticUniformGrid<FaceGridKernel>
{
public:
  std::pair<std::pair<Vec3d, double>, int> closest_point(const Vec3d& query);

  bool do_intersect(const ExactTriangle<Triangle>& tri);
public:
  /* interfaces for OpenMesh */
  FaceGrid(SMeshT& mesh);
  FaceGrid(SMeshT& mesh, const std::vector<FaceHandle>& faces);
  void set_from_mesh(SMeshT& mesh);
  void set_from_mesh(SMeshT& mesh, const std::vector<FaceHandle>& faces);
};

}// namespace Geometry
}// namespace Simplification
}// namespace GCLF