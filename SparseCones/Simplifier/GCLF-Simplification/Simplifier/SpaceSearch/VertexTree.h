#pragma once

#include "AABB/AABBTraits.h"
#include "Mesh/SurfaceMeshDefinition.h"

namespace GCLF
{
using namespace Geometry;
namespace Simplification
{
using namespace SMesh;
namespace Geometry
{

class VertexTree : public KdTree
{
public:
  VertexTree(SMeshT& mesh);

  SMeshT::VertexHandle closest_vertex(const Vec3d& query);
};
}// namespace Geometry
}// namespace Simplification
}// namespace GCLF