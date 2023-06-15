#pragma once

#include "SurfaceMeshDefinition.h"

using namespace OpenMesh;

namespace GCLF
{
namespace Simplification
{
namespace SMesh
{

class MeshConnectRegion
{
public:
  MeshConnectRegion() = default;
  MeshConnectRegion(SMeshT* _mesh) :mesh(_mesh) {}

  void find();
  inline const std::vector<std::vector<FaceHandle>>& get() { return regions; }
private:
  SMeshT* mesh;
  std::vector<std::vector<FaceHandle>> regions;
};

}// namespace SMesh
}// namespace Simplification
}// namespace GCLF