#pragma once

#include "DERemesher.h"
#include "DEDM.h"

namespace GCLF
{
namespace Simplification
{
using namespace Utils;
using namespace Geometry;
using namespace SMesh;
namespace DE
{

class DESimplifier
{
public:
  // input
  std::unique_ptr<SMeshT> om;
  std::unique_ptr<SMeshT> rm;
  // tree
  std::unique_ptr<VertexTree> vt;
  std::unique_ptr<DFaceTree> ot;
  std::unique_ptr<DFaceTree> rt;
  std::unique_ptr<FaceGrid> og;
  // parameters
  ParamDESimplifier* param;

  DESimplifier(ParamDESimplifier* _param): param(_param) {}

  void initialize(const SMeshT& original);
  void simplify();
private:
  std::unique_ptr<DERemesher> remesher;
  std::unique_ptr<DEDM> delaunay_mesher;

  // -- for error bound
  double original_diagonal_length;
  double max_distance_error_value;
  Vec3d p_epsilon;
  void calc_max_distance_error_value();
};

}// namespace DE
}// namespace Simplification
}// namespace GCLF