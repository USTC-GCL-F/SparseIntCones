#pragma once
#include "DE/DESimplifier.h"

namespace GCLF
{
namespace Simplification
{
using namespace DE;

inline void simplify(const SMeshT& src, SMeshT& dst, size_t target_vn)
{
  Logger::InitLogger();

  ParamDESimplifier param;
  param.errorRatio = 0.1;
  param.paramRemesher_DE.collapseIter = 5;
  param.paramRemesher_DE.splitIter = 5;
  param.paramRemesher_DE.smoothIter = 3;
  param.paramRemesher_DE.equalizeValenceIter = 3;
  param.paramRemesher_DE.enlargeTargetLengthRatio = 1.1;
  param.paramRemesher_DE.target_vn = target_vn;
  param.paramDM_DE.Np = 50;
  param.paramDM_DE.mutationScale = 0.5;
  param.paramDM_DE.crossoverRate = 0.9;
  param.paramDM_DE.target_vn = target_vn;

  DESimplifier simplifier(&param);

  simplifier.initialize(src);
  simplifier.simplify();

  dst = *simplifier.rm;
}

}// namespace Simplification
}// namespace GCLF