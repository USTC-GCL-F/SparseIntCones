#pragma once
#include <vector>
#include <string>
#include "boost/json.hpp"

namespace GCLF
{
namespace Simplification
{
namespace DE
{
/***********************/
/* Parameters for DE   */
/***********************/

struct DEParamRemesher
{
  // remeshing parameters
  size_t collapseIter;
  size_t splitIter;
  size_t equalizeValenceIter;
  size_t smoothIter;
  double enlargeTargetLengthRatio;
  size_t target_vn;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["collapseIter"] = collapseIter;
    jo["splitIter"] = splitIter;
    jo["equalizeValenceIter"] = equalizeValenceIter;
    jo["smoothIter"] = smoothIter;
    jo["enlargeTargetLengthRatio"] = enlargeTargetLengthRatio;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    collapseIter = jo.at("collapseIter").as_int64();
    splitIter = jo.at("splitIter").as_int64();
    equalizeValenceIter = jo.at("equalizeValenceIter").as_int64();
    smoothIter = jo.at("smoothIter").as_int64();
    enlargeTargetLengthRatio = jo.at("enlargeTargetLengthRatio").as_double();
  }
};

struct ParamDEDM
{
  // -- for DE
	size_t Np;
	double mutationScale;
	double crossoverRate;
  size_t target_vn;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["Np"] = Np;
    jo["mutationScale"] = mutationScale;
    jo["crossoverRate"] = crossoverRate;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    Np = jo.at("Np").as_int64();
    mutationScale = jo.at("mutationScale").as_double();
    crossoverRate = jo.at("crossoverRate").as_double();
  }
};

struct ParamDESimplifier
{
  double errorRatio;
  DEParamRemesher paramRemesher_DE;
  ParamDEDM paramDM_DE;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["errorRatio"] = errorRatio;
    jo["paramRemesher_DE"] = paramRemesher_DE.serialize();
    jo["paramDM_DE"] = paramDM_DE.serialize();
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    errorRatio = jo.at("errorRatio").as_double();
    paramRemesher_DE.deserialize(jo.at("paramRemesher_DE").as_object());
    paramDM_DE.deserialize(jo.at("paramDM_DE").as_object());
  }
};
}// namespace DE
}// namespace Simplification
}// namespace GCLF