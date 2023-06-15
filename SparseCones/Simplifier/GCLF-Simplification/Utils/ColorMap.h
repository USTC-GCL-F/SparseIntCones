#pragma once
#include "Basic/Types.h"

namespace GCLF
{
namespace Simplification
{
namespace Utils
{

using namespace Geometry;

class ColorMap
{
private:
	double maxValue, minValue;
public:
	ColorMap();
	ColorMap(double maxV, double minV);

	double GetMaxValue();
	double GetMinValue();
	void SetMaxMinValue(double maxV, double minV);

	Vec3d MapToColor(double value);
};
}// namespace Utils
}// namespace Simplification
}// namespace GCLF