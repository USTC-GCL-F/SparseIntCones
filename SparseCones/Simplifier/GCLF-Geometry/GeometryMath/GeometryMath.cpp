#include "GeometryMath.h"
#include "indirect_predicates.h"

namespace GCLF
{
namespace Geometry
{

/**********************/
/* Geometry threshold */
/**********************/

double GeomThr::length_de_thr = 1e-6;
double GeomThr::area_de_thr = 1e-10;
double GeomThr::cos_de_thr = 0.965;
double GeomThr::radian_de_thr = 1e-4;
double GeomThr::edge_length_ratio_de_thr = 0.2;
double GeomThr::relocate_length_thr = 1e-4;

/***********/
/*  Angle  */
/***********/

double angle_from_three_points(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3)
{
	if ((p1 - p2).length() < GeomThr::length_de_thr
		|| (p1 - p3).length() < GeomThr::length_de_thr
		|| (p2 - p3).length() < GeomThr::length_de_thr)
		return 0.0;
	Vec3d v1 = (p1 - p2).normalized();
	Vec3d v2 = (p3 - p2).normalized();
	double  cos_value = v1 | v2;
	cos_value = std::max(std::min(cos_value, 1.0), -1.0);
	return std::acos(cos_value);
}

int find_obtuse_angle(const Triangle& triangle)
{
  for (int i = 0; i < 3; i++)
  {
    if (dot_sign(
      triangle[(i + 1) % 3] - triangle[i],
      triangle[(i + 2) % 3] - triangle[i]) == -1)
    {
      return i;
    }
  }
  return -1;
}

/// @brief find which octant the vec lies in.
int find_octant(const Vec3d& vec)
{
  static int octant[2][2][2] = { {{6, 2}, {5, 1}}, {{7, 3},{4, 0}} };
  int xsign = vec[0] >= 0 ? 1 : 0;
  int ysign = vec[1] >= 0 ? 1 : 0;
  int zsign = vec[2] >= 0 ? 1 : 0;
  return octant[xsign][ysign][zsign];
}

/*****************/
/*     Area      */
/*****************/

double triangle_area(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c)
{
	return 0.5 * ((a - b) % (c - b)).length();
}

/*****************/
/*     Volume    */
/*****************/

/// calculate tet volume
double tet_volume(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c, CR_Vec3d d)
{
	// detect d is on the tri, or below the tri, or in the tri
	double temp[3][3];
	for (int i = 0; i < 3; ++i)
	{
		temp[0][i] = a[i] - d[i];
		temp[1][i] = b[i] - d[i];
		temp[2][i] = c[i] - d[i];
	}
	return std::abs(temp[0][0] * temp[1][1] * temp[2][2]
		+ temp[0][1] * temp[1][2] * temp[2][0]
		+ temp[0][2] * temp[1][0] * temp[2][1]
		- temp[0][2] * temp[1][1] * temp[2][0]
		- temp[0][1] * temp[1][0] * temp[2][2]
		- temp[0][0] * temp[1][2] * temp[2][1]);
}

/******************************/
/*     Coordinate system      */
/******************************/

Vec3d calc_vertical_vec2vec(const Vec3d& vec)
{
  if (vec[0] != 0.0)
  {
    return Vec3d(-(vec[1] + vec[2]) / vec[0], 1.0, 1.0).normalized();
  }
  else	// vec must be in y-z-plane or z-axis
  {
    return Vec3d(1.0, 0, 0);
  }
}

void make_coordinate_system(const Vec3d& axis_z, Vec3d& axis_x, Vec3d& axis_y)
{
	Vec3d temp_vec;
	if (abs(axis_z.x()) < 0.9)
		temp_vec = Vec3d(1., 0., 0.);
	else if (abs(axis_z.y()) < 0.9)
		temp_vec = Vec3d(0., 1., 0.);
	else if (abs(axis_z.z()) < 0.9)
		temp_vec = Vec3d(0., 0., 1.);

	axis_x = axis_z.cross(temp_vec).normalized();
	axis_y = axis_z.cross(axis_x).normalized();
}
}// namespace Geometry
}// namespace GCLF