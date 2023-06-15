#pragma once
#include <cmath>
#include "Basic/Types.h"
#include "Basic/Triangle.h"
#include "Exact/Predicates.h"

namespace GCLF
{
namespace Geometry
{
typedef const Vec3d& CR_Vec3d;

/**********************/
/* Geometry threshold */
/**********************/

struct GeomThr
{
  /// used in:
  /// * angle_from_three_points
  static double length_de_thr;
  /// used in:
  /// * check_wrinkle
  /// * flip_ok
  static double area_de_thr;
  /// used in:
  /// * eliminate_almost_degeneration
  static double cos_de_thr;
  /// used in:
  /// * check_wrinkle
  static double radian_de_thr;
  /// used in:
  /// * eliminate_almost_degeneration
  static double edge_length_ratio_de_thr;
  /// used in:
  /// * tangential relaxation
  static double relocate_length_thr;
};

/*****************/
/*    Angle      */
/*****************/
double angle_from_three_points(const Vec3d& p1, const Vec3d& p2, const Vec3d& p3);

int find_obtuse_angle(const Triangle& triangle);

int find_octant(const Vec3d& vec);

/*****************/
/*     Area      */
/*****************/

/// calculate tri area
double triangle_area(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c);

/*****************/
/*     Volume    */
/*****************/

/// calculate tet volume
double tet_volume(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c, CR_Vec3d d);

/******************************/
/*     Coordinate system      */
/******************************/

Vec3d calc_vertical_vec2vec(const Vec3d& vec);

void make_coordinate_system(const Vec3d& axis_z, Vec3d& axis_x, Vec3d& axis_y);

}// namespace Geometry
}// namespace GCLF