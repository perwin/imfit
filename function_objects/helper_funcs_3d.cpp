/* FILE: helper_funcs_3d.cpp ------------------------------------------- */
/* 
 * Compute local object coordinates x_obj,y_obj,z_obj corresponding to a location
 * (x_d,y_d,z_d) which is at line-of-sight distance s from start point (x_d0, y_d0, z_d0), 
 * where the midplane of the object is oriented at angle (90 - inclination) to the line 
 * of sight vector.
 *
 * x_d,y_d,z_d = 3d Cartesian coordinates in frame aligned with object midplane,
 * with x-axis = line of nodes for object midplane (i.e., intersection of object
 * midplane with the sky). Corresponds to position s along the line-of-sight ray for
 * the current pixel.
 *
 * x_d0,y_d0,z_d0 = same, but corresponding to intersection of line-of-sight ray
 * with the sky plane (s = 0).
 *
 * x_obj,y_obj,z_obj = 3d Cartesian coordinates in frame aligned with object midplane,
 * with x-axis = major axis of object.
*/

// Copyright 2011--2017 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.

#include <math.h>

#include "helper_funcs_3d.h"

void Compute3dObjectCoords( double s, double x_d0, double y_d0, double z_d0, 
							double sinInc, double cosInc, double cosObjPA, double sinObjPA,
							double & x_obj, double & y_obj, double & z_obj )
{
  // Given current value of s and the pre-defined parameters, determine our 
  // 3D location (x_d,y_d,z_d) [by construction, x_d = x_d0]
  double y_d = y_d0 + s*sinInc;
  double z_d = z_d0 - s*cosInc;
  
  // Convert 3D Cartesian coordinate to rotated x_obj,y_obj,z_obj coordinate,
  // where x_obj is along object's major axis, y_obj is perpendicular in object
  // midplane, and z_obj is perpendicular to object midplane.
  x_obj = x_d0*cosObjPA + y_d*sinObjPA;
  y_obj = -x_d0*sinObjPA + y_d*cosObjPA;
  z_obj = fabs(z_d);
}


/* END OF FILE: helper_funcs_3d.cpp ------------------------------------ */
