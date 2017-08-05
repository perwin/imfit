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
#include <tuple>

#include "helper_funcs_3d.h"



std::tuple<double, double, double> ImageCoordsTo3dStartCoords( double x_diff, double y_diff,
															double cosPA, double sinPA, 
															double cosInc, double sinInc )
{
  // Calculate x,y in component (projected sky) reference frame, corrected for
  // rotation of line of nodes
  double xp = x_diff*cosPA + y_diff*sinPA;
  double yp = -x_diff*sinPA + y_diff*cosPA;

  // Calculate (x,y,z)_start point in component's axisymmetric xyz reference frame, 
  // corresponding to intersection of line-of-sight ray with projected sky frame
  double x_d0 = xp;
  double y_d0 = yp * cosInc;
  double z_d0 = yp * sinInc;
  return std::make_tuple(x_d0, y_d0, z_d0);
}


std::tuple<double, double, double> Compute3dObjectCoords( double s, double x_d0, double y_d0, 
												double z_d0, double sinInc, double cosInc, 
												double cosObjPA, double sinObjPA )
{
  // Given the 3D start coordinates, the current value of s and the pre-defined 
  // parameters, determine our current 3D location (x_d,y_d,z_d) [by construction, 
  // x_d = x_d0]
  double y_d = y_d0 + s*sinInc;
  double z_d = z_d0 - s*cosInc;
  
  // Convert 3D "axisymmetric" coordinate to rotated x_obj,y_obj,z_obj coordinate,
  // where x_obj is defined by object's major axis, y_obj is perpendicular in object
  // midplane, and z_obj is perpendicular to object midplane.
  double x_obj = x_d0*cosObjPA + y_d*sinObjPA;
  double y_obj = -x_d0*sinObjPA + y_d*cosObjPA;
  double z_obj = fabs(z_d);
  return std::make_tuple(x_obj, y_obj, z_obj);
}


// Compute equivalent radius given position (xp, yp, zp) relative to center of
// object for triaxial ellipsoid with super-quadric isodensity contours
double CalculateTriaxEquivRadius( double xp, double yp, double zp, double q, double q_z, 
									double c_xy, double c_z )
{
  double  x_to_cxy, y_scaled_to_cxy, z_scaled_to_cz, r_cz;

  x_to_cxy = pow(fabs(xp), c_xy);
  y_scaled_to_cxy = pow(fabs(yp/q), c_xy);
  z_scaled_to_cz = pow(fabs(zp/q_z), c_z);
  r_cz = pow(x_to_cxy + y_scaled_to_cxy, c_z / c_xy) + z_scaled_to_cz;

  return pow(r_cz, 1.0/c_z);
}



// r2 = r^2
// twosigma_squared = 2 * sigma^2
double LumDensity_Gaussian( double r2, double J_0, double twosigma_squared )
{
  return J_0 * exp(-r2/twosigma_squared);
}


// Generalized Gaussian function (equivalent to Sersic function, but using
// radial scale length alpha instead of half-light radius r_e
//    beta = 1/n
//    inv_alpha = 1.0 / alpha
// For Gaussian, inv_alpha = 1/(sqrt(2)*sigma), beta = 2 (n = 1/2)
double LumDensity_GenGaussian( double r, double J_0, double inv_alpha, double beta )
{
  double X = pow(r * inv_alpha, beta);
  return J_0 * exp(-X);
}


/* END OF FILE: helper_funcs_3d.cpp ------------------------------------ */
