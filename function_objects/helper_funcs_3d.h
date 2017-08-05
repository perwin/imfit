
#ifndef _3D_HELPER_H_
#define _3D_HELPER_H_

#include <tuple>

// Coordinate-transformation functions

std::tuple<double, double, double> ImageCoordsTo3dStartCoords( double x_diff, double y_diff,
															double cosPA, double sinPA, 
															double cosInc, double sinInc );

std::tuple<double, double, double> Compute3dObjectCoords( double s, double x_d0, double y_d0, 
												double z_d0, double sinInc, double cosInc, 
												double cosObjPA, double sinObjPA );

double CalculateTriaxEquivRadius( double xp, double yp, double zp, double q, double q_z, 
									double c_xy, double c_z );


// Luminosity-density functions

double LumDensity_Gaussian( double r2, double J_0, double twosigma_squared );

double LumDensity_GenGaussian( double r, double J_0, double inv_alpha, double beta );


#endif  // _3D_HELPER_H_
