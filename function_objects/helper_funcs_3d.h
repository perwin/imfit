
#ifndef _3D_HELPER_H_
#define _3D_HELPER_H_

#include <tuple>

// void Compute3dObjectCoords( double s, double x_d0, double y_d0, double z_d0, 
// 							double sinInc, double cosInc, double cosObjPA, double sinObjPA,
// 							double & x_obj, double & y_obj, double & z_obj );
std::tuple<double, double, double> Compute3dObjectCoords( double s, double x_d0, double y_d0, 
												double z_d0, double sinInc, double cosInc, 
												double cosObjPA, double sinObjPA );


#endif  // _3D_HELPER_H_
