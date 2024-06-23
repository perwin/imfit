/* FILE: utilities_multimfit.cpp ----------------------------------- */
/*   Utility routines used by multimfit.
 */

// Copyright 2018 by Peter Erwin.
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

#include <stdio.h>
#include <vector>
#include <tuple>

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"
#include "model_object_multimage.h"
#include "utilities_multimfit.h"


/// Corrects X0,Y0 parameter initial values and limits for image-subsection-based
/// offsets (including both image-description parameters and main model parameters).
void CorrectForImageOffsets( double *paramVect, vector<mp_par> &paramLimits, 
							ModelObjectMultImage *multImageModel )
{
  int  nParamsTot = multImageModel->GetNParams();
  int  nImages = multImageModel->GetNImages();
  int  nImageParams = N_IMAGE_PARAMS * (nImages - 1);
  int  xOffset_this, yOffset_this, xOffset_ref, yOffset_ref;
  int  i_x0, i_y0;
  
//   for(int i = 0; i < nParamsTot; i++)
//     printf("i = %d: name = %s\n", i, multImageModel->GetParameterName(i).c_str());
    
  // Apply offset fixes to image-description parameter X0,Y0 pairs, for each 
  // extra image
  for (int n = 1; n < nImages; n++) {
    std::tie(xOffset_this,yOffset_this) = multImageModel->GetOffsetsForImage(n);
//     printf("Offsets for image %d: %d,%d\n", n, xOffset_this,yOffset_this);
    // assume that X0,Y0 are last two parameters for each image
    i_x0 = n*N_IMAGE_PARAMS - 2;
//     printf("  i_x0 = %d: paramVect[i_x0][i] = %.1f\n", i_x0, paramVect[i_x0]);
    paramVect[i_x0] -= xOffset_this;
    paramLimits[i_x0].offset = xOffset_this;
    if (paramLimits[i_x0].limited[0] == 1)
      paramLimits[i_x0].limits[0] -= xOffset_this;
    if (paramLimits[i_x0].limited[1] == 1)
      paramLimits[i_x0].limits[1] -= xOffset_this;
//     printf("     paramVect[i_x0][i] = %.1f\n", paramVect[i_x0]);
    i_y0 = i_x0 + 1;
//     printf("  i_y0 = %d: paramVect[i_y0][i] = %.1f\n", i_y0, paramVect[i_y0]);
    paramVect[i_y0] -= yOffset_this;
    paramLimits[i_y0].offset = yOffset_this;
    if (paramLimits[i_y0].limited[0] == 1)
      paramLimits[i_y0].limits[0] -= yOffset_this;
    if (paramLimits[i_y0].limited[1] == 1)
      paramLimits[i_y0].limits[1] -= yOffset_this;
//     printf("     paramVect[i_y0][i] = %.1f\n", paramVect[i_y0]);
  }
  
  // Apply reference-image offsets to main model's function-block X0,Y0 sets
  // (reference-image coordinates)
  std::tie(xOffset_ref,yOffset_ref) = multImageModel->GetOffsetsForImage(0);
//   printf("Offsets for reference image: %d,%d\n", xOffset_ref,yOffset_ref);
  for (int i = nImageParams; i < nParamsTot; i++) {
    if (multImageModel->GetParameterName(i) == X0_string) {
//       printf("[X0:]     paramVect[%d] = %.1f\n", i, paramVect[i]);
      paramVect[i] -= xOffset_ref;
//       printf("     updated paramVect[%d] = %.1f\n", i, paramVect[i]);
      paramLimits[i].offset = xOffset_ref;
      if (paramLimits[i].limited[0] == 1)
        paramLimits[i].limits[0] -= xOffset_ref;
      if (paramLimits[i].limited[1] == 1)
        paramLimits[i].limits[1] -= xOffset_ref;
    } else if (multImageModel->GetParameterName(i) == Y0_string) {
//       printf("[Y0:]     paramVect[%d] = %.1f\n", i, paramVect[i]);
      paramVect[i] -= yOffset_ref;
//       printf("     updated paramVect[%d] = %.1f\n", i, paramVect[i]);
      paramLimits[i].offset = yOffset_ref;
      if (paramLimits[i].limited[0] == 1)
        paramLimits[i].limits[0] -= yOffset_ref;
      if (paramLimits[i].limited[1] == 1)
        paramLimits[i].limits[1] -= yOffset_ref;
    }
  }
}



/* END OF FILE: utilities_multimfit.cpp ---------------------------- */

