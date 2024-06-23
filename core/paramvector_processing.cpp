// Code for parsing and re-arranging parameter vectors in ModelObjectMultImage

// Copyright 2017-2020 by Peter Erwin.
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
#include <stdio.h>
#include <stdlib.h>
#include <tuple>

#include "definitions_multimage.h"
#include "paramvector_processing.h"

using namespace std;

const double  DEG2RAD = 0.017453292519943295;



// Given a global parameter vector (e.g., as supplied to ModelObjectMultImage
// by an external caller), extract the image parameters for the specified
// imageNumber (0 = first [reference] image, 1 = second image, etc.).
// Note that this is actually only valid for the second and subsequent images!
//
// Currently used in param_holder.cpp *and* model_object_multimage.cpp
// Python prototype: parameter_munging.ExtractImageParams
int ExtractImageParams( double inputParamsVector[], int imageNumber, 
						double& pixScale, double& rotation, double& intensityScale,
    					double& X0, double& Y0 )
{
  // if imageNumber == 0 [first image], we return without modifying any of
  // the output parameters (pixScale, etc.)
  if (imageNumber <= 0)
    return -1;
    
  // image *numbering* starts with 0, but image parameters in inputParamsVector
  // only exist for 2nd and subsequent images
  int nn = imageNumber - 1;
  
  pixScale = inputParamsVector[N_IMAGE_PARAMS*nn];
  rotation = inputParamsVector[N_IMAGE_PARAMS*nn + 1];
  intensityScale = inputParamsVector[N_IMAGE_PARAMS*nn + 2];
  X0 = inputParamsVector[N_IMAGE_PARAMS*nn + 3];
  Y0 = inputParamsVector[N_IMAGE_PARAMS*nn + 4];
  
  return 0;
}



// 	(non-reference) image: X0_n_im, Y0_n_im, given:
//
// 	X0_0_ref, Y0_0_ref = X0,Y0 of first function block in reference image
// 	X0_n_ref, Y0_n_ref = X0,Y0 of function block n in reference image
// 	X0_0_im, Y0_0_im = X0,Y0 of first function block in current image
// 	pixScale_im = pixel scale in current image *relative to* reference image
// 	rotation_im = rotation of current image *relative to* reference image
// 		(i.e., rotation of current-image y-axis relative to reference-image
// 		y-axis, in degrees CCW)
// 		
// 	Return values X0_out, Y0_out, status = X0_n_im,Y0_n_im
//
// Currently used in param_holder.cpp
// Python prototype: parameter_munging.CalculateOffsetX0Y0
std::tuple<double, double, int> CalculateOffset_X0Y0( double X0_0_ref, double Y0_0_ref, 
								double X0_n_ref, double Y0_n_ref, double X0_0_im, 
								double Y0_0_im, double pixScale_im, double rotation_im )
{
  if (pixScale_im <= 0) {
    fprintf(stderr, "   CalculateOffset_X0Y0: Zero or negative pixel-scale parameter!\n");
    return std::make_tuple(0.0, 0.0, -1);
  }
  
  double rotation_rad = DEG2RAD*rotation_im;
  double cosTheta = cos(rotation_rad);
  double sinTheta = sin(rotation_rad);
  // calculate offsets in reference image
  double deltaX_n_ref = X0_n_ref - X0_0_ref;
  double deltaY_n_ref = Y0_n_ref - Y0_0_ref;
  // transform to rotated coordinate system with same center
  double deltaX_n_im = deltaX_n_ref*cosTheta + deltaY_n_ref*sinTheta;
  double deltaY_n_im = -deltaX_n_ref*sinTheta + deltaY_n_ref*cosTheta;

  // transform scaled offsets to current image
  double X0_out = X0_0_im + pixScale_im*deltaX_n_im;
  double Y0_out = Y0_0_im + pixScale_im*deltaY_n_im;
  
  return std::make_tuple(X0_out, Y0_out, 0);
}


// FIXME: Is this superfluous? (duplicated by method in ParamHolder?)
// Python prototype: parameter_munging.AssembleParams
void AssembleParamsForImage( double externalInputParamsVect[], int nInputParamsTot,
							int nImagesTot, int imageNumber, int nFunctions, 
							vector<int>& paramSizes, bool fblockStartFlags[], 
							int nFuncBlocks, double outputModelParams[] )
{
  int i;
  double pixScale_im, rot_im, iScale;
  double X0_0_im, Y0_0_im;
  // 5 image parameters for each of the *second and subsequent* images
  int nImageParams = N_IMAGE_PARAMS*(nImagesTot - 1);
  int nModelParams = nInputParamsTot - nImageParams;  // number of params for ModelObject

  // Note that model parameters start at externalInputParamsVect[nImageParams]
  // (i.e., skipping over the initial image-specific parameters)
  for (i = 0; i < nModelParams; i++)
    outputModelParams[i] = externalInputParamsVect[nImageParams + i];

  if (imageNumber > 0) {
    ExtractImageParams(externalInputParamsVect, imageNumber, pixScale_im, rot_im, iScale, 
    					X0_0_im, Y0_0_im);
    // start updating copy of main vector with proper X0,Y0 for this image
    outputModelParams[0] = X0_0_im;
    outputModelParams[1] = Y0_0_im;
    // now compute offsets and second and subsequent function-block X0,Y0 values
    if (nFuncBlocks > 1) {
      double X0_0_ref = externalInputParamsVect[nImageParams];
      double Y0_0_ref = externalInputParamsVect[nImageParams + 1];
      // skip the first function (X0,Y0 + paramSizes[0])
      int offset = paramSizes[0] + 2;
      for (int n = 1; n < nFunctions; n++) {
        if (fblockStartFlags[n] == true) {
          // new function block, so we calculate and store offset/rotated/scaled
          // X0,Y0 for this function block in this image
          double X0_n_ref = externalInputParamsVect[nImageParams + offset];
          double Y0_n_ref = externalInputParamsVect[nImageParams + offset + 1];
          double X0_n_im, Y0_n_im;
          std::tie(X0_n_im, Y0_n_im, std::ignore) = CalculateOffset_X0Y0(X0_0_ref, Y0_0_ref, 
          											X0_n_ref, Y0_n_ref, X0_0_im, Y0_0_im,
          											pixScale_im, rot_im);
          outputModelParams[offset] = X0_n_im;
          outputModelParams[offset + 1] = Y0_n_im;
//           printf("  n = %d: offset = %d, X0_n_im = %.1f, Y0_n_im = %.1f\n",
//           		n, offset, X0_n_im, Y0_n_im);
          offset += 2;
        }
        offset += paramSizes[n];
      }
    }
  }

}
