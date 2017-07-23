// Copyright 2017 by Peter Erwin.
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
#include <string>
#include <tuple>

#include "image_io.h"
#include "options_base.h"
#include "getimages.h"
#include "fftw3.h"  // so we can call fftw_free()


// Code which retrieves an image from a (FITS) file and checks its dimensions
// against a reference
std::tuple<double *, int> GetAndCheckImage( const string imageName, const string imageType,
											int nColumns_ref, int nRows_ref )
{
  int  nColumns = 0;
  int  nRows = 0;
  double *imagePixels = nullptr;

  /* Get and check mask image */
  imagePixels = ReadImageAsVector(imageName, &nColumns, &nRows);
  if (imagePixels == NULL) {
    fprintf(stderr,  "\n*** ERROR: Unable to read %s file \"%s\"!\n\n", 
    			imageType.c_str(), imageName.c_str());
    return std::make_tuple(imagePixels, -1);
  }
  if ( ((nColumns_ref > 0) && (nRows_ref > 0)) && ((nColumns != nColumns_ref) ||
  		(nRows != nRows_ref)) ) {
    fprintf(stderr, "\n*** ERROR: Dimensions of %s image (%s: %d columns, %d rows)\n",
            imageType.c_str(), imageName.c_str(), nColumns, nRows);
    fprintf(stderr, "do not match dimensions of data image (%d columns, %d rows)!\n\n",
            nColumns_ref, nRows_ref);
    fftw_free(imagePixels);
    return std::make_tuple(imagePixels, -2);
  }

  return std::make_tuple(imagePixels, 0);
}


// Function which retrieves and checks dimensions for mask and/or noise/error images.
// Return values:
// 		status = 1: mask image loaded, but no error image specified
// 		status = 2: error image loaded, but no mask image was specified
// 		status = 3: both images specified & loaded
std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
										OptionsBase *options, bool &maskPixelsAllocated, 
										bool &errorPixelsAllocated )
{
  int  status = 0;
  int  returnVal = 0;
  double *maskPixels = nullptr;
  double *errorPixels = nullptr;
  
  maskPixelsAllocated = false;
  errorPixelsAllocated = false;

  /* Get and check mask image */
  if (options->maskImagePresent) {
    printf("Reading mask image (\"%s\") ...\n", options->maskFileName.c_str());
    std::tie(maskPixels, status) = GetAndCheckImage(options->maskFileName, "mask",
    												nColumns, nRows);
    if (status < 0)
      return std::make_tuple(maskPixels, errorPixels, -1);
    maskPixelsAllocated = true;
    returnVal += 1;
  }
  /* Get and check error image, if supplied */
  if (options->noiseImagePresent) {
    printf("Reading noise image (\"%s\") ...\n", options->noiseFileName.c_str());
    std::tie(errorPixels, status) = GetAndCheckImage(options->noiseFileName, "noise",
    												nColumns, nRows);
    if (status < 0)
      return std::make_tuple(maskPixels, errorPixels, -1);
    errorPixelsAllocated = true;
    returnVal += 2;
  }
  
  return std::make_tuple(maskPixels, errorPixels, returnVal);
}


// Function which reads and returns data corresponding to requested PSF image
std::tuple<double *, int, int, int> GetPsfImage( OptionsBase *options )
{
  int  status;
  int  nColumns_psf, nRows_psf;
  double *psfPixels = nullptr;
  
  // Read in PSF image
  printf("Reading PSF image (\"%s\") ...\n", options->psfFileName.c_str());
  std::tie(psfPixels, status) = GetAndCheckImage(options->psfFileName.c_str(), "PSF", 0,0);
  if (status < 0)
    return std::make_tuple(psfPixels, 0,0, -1);

  GetImageSize(options->psfFileName, &nColumns_psf, &nRows_psf);
  long nPixels_psf = (long)nColumns_psf * (long)nRows_psf;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
         nColumns_psf, nRows_psf, nPixels_psf);

  return std::make_tuple(psfPixels, nColumns_psf, nRows_psf, 0);
}

