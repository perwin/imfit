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
#include "psf_oversampling_info.h"
#include "fftw3.h"  // so we can call fftw_free()


/// Main utility function: reads in image from FITS file imageName and checks dimensions
/// against reference values (if latter are nonzero). Returns (nullptr, -1) if unable
/// to read image-data from file; returns -2 (and frees allocated memory) if image 
/// dimensions do not match references dimensions (unless latter are both 0).
std::tuple<double *, int> GetAndCheckImage( const string imageName, const string imageType,
											int nColumns_ref, int nRows_ref )
{
  int  nColumns = 0;
  int  nRows = 0;
  double *imagePixels = nullptr;

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
    imagePixels = nullptr;
    return std::make_tuple(imagePixels, -2);
  }

  return std::make_tuple(imagePixels, 0);
}



/// Function which retrieves and checks dimensions for mask and/or noise/error images.
/// Returns tuple of (maskPixels, errorPixels, status), where maxPixels and errorPixels
/// are double * (and are = nullptr when image in question was not requested).
/// In case of errors in retrieving image data -- or if image dimensions do not
/// match reference dimensions nColumns, nRows -- then (nullptr, nullptr, -1) is
/// returned.
/// Return values:
/// 		status = 1: mask image loaded, but no error image specified
/// 		status = 2: error image loaded, but no mask image was specified
/// 		status = 3: both images specified & loaded
std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
										string &maskFileName, string &errorFileName, 
										bool &maskPixelsAllocated, bool &errorPixelsAllocated )
{
  int  status = 0;
  int  returnVal = 0;
  double *maskPixels = nullptr;
  double *errorPixels = nullptr;
  
  maskPixelsAllocated = false;
  errorPixelsAllocated = false;

  /* Get and check mask image */
  if (maskFileName.size() > 0) {
    printf("Reading mask image (\"%s\") ...\n", maskFileName.c_str());
    std::tie(maskPixels, status) = GetAndCheckImage(maskFileName, "mask", nColumns, nRows);
    if (status < 0)
      return std::make_tuple(maskPixels, errorPixels, -1);
    maskPixelsAllocated = true;
    returnVal += 1;
  }
  /* Get and check error image, if supplied */
  if (errorFileName.size() > 0) {
    printf("Reading noise image (\"%s\") ...\n", errorFileName.c_str());
    std::tie(errorPixels, status) = GetAndCheckImage(errorFileName, "noise", nColumns, nRows);
    if (status < 0)
      return std::make_tuple(maskPixels, errorPixels, -1);
    errorPixelsAllocated = true;
    returnVal += 2;
  }
  
  return std::make_tuple(maskPixels, errorPixels, returnVal);
}

// std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
// 										const OptionsBase *options, bool &maskPixelsAllocated, 
// 										bool &errorPixelsAllocated )
// {
//   int  status = 0;
//   int  returnVal = 0;
//   double *maskPixels = nullptr;
//   double *errorPixels = nullptr;
//   
//   maskPixelsAllocated = false;
//   errorPixelsAllocated = false;
// 
//   /* Get and check mask image */
//   if (options->maskImagePresent) {
//     printf("Reading mask image (\"%s\") ...\n", options->maskFileName.c_str());
//     std::tie(maskPixels, status) = GetAndCheckImage(options->maskFileName, "mask",
//     												nColumns, nRows);
//     if (status < 0)
//       return std::make_tuple(maskPixels, errorPixels, -1);
//     maskPixelsAllocated = true;
//     returnVal += 1;
//   }
//   /* Get and check error image, if supplied */
//   if (options->noiseImagePresent) {
//     printf("Reading noise image (\"%s\") ...\n", options->noiseFileName.c_str());
//     std::tie(errorPixels, status) = GetAndCheckImage(options->noiseFileName, "noise",
//     												nColumns, nRows);
//     if (status < 0)
//       return std::make_tuple(maskPixels, errorPixels, -1);
//     errorPixelsAllocated = true;
//     returnVal += 2;
//   }
//   
//   return std::make_tuple(maskPixels, errorPixels, returnVal);
// }


/// Function which reads and returns data corresponding to requested PSF image,
/// along with PSF image dimensions:
/// tuple of (psfPixels, nColumns_psf, nRows_psf, status).
/// In case of errors in retrieving image data, return value is
/// (nullptr, 0, 0, -1)
std::tuple<double *, int, int, int> GetPsfImage( const string &psfFileName )
{
  int  status;
  int  nColumns_psf, nRows_psf;
  double *psfPixels = nullptr;
  
  // Read in PSF image
  printf("Reading PSF image (\"%s\") ...\n", psfFileName.c_str());
  std::tie(psfPixels, status) = GetAndCheckImage(psfFileName.c_str(), "PSF", 0,0);
  if (status < 0)
    return std::make_tuple(psfPixels, 0,0, -1);

  GetImageSize(psfFileName, &nColumns_psf, &nRows_psf);
  long nPixels_psf = (long)nColumns_psf * (long)nRows_psf;
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
         nColumns_psf, nRows_psf, nPixels_psf);

  return std::make_tuple(psfPixels, nColumns_psf, nRows_psf, 0);
}

// std::tuple<double *, int, int, int> GetPsfImage( const OptionsBase *options )
// {
//   int  status;
//   int  nColumns_psf, nRows_psf;
//   double *psfPixels = nullptr;
//   
//   // Read in PSF image
//   printf("Reading PSF image (\"%s\") ...\n", options->psfFileName.c_str());
//   std::tie(psfPixels, status) = GetAndCheckImage(options->psfFileName.c_str(), "PSF", 0,0);
//   if (status < 0)
//     return std::make_tuple(psfPixels, 0,0, -1);
// 
//   GetImageSize(options->psfFileName, &nColumns_psf, &nRows_psf);
//   long nPixels_psf = (long)nColumns_psf * (long)nRows_psf;
//   printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
//          nColumns_psf, nRows_psf, nPixels_psf);
// 
//   return std::make_tuple(psfPixels, nColumns_psf, nRows_psf, 0);
// }



// Function which reads and returns data corresponding to requested oversampled PSF
// images
int GetOversampledPsfInfo( const OptionsBase *options, int xOffset, int yOffset, 
							vector<PsfOversamplingInfo *> &psfOversamplingInfoVect )
{
  PsfOversamplingInfo * psfOversampleInfo;
  double * psfOversampledPixels;
  int  nColumns_psf_oversampled, nRows_psf_oversampled;
  long  nPixels_psf_oversampled;
  int psfOversamplingScale;
  
  int nOversampledPsfImages = (int)options->psfOversampledFileNames.size();
  int nOversampledScales = (int)options->psfOversamplingScales.size();
  if (nOversampledPsfImages != nOversampledScales) {
    fprintf(stderr, "\n*** ERROR: number of oversampling scales (%d) is not the same\n", 
   					nOversampledScales);
    fprintf(stderr, "           as number of oversampled-PSF images (%d)!\n\n",
   					nOversampledPsfImages);
    return -1;
  }
  if ((nOversampledPsfImages > 1) && (nOversampledPsfImages != options->nOversampleRegions)) {
    fprintf(stderr, "\n*** ERROR: number of oversampled-PSF images (%d) must be = 1 OR\n", 
   					nOversampledPsfImages);
    fprintf(stderr, "           must be same as number of oversampled-PSF regions (%d)!\n\n",
   					options->nOversampleRegions);
    return -1;
  }
  for (int nn = 0; nn < options->nOversampleRegions; nn++) {
    psfOversampleInfo = new PsfOversamplingInfo();
    psfOversampleInfo->AddImageOffset(xOffset, yOffset);
    psfOversampleInfo->AddRegionString(options->psfOversampleRegions[nn]);
    bool newPsfOversampledPixelsFlag = false;
    if ( (nn == 0) || ((nn > 0) && (nOversampledPsfImages > 1)) ) {
      // Always read PSF image and get oversampling scale from options object the
      // first time through; do it again if user supplied more than one image and
      // scale (otherwise, we reuse the same image and scale for subsequent regions)
      printf("Reading oversampled PSF image (\"%s\") ...\n", options->psfOversampledFileNames[nn].c_str());
      psfOversampledPixels = ReadImageAsVector(options->psfOversampledFileNames[nn], 
	    							&nColumns_psf_oversampled, &nRows_psf_oversampled);
      if (psfOversampledPixels == NULL) {
        fprintf(stderr, "\n*** ERROR: Unable to read oversampled PSF image file \"%s\"!\n\n", 
		      			options->psfOversampledFileNames[nn].c_str());
        return -1;
      }
      newPsfOversampledPixelsFlag = true;
      nPixels_psf_oversampled = (long)nColumns_psf_oversampled * (long)nRows_psf_oversampled;
      printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
             nColumns_psf_oversampled, nRows_psf_oversampled, nPixels_psf_oversampled);
      psfOversamplingScale = options->psfOversamplingScales[nn];
    }
    psfOversampleInfo->AddPsfPixels(psfOversampledPixels, nColumns_psf_oversampled,
   									nRows_psf_oversampled, newPsfOversampledPixelsFlag);
    psfOversampleInfo->AddOversamplingScale(psfOversamplingScale);
    psfOversamplingInfoVect.push_back(psfOversampleInfo);
  }
  
  return 0;
}
