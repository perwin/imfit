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


// std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
// 										OptionsBase *options, bool &maskPixelsAllocated, 
// 										bool &errorPixelsAllocated )
// {
//   int  nMaskColumns, nMaskRows;
//   int  nErrColumns, nErrRows;
//   double *maskPixels = nullptr;
//   double *errorPixels = nullptr;
//   
//   maskPixelsAllocated = false;
//   errorPixelsAllocated = false;
//   
//   /* Get and check mask image */
//   if (options->maskImagePresent) {
//     printf("Reading mask image (\"%s\") ...\n", options->maskFileName.c_str());
//     maskPixels = ReadImageAsVector(options->maskFileName, &nMaskColumns, &nMaskRows);
//     if (maskPixels == NULL) {
//       fprintf(stderr,  "\n*** ERROR: Unable to read mask file \"%s\"!\n\n", 
//     			options->maskFileName.c_str());
//       return std::make_tuple(maskPixels, errorPixels, -1);
//     }
//     if ((nMaskColumns != nColumns) || (nMaskRows != nRows)) {
//       fprintf(stderr, "\n*** ERROR: Dimensions of mask image (%s: %d columns, %d rows)\n",
//              options->maskFileName.c_str(), nMaskColumns, nMaskRows);
//       fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
//              options->imageFileName.c_str(), nColumns, nRows);
//       return std::make_tuple(maskPixels, errorPixels, -1);
//     }
//     maskPixelsAllocated = true;
//   }
//            
//   /* Get and check error image, if supplied */
//   if (options->noiseImagePresent) {
//     printf("Reading noise image (\"%s\") ...\n", options->noiseFileName.c_str());
//     errorPixels = ReadImageAsVector(options->noiseFileName, &nErrColumns, &nErrRows);
//     if (errorPixels == NULL) {
//       fprintf(stderr,  "\n*** ERROR: Unable to read noise-image file \"%s\"!\n\n", 
//     			options->noiseFileName.c_str());
//       return std::make_tuple(maskPixels, errorPixels, -1);
//     }
//     if ((nErrColumns != nColumns) || (nErrRows != nRows)) {
//       fprintf(stderr, "\n*** ERROR: Dimensions of error image (%s: %d columns, %d rows)\n",
//              options->noiseFileName.c_str(), nErrColumns, nErrRows);
//       fprintf(stderr, "do not match dimensions of data image (%s: %d columns, %d rows)!\n\n",
//              options->imageFileName.c_str(), nColumns, nRows);
//       return std::make_tuple(maskPixels, errorPixels, -1);
//     }
//     errorPixelsAllocated = true;
//   }
//   
//   return std::make_tuple(maskPixels, errorPixels, 0);
// }
