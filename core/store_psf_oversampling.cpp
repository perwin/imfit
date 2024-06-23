// Code for consolidating the processing and aggregating oversampled-PSF info and
// data, starting from ImageInfo objects -- meant to be called from 
// makemultimages/multimfit/etc. main().

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

// Inputs:
// 	imageInfoVect[i] = single ImageInfo instance
// 
// Outputs:
// 	status indicating success or failure
// 	Add new PsfOversamplingInfo instance to psfOversamplingInfoVect
// 
// 
// CALL AS:
//   status = ExtractAndStorePsfOversampling(imageInfoVect[i], i, psfOversamplingInfoVect);


#include <stdio.h>
#include <vector>

using namespace std;

#include "image_io.h"
#include "definitions_multimage.h"
#include "psf_oversampling_info.h"


int ExtractAndStorePsfOversampling( const ImageInfo &imageInfo, int i,
									vector<PsfOversamplingInfo *>  &oversamplingInfoVect )
{
  int  nColumns_psf_oversampled, nRows_psf_oversampled;
  long  nPixels_psf_oversampled;
  double *psfOversampledPixels;
  PsfOversamplingInfo  *newPsfOversamplingInfo;
  int  psfOversamplingScale;
  int  nOversampledPsfImages = imageInfo.nOversampledFileNames;
  int  nOversampledScales = imageInfo.nOversamplingScales;
  int  nOversampledRegions = imageInfo.nOversampleRegions;
  
  if (nOversampledPsfImages != nOversampledScales) {
    fprintf(stderr, "\n*** ERROR: image #%d: number of oversampling scales (%d) is not the same\n", 
  					i + 1, nOversampledScales);
    fprintf(stderr, "           as number of oversampled-PSF images (%d)!\n\n",
  					nOversampledPsfImages);
    return -1;
  }
  if ((nOversampledPsfImages > 1) && (nOversampledPsfImages != nOversampledRegions)) {
    fprintf(stderr, "\n*** ERROR: image #%d: number of oversampled-PSF images (%d) must be = 1 OR\n", 
  					i + 1, nOversampledPsfImages);
    fprintf(stderr, "           must be same as number of oversampled-PSF regions (%d)!\n\n",
  					nOversampledRegions);
    return -1;
  }
  if (! imageInfo.oversampleRegionSet) {
    fprintf(stderr, "\n*** ERROR: image #%d: no oversampling region within the main image was defined!\n",
    				i + 1);
    fprintf(stderr, "           (use \"OVERSAMPLED_REGION x1:x2,y1:y2\" to specify the region)\n\n");
    return -1;
  }

  // Define and populate new instances of PsfOversamplingInfo (one per user-specified
  // oversampling region), then add them to oversamplingInfoVect
  for (int nn = 0; nn < nOversampledRegions; nn++) {
    newPsfOversamplingInfo = new PsfOversamplingInfo();
    newPsfOversamplingInfo->AddRegionString(imageInfo.psfOversampleRegions[nn]);
    bool newPsfOversampledPixelsFlag = false;
    if ( (nn == 0) || ((nn > 0) && (nOversampledPsfImages > 1)) ) {
      // Always read PSF image and get oversampling scale from options object the
      // first time through; do it again if user supplied more than one image and
      // scale (otherwise, we reuse the same image and scale for subsequent regions)
      printf("Reading oversampled PSF image (\"%s\") ...\n", imageInfo.psfOversampledFileNames[nn].c_str());
      psfOversampledPixels = ReadImageAsVector(imageInfo.psfOversampledFileNames[nn], 
	    							&nColumns_psf_oversampled, &nRows_psf_oversampled);
      if (psfOversampledPixels == nullptr) {
        fprintf(stderr, "\n*** ERROR: image #%d: Unable to read oversampled PSF image file \"%s\"!\n\n", 
		      			i + 1, imageInfo.psfOversampledFileNames[nn].c_str());
        return -1;
      }
      newPsfOversampledPixelsFlag = true;
      nPixels_psf_oversampled = (long)nColumns_psf_oversampled * (long)nRows_psf_oversampled;
      printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
           nColumns_psf_oversampled, nRows_psf_oversampled, nPixels_psf_oversampled);
      psfOversamplingScale = imageInfo.psfOversamplingScales[nn];
    }
    newPsfOversamplingInfo->AddPsfPixels(psfOversampledPixels, nColumns_psf_oversampled,
  									nRows_psf_oversampled, newPsfOversampledPixelsFlag);
    newPsfOversamplingInfo->AddOversamplingScale(psfOversamplingScale);
    oversamplingInfoVect.push_back(newPsfOversamplingInfo);
  }
  
  return 0;
}
