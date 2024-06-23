// Code for consolidating the processing and aggregating oversampled-PSF info and
// data, starting from ImageInfo objects -- meant to be called from 
// makemultimages/multimfit/etc. main().

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


#ifndef _STORE_PSF_OVERSAMPLING_H_
#define _STORE_PSF_OVERSAMPLING_H_

#include <vector>

using namespace std;

#include "definitions_multimage.h"
#include "psf_oversampling_info.h"


int ExtractAndStorePsfOversampling( const ImageInfo &imageInfoObj, int i,
									vector<PsfOversamplingInfo *>  &oversamplingInfoVect );

#endif /* _STORE_PSF_OVERSAMPLING_H_ */
