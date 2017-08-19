/** @file
    \brief Functions for reading in (and checking) various images (other than
    main data image): mask, error, PSF, and oversampled PSF images.
    For use with makeimage, imfit, and imfit-mcmc.
 *
 */
#ifndef _GETIMAGES_MASKERROR_H_
#define _GETIMAGES_MASKERROR_H_

#include <tuple>
#include <string>

#include "options_base.h"
#include "psf_oversampling_info.h"


/// Main utility function: reads in image from FITS file imageName and checks dimensions
/// against reference values (if latter are nonzero)
std::tuple<double *, int> GetAndCheckImage( const string imageName, const string imageType,
											int nColumns_ref, int nRows_ref );

/// Reads in mask and/or noise/error images
std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
										string &maskFileName, string &errorFileName, 
										bool &maskPixelsAllocated, bool &errorPixelsAllocated );
// std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
// 										const OptionsBase *options,  bool &maskPixelsAllocated, 
// 										bool &errorPixelsAllocated );

/// Reads in PSF image, returning dimensions as well
// std::tuple<double *, int, int, int> GetPsfImage( const OptionsBase *options );
std::tuple<double *, int, int, int> GetPsfImage( const string &psfFileName );

/// Reads in multiple oversampled PSF images, storing them (along with user-specified
/// info like oversampling scale and image regions for oversampling) in a vector
/// of (pointers to) PsfOversamplingInfo objects.
int GetOversampledPsfInfo( const OptionsBase *options, int xOffset, int yOffset, 
							vector<PsfOversamplingInfo *> &psfOversamplingInfoVect );


#endif  // _GETIMAGES_MASKERROR_H_
