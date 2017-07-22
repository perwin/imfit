/** @file
    \brief Function for reading in (and checking) mask and/or error/noise
    image; for use with imfit and imfit-mcmc.
 *
 */
#ifndef _GETIMAGES_MASKERROR_H_
#define _GETIMAGES_MASKERROR_H_

#include <tuple>
#include <string>

#include "options_base.h"


std::tuple<double *, int> GetAndCheckImage( const string imageName, const string imageType,
											int nColumns_ref, int nRows_ref );

std::tuple<double *, double *, int> GetMaskAndErrorImages( int nColumns, int nRows, 
										OptionsBase *options,  bool &maskPixelsAllocated, 
										bool &errorPixelsAllocated );


#endif  // _GETIMAGES_MASKERROR_H_
