/** @file
    \brief Public interfaces for the FITS image I/O routines.
 *
 */

#ifndef _IMAGE_IO_H
#define _IMAGE_IO_H

#include <string>
#include <vector>
#include "fitsio.h"


/// Checks to see if current HDU in FITS file has proper 2D image
bool CheckHDUForImage( fitsfile *imageFile_ptr, int hduNum, int *status_ptr );

/// Checks to see if FITS file has proper 2D image in first HDU (or second, if
/// first is empty)
int CheckForImage( const std::string filename, const bool verbose=false );

/// Gets dimensions (nColumns, nRows) of specified FITS image
std::tuple<int, int, int> GetImageSize( const std::string filename, const bool verbose=false );

/// \brief Reads image data from specified FITS image, returning it as 1D array
///        (with image dimensions stored in nColumns, nRows)
double * ReadImageAsVector( const std::string filename, int *nColumns, int *nRows,
							const bool verbose=false );

/// \brief Saves image data (1D array, logical dimensions nColumns x nRows) as
///        FITS file, with comments added to header.
int SaveVectorAsImage( double *pixelVector, const std::string filename, int nColumns,
                         int nRows, std::vector<std::string> comments );

int CountHeaderDataUnits( fitsfile  *imfile_ptr );

#endif  // _IMAGE_IO_H
