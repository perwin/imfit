/** @file
    \brief Public interfaces for the FITS image I/O routines.
 *
 */

#ifndef _IMAGE_IO_H
#define _IMAGE_IO_H

#include <string>
#include <vector>


/// Gets dimensions (nColumns, nRows) of specified FITS image
int GetImageSize( const std::string filename, int *nColumns, int *nRows, const bool verbose=false );

/// \brief Reads image data from specified FITS image, returning it as 1D array
///        (with image dimensions stored in nColumns, nRows)
double * ReadImageAsVector( const std::string filename, int *nColumns, int *nRows,
							const bool verbose=false );

/// \brief Saves image data (1D array, logical dimensions nColumns x nRows) as
///        FITS file, with comments added to header.
int SaveVectorAsImage( double *pixelVector, const std::string filename, const int nColumns,
                         const int nRows, std::vector<std::string> comments );

#endif  // _IMAGE_IO_H
