/*   Public interfaces for FITS image wrapper routines.
 *   INTERFACE VERSION 0.2
 */

#ifndef _IMAGE_IO_H
#define _IMAGE_IO_H

#include <string>
#include <vector>


int GetImageSize( std::string filename, int *nColumns, int *nRows, bool verbose=false );

double * ReadImageAsVector( std::string filename, int *nColumns, int *nRows,
														bool verbose=false );

int SaveVectorAsImage( double *pixelVector, std::string filename, int nColumns,
                         int nRows, std::vector<std::string> comments );

#endif  // _IMAGE_IO_H
