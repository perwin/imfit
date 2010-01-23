/*   Public interfaces for FITS image wrapper routines.
 *   INTERFACE VERSION 0.2
 */

#ifndef _IMAGE_IO_H
#define _IMAGE_IO_H

#include <string>


double * ReadImageAsVector( std::string filename, int *nColumns, int *nRows,
														bool verbose=false );

void SaveVectorAsImage( double *pixelVector, std::string filename, int nColumns,
                         int nRows );

#endif  // _IMAGE_IO_H
