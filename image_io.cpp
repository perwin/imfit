/* FILE: image_io.cpp -------------------------------------------------- */
/* VERSION 0.3
 *
 *   Function for dealing with FITS files, using cfitsio routines:
 *   1. Read in a FITS image and store it in a 1-D array
 *   2. Given a 1-D array (and # rows, columns specification), save it as
 *      a FITS image.
 *
 *   Based on fitsimage_readwrite.cpp.
 * 
 * The proper translations are:
 * NAXIS1 = naxes[0] = nColumns = sizeX;
 * NAXIS2 = naxes[1] = nRows = sizeY.
 *
 *   Must be linked with the cfitsio library.
 *
 *   MODIFICATION HISTORY:
 *     [version 0.2:] 20 Aug 2010: Added writing of (optional) comments to
 # output FITS header.
 *     [version 0.15:] 27 Mar 2010: Added writing of DATE header and saving
 * of image in single-precision format to SaveVectorAsImage().
 *     [version 0.10:] 17 Nov 2009: Created by extending readimage.cpp to
 * include SaveVectorAsImage().
 */

// Copyright 2010, 2011, 2012, 2014 by Peter Erwin.
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



/* ------------------------ Include Files (Header Files )--------------- */

#include <stdlib.h>
#include <string>
#include <vector>

#include "fitsio.h"

#include "image_io.h"


/* ---------------- Definitions ---------------------------------------- */



/* ------------------- Function Prototypes ----------------------------- */
static void PrintError( int status );


/* ------------------------ Global Variables --------------------------- */


/* ------------------------ Module Variables --------------------------- */



/* ---------------- FUNCTION: GetImageSize ----------------------------- */
/*    Given a filename, it opens the file, reads the size of the image and
 * stores that size in *nRows and *nColumns.
 *
 *    Returns 0 for successful operation, -1 if a CFITSIO-related error occurred.
 */
int GetImageSize( std::string filename, int *nColumns, int *nRows, bool verbose )
{
  fitsfile  *imfile_ptr;
  int  status, nfound;
  int  problems;
  long  naxes[2];
  int  n_rows, n_columns;

  status = problems = 0;
  
   /* Open the FITS file: */
  problems = fits_open_file(&imfile_ptr, filename.c_str(), READONLY, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems opening FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return -1;
  }

  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  problems = fits_read_keys_lng(imfile_ptr, "NAXIS", 1, 2, naxes, &nfound,
				  &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems reading FITS keywords from file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return -1;
  }
  if (verbose)
    printf("GetImageSize: Image keywords: NAXIS1 = %ld, NAXIS2 = %ld\n", naxes[0], naxes[1]);

  n_columns = naxes[0];      // FITS keyword NAXIS1 = # columns
  *nColumns = n_columns;
  n_rows = naxes[1];         // FITS keyword NAXIS2 = # rows
  *nRows = n_rows;

  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems closing FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return -1;
  }
  
  return 0;
}



/* ---------------- FUNCTION: ReadImageAsVector ------------------------ */
/*    Given a filename, it opens the file, reads the size of the image and
 * stores that size in *nRows and *nColumns, then allocates memory for a 1-D
 * array to hold the image and reads the image from the file into the
 * array.  Finally, it returns the image array -- or, more precisely, it
 * returns a pointer to the array; it also stores the image dimensions
 * in the pointer-parameters nRows and nColumns.
 *    Note that this function does *not* use Numerical Recipes functions; instead
 * it allocates a standard 1-D C vector [this means that the first index will
 * be 0, not 1].
 *
 *    Returns 0 for successful operation, -1 if a CFITSIO-related error occurred.
 *
 */
double * ReadImageAsVector( std::string filename, int *nColumns, int *nRows,
														bool verbose )
{
  fitsfile  *imfile_ptr;
  double  *imageVector;
  int  status, nfound;
  int  problems;
  long  naxes[2];
  int  nPixelsTot;
  long  firstPixel[2] = {1, 1};
  int  n_rows, n_columns;
  
  status = problems = 0;
  
   /* Open the FITS file: */
  problems = fits_open_file(&imfile_ptr, filename.c_str(), READONLY, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems opening FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return NULL;
  }

  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  problems = fits_read_keys_lng(imfile_ptr, "NAXIS", 1, 2, naxes, &nfound,
				  &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems reading FITS keywords from file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return NULL;
  }
  if (verbose)
    printf("ReadImageAsVector: Image keywords: NAXIS1 = %ld, NAXIS2 = %ld\n", naxes[0], naxes[1]);

  n_columns = naxes[0];      // FITS keyword NAXIS1 = # columns
  *nColumns = n_columns;
  n_rows = naxes[1];         // FITS keyword NAXIS2 = # rows
  *nRows = n_rows;
  nPixelsTot = n_columns * n_rows;      // number of pixels in the image
  
  // Allocate memory for the image-data vector:
  imageVector = (double *) malloc(nPixelsTot * sizeof(double));
  // Read in the image data
  problems = fits_read_pix(imfile_ptr, TDOUBLE, firstPixel, nPixelsTot, NULL, imageVector,
                            NULL, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems reading pixel data from FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return NULL;
  }

  if (verbose)
    printf("\nReadImageAsVector: Image read.\n");

  problems = fits_close_file(imfile_ptr, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems closing FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return NULL;
  }

  return imageVector;
}



/* ---------------- FUNCTION: SaveVectorAsImage ------------------------ */
/*
 *    Returns 0 for successful operation, -1 if a CFITSIO-related error occurred.
 */
int SaveVectorAsImage( double *pixelVector, std::string filename, int nColumns,
                         int nRows, std::vector<std::string> comments )
{
  fitsfile  *imfile_ptr;
  std::string  finalFilename = "!";   // starting filename with "!" ==> clobber any existing file
  int  status, problems;
  long  naxes[2];
  int  nPixels;
  long  firstPixel[2] = {1, 1};

  status = problems = 0;
  
  naxes[0] = nColumns;
  naxes[1] = nRows;
  nPixels = nColumns * nRows;
  
  /* Create the FITS file: */
  //    NOTE: need to prefix filename with "!" if we want to clobber existing file...
  finalFilename += filename;
  fits_create_file(&imfile_ptr, finalFilename.c_str(), &status);
  /* Create the primary image (single-precision floating-point format) */
  fits_create_img(imfile_ptr, FLOAT_IMG, 2, naxes, &status);
  
  // Insert keyword writing here ...
  if (comments.size() > 0) {
    for (int i = 0; i < comments.size(); i++)
      fits_write_comment(imfile_ptr, comments[i].c_str(), &status);
  }
  fits_write_date(imfile_ptr, &status);

  /* Write vector of pixel values to the image (note that cfitsio automatically handles
   * the conversion from double-precision (pixelVector values) to single-precision
   * output image format) */
  problems = fits_write_pix(imfile_ptr, TDOUBLE, firstPixel, nPixels, pixelVector,
                            &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems writing pixel data to FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return -1;
  }

  problems = fits_close_file(imfile_ptr, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems closing FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return -1;
  }
  
  return 0;
}




/* ---------------- FUNCTION: PrintError --------------------------- */

static void PrintError( int status )
{

  if ( status ) {
    fits_report_error(stderr, status);
    fprintf(stderr, "\n");
//    exit(status);
  }
}



/* END OF FILE: readimage.cpp ------------------------------------------ */
