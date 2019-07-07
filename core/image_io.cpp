/* FILE: image_io.cpp -------------------------------------------------- */
/*
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
 * The public filename interface for specifying HDUs/extensions is 0-based:
 *    filename.fits[0] = primary HDU
 *    filename.fits[1] = second HDU (first extension)
 *    etc.
 * HOWEVER, the *internal* number scheme for HDU numbers is 1-based:
 *    1 = primary HDU
 *    2 = second HDU (first extension)
 *    etc.
 *
 *
 *   Must be linked with the cfitsio library.
 *
 * NOTES AND WARNINGS ABOUT CFITSIO:
 *   CFITSIO functions typically take as input a pointer-to-int called "status".
 * Although this is generally meant to be an *output* parameter (because CFITSIO
 * is based on old Fortran code, where subroutines can't return values), quite a
 * few internal functions will actually check the value of (*status) and, if
 * it's nonzero, will then exit as though an error had been encountered (or
 * if it's < 0, will behave in a different fashion; there are about eight or so
 * negative "error codes" defined so as to control internal behavior).
 *
 *   Confusingly, the value of *status is usually *also* returned as the return
 * value of the function. Some CFITSIO functions have code like this:
 *       *status = somethingorother;
 *       return (*status);
 * or even:
 *       *status = fits_close_file(tmpfptr,status);
 *
 *   Quoting from the "Examples Programs" part of the documentation:
 *
 *   "Almost every CFITSIO routine has a status parameter as the last
 *   argument. The status value is also usually returned as the value of the
 *   function itself. Normally status = 0, and a positive status value
 *   indicates an error of some sort. The status variable must always be
 *   initialized to zero before use, because if status is greater than zero
 *   on input then the CFITSIO routines will simply return without doing
 *   anything. This `inherited status' feature, where each CFITSIO routine
 *   inherits the status from the previous routine, makes it unnecessary to
 *   check the status value after every single CFITSIO routine call.
 *   Generally you should check the status after an especially important or
 *   complicated routine has been called, or after a block of closely related
 *   CFITSIO calls. This example program has taken this feature to the
 *   extreme and only checks the status value at the very end of the program."
 *
 *   [Comment: nasty, because the `inherited status' feature amounts to hidden
 *   state...]
 *
 *   MODIFICATION HISTORY:
 *     [version 0.2:] 20 Aug 2010: Added writing of (optional) comments to
 # output FITS header.
 *     [version 0.15:] 27 Mar 2010: Added writing of DATE header and saving
 * of image in single-precision format to SaveVectorAsImage().
 *     [version 0.10:] 17 Nov 2009: Created by extending readimage.cpp to
 * include SaveVectorAsImage().
 */

// Copyright 2010--2018 by Peter Erwin.
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
#include <tuple>
#include "fftw3.h"
#include "fitsio.h"

#include "image_io.h"


/* ---------------- Definitions ---------------------------------------- */



/* ------------------- Function Prototypes ----------------------------- */
static void PrintError( int status );




// Possible FITS file situations:
// single HDU is valid image -- OK
// single HDU is not image -- REJECT
// user-specified HDU
//    check if that is valid or not

/* ---------------- FUNCTION: CheckHDUForImage ------------------------- */
bool CheckHDUForImage( fitsfile *imageFile_ptr, int hduNum, int *status_ptr )
{
  bool validImageFlag = false;
  int  hduType, naxis;

  fits_movabs_hdu(imageFile_ptr, hduNum, &hduType, status_ptr);
  if (hduType == IMAGE_HDU) {   /* primary array or image HDU */
    fits_get_img_dim(imageFile_ptr, &naxis, status_ptr);
    if (naxis == 2)
      validImageFlag = true;
  }
  return validImageFlag;
}


/* ---------------- FUNCTION: CheckForImage ---------------------------- */
///    Given the filename of a FITS image, this function opens the file and
/// checks the first header-data unit to see if it is a 2D image. If it is
/// *not* and a second HDU exists, then that is also checked.
///
///   Returns -1 if something goes wrong opening the file or if primary HDU is not
/// a proper 2D image; otherwise, returns 1 if the specified FITS file (including 
/// implicitly specified extension number, if that was part of filename) is a proper2D image.
///
///   Currently uses int as return value for possible case of finding and returning
/// a different HDU. FIXME: probably good idea to make this bool return value.
int CheckForImage( const std::string filename, const bool verbose )
{
  fitsfile  *imfile_ptr;
  int  problems = 0;
  int  status = 0;
  int  validHDU = -1;   // meaningless bad value (for possible test purposes)
  int  currentHDU = -1;

  problems = fits_open_file(&imfile_ptr, filename.c_str(), READONLY, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems opening FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return -1;
  }

  // Get the current HDU (= 1 = primaryHDU if use didn't specify anything, = n+1 if
  // user specified "filename.fits[n]"
  fits_get_hdu_num(imfile_ptr, &currentHDU);
  
  if (CheckHDUForImage(imfile_ptr, currentHDU, &status))
    validHDU = 1;
  
  fits_close_file(imfile_ptr, &status);
  if (validHDU < 0)
    fprintf(stderr, "\n*** WARNING: Unable to find valid 2D image in extension %d (HDU %d) of %s!\n", 
    		currentHDU - 1, currentHDU, filename.c_str());
  return validHDU;
}
  
  
/* ---------------- FUNCTION: GetImageSize ----------------------------- */
///    Given the filename of a FITS image, this function opens the file, reads the 
/// size of the image and returns the dimensions in nRows and nColumns.
///
///   Returns 0 for successful operation, -1 if a CFITSIO-related error occurred.
std::tuple<int, int, int> GetImageSize( const std::string filename, const bool verbose )
{
  fitsfile  *imfile_ptr;
  int  status = 0;
  int  problems = 0;
  int  nfound;
  long  naxes[2];
  int  n_columns = 0;
  int  n_rows = 0;
  
  problems = fits_open_file(&imfile_ptr, filename.c_str(), READONLY, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems opening FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return std::make_tuple(n_columns, n_rows, -1);
  }

  // Check to make sure primary HDU is an image, and is actually a 2D array
  int  currentHDU = -1;
  fits_get_hdu_num(imfile_ptr, &currentHDU);
  if (! CheckHDUForImage(imfile_ptr, currentHDU, &status)) {
    fprintf(stderr, "\n*** WARNING: Unable to find valid 2D image in extension %d (HDU %d) of %s!\n", 
    		currentHDU - 1, currentHDU, filename.c_str());
    return std::make_tuple(n_columns, n_rows, -1);
  }  

  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  problems = fits_read_keys_lng(imfile_ptr, "NAXIS", 1, 2, naxes, &nfound,
				  &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems reading FITS keywords from file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return std::make_tuple(n_columns, n_rows, -1);
  }
  if (verbose)
    printf("GetImageSize: Image keywords: NAXIS1 = %ld, NAXIS2 = %ld\n", naxes[0], naxes[1]);

  problems = fits_close_file(imfile_ptr, &status);
  if ( problems ) {
    fprintf(stderr, "\n*** WARNING: Problems closing FITS file \"%s\"!\n    FITSIO error messages follow:", filename.c_str());
    PrintError(status);
    return std::make_tuple(n_columns, n_rows, -1);
  }
  
  n_columns = naxes[0];      // FITS keyword NAXIS1 = # columns
  n_rows = naxes[1];         // FITS keyword NAXIS2 = # rows

  return std::make_tuple(n_columns, n_rows, 0);
}



/* ---------------- FUNCTION: ReadImageAsVector ------------------------ */
///    Given a filename, it opens the file, reads the size of the image and
/// stores that size in nRows and nColumns, then allocates memory for a 1-D
/// array to hold the image and reads the image from the file into the
/// array.  Finally, it returns the image array -- or, more precisely, it
/// returns a pointer to the array; it also stores the image dimensions
/// in the pointer-parameters nRows and nColumns.
///
///    Returns NULL (and prints error message) if a CFITSIO-related error occurred.
double * ReadImageAsVector( const std::string filename, int *nColumns, int *nRows,
							const bool verbose )
{
  fitsfile  *imfile_ptr;
  double  *imageVector;
  int  status = 0;
  int  problems = 0;
  int  validHDU_flag = 0;
  int  nfound;
  long  naxes[2];
  long  nPixelsTot;
  long  firstPixel[2] = {1, 1};
  int  n_rows, n_columns;
  
  status = problems = 0;
  
  // Check to make sure this is a valid image, then open the FITS file for 
  // further operations
  validHDU_flag = CheckForImage(filename);
  if (validHDU_flag > 0)
    fits_open_file(&imfile_ptr, filename.c_str(), READONLY, &status);
  else
    return NULL;

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
  nPixelsTot = (long)n_columns * (long)n_rows;      // number of pixels in the image
  
  // Read in the image data
  imageVector = fftw_alloc_real(nPixelsTot);
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
    fftw_free(imageVector);
    return NULL;
  }

  return imageVector;
}



/* ---------------- FUNCTION: SaveVectorAsImage ------------------------ */
///    Saves specified 1D array (representing an image with size nColumns x nRows)
/// in specified filename as a FITS image, with strings stored in comments
/// vector written as comments to the FITS header.
///
///    Returns 0 for successful operation, -1 if a CFITSIO-related error occurred.
int SaveVectorAsImage( double *pixelVector, const std::string filename, int nColumns,
                         int nRows, std::vector<std::string> comments )
{
  fitsfile  *imfile_ptr;
  std::string  finalFilename = "!";   // starting filename with "!" ==> clobber any existing file
  int  status = 0;
  int  problems = 0;
  long  naxes[2];
  long  nPixels;
  long  firstPixel[2] = {1, 1};

  status = problems = 0;
  
  // Check for bad input
  if (pixelVector == NULL) {
    fprintf(stderr, "\n*** WARNING: input image array to SaveVectorAsImage is NULL!\n");
    return -2;
  }
  
  naxes[0] = nColumns;
  naxes[1] = nRows;
  nPixels = (long)nColumns * (long)nRows;
  
  /* Create the FITS file: */
  //    NOTE: need to prefix filename with "!" if we want to clobber existing file...
  finalFilename += filename;
  fits_create_file(&imfile_ptr, finalFilename.c_str(), &status);
  /* Create the primary image (single-precision floating-point format) */
  fits_create_img(imfile_ptr, FLOAT_IMG, 2, naxes, &status);
  
  // Insert keyword writing here ...
  if (comments.size() > 0) {
    for (int i = 0; i < (int)comments.size(); i++)
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
  }
}



/* END OF FILE: image_io.cpp ------------------------------------------- */
