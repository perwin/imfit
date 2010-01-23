/* FILE: readimage.cpp ------------------------------------------------- */
/* VERSION 0.10
 *
 *   Function to read in a FITS image and store it in a 1-D array, using
 * cfitsio routines.
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
 *     [version 0.10:] 11 Nov 2009: Created.
 */


/* ------------------------ Include Files (Header Files )--------------- */

#include <stdlib.h>
#include <string>
#include "fitsio.h"
#include "readimage.h"


/* ---------------- Definitions ---------------------------------------- */



/* ------------------- Function Prototypes ----------------------------- */
static void PrintError( int status );


/* ------------------------ Global Variables --------------------------- */


/* ------------------------ Module Variables --------------------------- */



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
 */
double * ReadImageAsVector( std::string filename, long *nColumns, long *nRows )
{
  fitsfile  *imfile_ptr;
  double  *imageVector;
  int  status, nfound;
  int  problems;
  long  naxes[2], nPixelsTot;
  long  firstPixel[2] = {1, 1};
  long  n_rows, n_columns;
  
  status = problems = 0;
  
   /* Open the FITS file: */
  problems = fits_open_file(&imfile_ptr, filename.c_str(), READONLY, &status);
  if ( problems )
    PrintError(status);

  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  problems = fits_read_keys_lng(imfile_ptr, "NAXIS", 1, 2, naxes, &nfound,
				  &status);
  if ( problems )
    PrintError(status);
  printf("Image keywords: NAXIS1 = %ld, NAXIS2 = %ld\n", naxes[0], naxes[1]);

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
  if ( problems )
    PrintError(status);   // Calling PrintError() will exit program...

  printf("\nImage read.\n");

  problems = fits_close_file(imfile_ptr, &status);
  if ( problems )
    PrintError(status);

  return imageVector;
}




/* ---------------- FUNCTION: PrintError --------------------------- */

static void PrintError( int status )
{

  if ( status ) {
    fits_report_error(stderr, status);
    exit(status);
  }
}



/* END OF FILE: readimage.cpp ------------------------------------------ */
