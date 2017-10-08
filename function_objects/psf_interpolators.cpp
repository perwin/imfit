#include <stdlib.h>

#include "gsl/gsl_spline2d.h"
#include "psf_interpolators.h"



/* ---------------- Definitions ---------------------------------------- */



// DERIVED CLASS: PsfInterpolator_bicubic -- uses GNU Scientific Library's
// 2D bicubic interpolation

/* ---------------- CONSTRUCTOR ---------------------------------------- */

PsfInterpolator_bicubic::PsfInterpolator_bicubic( double *inputImage, int nCols_image, int nRows_image )
{
  nColumns = nCols_image;
  nRows = nRows_image;
  nPixelsTot = (long)(nColumns * nRows);
  xBound = (nColumns - 1) / 2.0;
  yBound = (nRows - 1) / 2.0;
  xArray = (double *)calloc((size_t)nColumns, sizeof(double));
  yArray = (double *)calloc((size_t)nRows, sizeof(double));
  for (int n = 0; n < nColumns; n++)
    xArray[n] = n - xBound;
  for (int n = 0; n < nRows; n++)
    yArray[n] = n - yBound;
  deltaXMin = -xBound;
  deltaXMax = xBound;
  deltaYMin = -yBound;
  deltaYMax = yBound;
  
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
  splineInterp = gsl_spline2d_alloc(gsl_interp2d_bicubic, nColumns, nRows);
  int result = gsl_spline2d_init(splineInterp, xArray, yArray, inputImage, nColumns, nRows);
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

PsfInterpolator_bicubic::~PsfInterpolator_bicubic( )
{
  gsl_spline2d_free(splineInterp);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  free(xArray);
  free(yArray);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double PsfInterpolator_bicubic::GetValue( double x, double y )
{
  double newVal;
  if ((x < deltaXMin) || (x > deltaXMax) || (y < deltaYMin) || (y > deltaYMax))
    newVal = 0.0;
  else
    newVal = gsl_spline2d_eval(splineInterp, x, y, xacc, yacc);
  return newVal;
}
