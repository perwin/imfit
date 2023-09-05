// Copyright 2020--2023 by Peter Erwin.
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// The following requires GSL version 2.0 or later
#include "gsl/gsl_spline2d.h"
#include "psf_interpolators.h"



/* ---------------- Definitions ---------------------------------------- */

const double PI = 3.14159265358979;
const double PI_SQUARED = 9.86960440108936;



// DERIVED CLASS: PsfInterpolator_bicubic -- uses GNU Scientific Library's
// 2D bicubic interpolation

// Internally, we work in PSF-center-relative coordinates, which run from
// (xMin,yMin) = (-halfXsize,-halfYsize) to (xMax,yMax) = (+halfXsize,+halfYsize)
// and the center of the PSF is at (0,0)

/* ---------------- CONSTRUCTOR ---------------------------------------- */

PsfInterpolator_bicubic::PsfInterpolator_bicubic( double *inputImage, int nCols_image, 
													int nRows_image )
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
  
  interpolatorType = kInterpolator_bicubic;
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
// This function calculates and returns the value of the bicubic spline 
// interpolation kernel at x_diff,y_diff, with those coordinates being 
// relative to the center of the PSF. The corresponding calculations and
// call in PointSource::GetValue are
//    x_diff = x - x0;
//    y_diff = y - y0;
//    normalizedIntensity = psfInterpolator->GetValue(x_diff, y_diff);


double PsfInterpolator_bicubic::GetValue( double x, double y )
{
  double newVal;
  if ((x < deltaXMin) || (x > deltaXMax) || (y < deltaYMin) || (y > deltaYMax))
    newVal = 0.0;
  else
    newVal = gsl_spline2d_eval(splineInterp, x, y, xacc, yacc);
  return newVal;
}



// DERIVED CLASS: PsfInterpolator_lanczos2 -- uses Lanczos2 interpolation

/* ---------------- CONSTRUCTOR ---------------------------------------- */

PsfInterpolator_lanczos2::PsfInterpolator_lanczos2( double *inputImage, int nCols_image, 
													int nRows_image )
{
  nColumns = nCols_image;
  nRows = nRows_image;
  nPixelsTot = (long)(nColumns * nRows);
  psfDataArray = inputImage;
  
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

  interpolatorType = kInterpolator_lanczos2;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

PsfInterpolator_lanczos2::~PsfInterpolator_lanczos2( )
{
  free(xArray);
  free(yArray);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the value of the Lanczos2
// interpolation kernel, convolved with the PSF image, at x_diff,y_diff, with 
// those coordinates being relative to the center of the PSF. The corresponding 
// calculations and call in PointSource::GetValue are
//    x_diff = x - x0;
//    y_diff = y - y0;
//    normalizedIntensity = psfInterpolator->GetValue(x_diff, y_diff);

double PsfInterpolator_lanczos2::GetValue( double x, double y )
{
  double newVal;
  double lanczosScaling;
  int i_data_mid_x, i_data_mid_y;
  int i_data_x, i_data_y;
  double x_dat, y_dat;
  
  if ((x < deltaXMin) || (x > deltaXMax) || (y < deltaYMin) || (y > deltaYMax))
    newVal = 0.0;
  else {
    // FIXME: do Lanczos2 interpolation
    // reminder: xArray runs from [-halfXwidth, .., +halfXwidth], etc.
    i_data_mid_x = FindIndex(xArray, x);
    i_data_mid_y = FindIndex(yArray, y);
    
    newVal = 0.0;
    // loop over columns in PSF image
    for (int i = -2; i <= 2; i++) {
      i_data_x = i_data_mid_x + i;
      if ((i_data_x < 0) || (i_data_x >= nColumns))
        newVal += 0.0;  // outside PSF image (in x)
      else {
        x_dat = xArray[i_data_x];  // current x-value for PSF pixels
        // loop over rows in PSF image (for this column)
        for (int j = -2; j <= 2; j++) {
          i_data_y = i_data_mid_y + j;
          if ((i_data_y < 0) || (i_data_y >= nRows))
            newVal += 0.0;  // outside PSF image (in y)
          else {
            y_dat = yArray[i_data_y];  // current y-value for PSF pixels
            lanczosScaling = Lanczos(x - x_dat, 2) * Lanczos(y - y_dat, 2);
            newVal += lanczosScaling * psfDataArray[i_data_x*nColumns + i_data_y];
          }
        }
      }
    }
  }
  return newVal;
}



// DERIVED CLASS: PsfInterpolator_lanczos3 -- uses Lanczos3 interpolation

/* ---------------- CONSTRUCTOR ---------------------------------------- */

PsfInterpolator_lanczos3::PsfInterpolator_lanczos3( double *inputImage, int nCols_image, 
													int nRows_image )
{
  nColumns = nCols_image;
  nRows = nRows_image;
  nPixelsTot = (long)(nColumns * nRows);
  psfDataArray = inputImage;
  
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

  interpolatorType = kInterpolator_lanczos3;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

PsfInterpolator_lanczos3::~PsfInterpolator_lanczos3( )
{
  free(xArray);
  free(yArray);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the value of the Lanczos3
// interpolation kernel, convolved with the PSF image, at x_diff,y_diff, with 
// those coordinates being relative to the center of the PSF. The corresponding 
// calculations and call in PointSource::GetValue are
//    x_diff = x - x0;
//    y_diff = y - y0;
//    normalizedIntensity = psfInterpolator->GetValue(x_diff, y_diff);

double PsfInterpolator_lanczos3::GetValue( double x, double y )
{
  double newVal;
  double lanczosScaling;
  int i_data_mid_x, i_data_mid_y;
  int i_data_x, i_data_y;
  double x_dat, y_dat;
  
  if ((x < deltaXMin) || (x > deltaXMax) || (y < deltaYMin) || (y > deltaYMax))
    newVal = 0.0;
  else {
    // FIXME: do Lanczos3 interpolation
    // reminder: xArray runs from [-halfXwidth, .., +halfXwidth], etc.
    i_data_mid_x = FindIndex(xArray, x);
    i_data_mid_y = FindIndex(yArray, y);
    
    newVal = 0.0;
    // loop over columns in PSF image
    for (int i = -3; i <= 3; i++) {
      i_data_x = i_data_mid_x + i;
      if ((i_data_x < 0) || (i_data_x >= nColumns))
        newVal += 0.0;  // outside PSF image (in x)
      else {
        x_dat = xArray[i_data_x];  // current x-value for PSF pixels
        // loop over rows in PSF image (for this column)
        for (int j = -3; j <= 3; j++) {
          i_data_y = i_data_mid_y + j;
          if ((i_data_y < 0) || (i_data_y >= nRows))
            newVal += 0.0;  // outside PSF image (in y)
          else {
            y_dat = yArray[i_data_y];  // current y-value for PSF pixels
            lanczosScaling = Lanczos(x - x_dat, 3) * Lanczos(y - y_dat, 3);
            newVal += lanczosScaling * psfDataArray[i_data_x*nColumns + i_data_y];
          }
        }
      }
    }
  }
  return newVal;
}



// Extra non-method functions

// Find the index i into monotonically inreasing, evenly spaced array
// xArray for which xArray[i] is the largest value < xVal.
int FindIndex( double xArray[], double xVal )
{
  double deltaX, deltaX_relative;
  deltaX = xArray[1] - xArray[0];  // spacing between elements of xArray
  deltaX_relative = xVal - xArray[0];  // distance of xVal from start value of xArray
  return (int)(floor(deltaX_relative / deltaX));
}


// Generalized Lanczsos function (a.k.a. Lanczos-windowed sinc function)
double Lanczos( double x, int n )
{
  double x_abs = abs(x);
  if (x_abs < 1.0e-6)
    return 1.0;
  else if (x_abs > n)
    return 0.0;
  else
    return (n * sin(PI*x_abs) * sin(PI*x_abs/n)) / (PI_SQUARED*x_abs*x_abs);
}
