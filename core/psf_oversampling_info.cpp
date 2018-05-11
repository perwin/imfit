/* FILE: psf_oversampling_info.cpp ------------------------------------- */
/* 
 * Class definition for PsfOversamplingInfo.
 */

// Copyright 2017-2018 by Peter Erwin.
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

#include <string>

using namespace std;

#include "fftw3.h"  // so we can call fftw_free()
#include "psf_oversampling_info.h"
#include "utilities_pub.h"


/* ---------------- Definitions ---------------------------------------- */



/* ---------------- CONSTRUCTOR (default) ------------------------------ */
/// Default constructor for PsfOversamplingInfo class
PsfOversamplingInfo::PsfOversamplingInfo( )
{
  psfPixels = NULL;
  nColumns_psf = nRows_psf = 0;
  oversamplingScale = 1;
  regionString = "";
  X0_offset = Y0_offset = 0;
  pixelsArrayIsUnique = true;
  normalizePSF = true;
}


/* ---------------- CONSTRUCTOR ---------------------------------------- */
/// Constructor (with inputs) for PsfOversamplingInfo class
PsfOversamplingInfo::PsfOversamplingInfo( double *inputPixels, int nCols, int nRows, 
										int scale, string inputRegionString,
										int xOffset, int yOffset, bool isUnique,
										bool normalize )
{
  psfPixels = inputPixels;
  nColumns_psf = nCols;
  nRows_psf = nRows;
  oversamplingScale = scale;
  regionString = inputRegionString;
  X0_offset = xOffset;
  Y0_offset = yOffset;
  pixelsArrayIsUnique = isUnique;
  normalizePSF = normalize;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */
// Free the pixel-data array if we're an instance pointing to a unique
// array (or if we're the first instance pointing to a shared array)

/// Destructor for PsfOversamplingInfo class
PsfOversamplingInfo::~PsfOversamplingInfo( )
{
  if (pixelsArrayIsUnique)
    fftw_free(psfPixels);
}



/* ---------------- AddRegionString ------------------------------------ */
/// Add a region string (e.g., "[500:600,1080:1422]") defining a subsection
/// of an image that will be oversampled
void PsfOversamplingInfo::AddRegionString( string inputRegionString )
{
  regionString = inputRegionString;
}


/* ---------------- AddPsfPixels --------------------------------------- */
/// Add data describing a PSF image
void PsfOversamplingInfo::AddPsfPixels( double *inputPixels, int nCols, int nRows,
										bool isUnique )
{
  psfPixels = inputPixels;
  nColumns_psf = nCols;
  nRows_psf = nRows;
  pixelsArrayIsUnique = isUnique;
}


/* ---------------- AddOversamplingScale ------------------------------- */
/// Specify the pixel oversampling scale for an oversampling region (must
/// be an integer >= 1)
void PsfOversamplingInfo::AddOversamplingScale( int scale )
{
  oversamplingScale = scale;
}


/* ---------------- AddImageOffset ------------------------------------- */
/// Specify the offset of an overampling region (the location of its lower-left
/// corner within the full image)
void PsfOversamplingInfo::AddImageOffset( int X0, int Y0 )
{
  X0_offset = X0;
  Y0_offset = Y0;
}


/* ---------------- SetNormalizationFlag ------------------------------- */
/// Whether PSF image should be normalized or not
void PsfOversamplingInfo::SetNormalizationFlag( bool normalize )
{
  normalizePSF = normalize;
}


/* ---------------- GetNColumns ---------------------------------------- */
/// Returns x-size of PSF image
int PsfOversamplingInfo::GetNColumns( )
{
  return nColumns_psf;
}


/* ---------------- GetNRows ------------------------------------------- */
/// Returns y-size of PSF image
int PsfOversamplingInfo::GetNRows( )
{
  return nRows_psf;
}


/* ---------------- pixelsArrayIsUnique -------------------------------- */
/// Returns true if PSF image is unique to this oversampling region; returns
/// false if PSF image is shared with other oversampling regions
bool PsfOversamplingInfo::PixelsArrayIsUnique( )
{
  return pixelsArrayIsUnique;
}

// 
//     bool pixelsArrayIsUnique( );

/* ---------------- GetPsfPixels --------------------------------------- */
/// Returns a pointer to the PSF image pixel array
double * PsfOversamplingInfo::GetPsfPixels( )
{
  return psfPixels;
}


/* ---------------- GetRegionString ------------------------------------ */
string PsfOversamplingInfo::GetRegionString( )
{
  return regionString;
}


/* ---------------- GetOversamplingScale ------------------------------- */
int PsfOversamplingInfo::GetOversamplingScale( )
{
  return oversamplingScale;
}


/* ---------------- GetImageOffset ------------------------------------- */
void PsfOversamplingInfo::GetImageOffset( int &x0, int &y0 )
{
  x0 = X0_offset;
  y0 = Y0_offset;
}


/* ---------------- GetCorrectedRegionCoords --------------------------- */
std::tuple<int, int, int, int> PsfOversamplingInfo::GetCorrectedRegionCoords( )
{
  int  x1_region, x2_region, y1_region, y2_region;
  int  x1, x2, y1, y2;
  std::tie(x1_region, x2_region, y1_region, y2_region) = GetAllCoordsFromBracket(regionString);
  x1 = x1_region - X0_offset;
  x2 = x2_region - X0_offset;
  y1 = y1_region - Y0_offset;
  y2 = y2_region - Y0_offset;
  return std::make_tuple(x1, x2, y1, y2);
}


/* ---------------- GetNormalizationFlag ------------------------------- */
/// Returns true if PSF is meant to be normalized
bool PsfOversamplingInfo::GetNormalizationFlag( )
{
  return normalizePSF;
}



/* END OF FILE: psf_oversampling_info.cpp ------------------------------ */
