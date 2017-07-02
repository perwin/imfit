/* FILE: psf_oversampling_info.cpp ------------------------------------- */
/* 
 * Class definition for PsfOversamplingInfo.
 */

// Copyright 2017 by Peter Erwin.
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
#include <vector>

using namespace std;

#include "psf_oversampling_info.h"


/* ---------------- Definitions ---------------------------------------- */



/* ---------------- CONSTRUCTOR (default) ------------------------------ */

PsfOversamplingInfo::PsfOversamplingInfo( )
{
  nColumns_psf = nRows_psf = 0;
}


/* ---------------- CONSTRUCTOR ---------------------------------------- */

PsfOversamplingInfo::PsfOversamplingInfo( double *inputPixels, int nCols, int nRows, 
										int scale, string inputRegionString )
{
  psfPixels = inputPixels;
  nColumns_psf = nCols;
  nRows_psf = nRows;
  oversamplingScale = scale;
  regionString = inputRegionString;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

PsfOversamplingInfo::~PsfOversamplingInfo( )
{
  ;
}



/* ---------------- AddRegionString ------------------------------------ */
void PsfOversamplingInfo::AddRegionString( string inputRegionString )
{
  regionString = inputRegionString;
}


/* ---------------- AddPsfPixels --------------------------------------- */
void PsfOversamplingInfo::AddPsfPixels( double *inputPixels, int nCols, int nRows )
{
  psfPixels = inputPixels;
  nColumns_psf = nCols;
  nRows_psf = nRows;
}


/* ---------------- AddOversamplingScale ------------------------------- */
void PsfOversamplingInfo::AddOversamplingScale( int scale )
{
  oversamplingScale = scale;
}


/* ---------------- GetNColumns ---------------------------------------- */
int PsfOversamplingInfo::GetNColumns( )
{
  return nColumns_psf;
}


/* ---------------- GetNRows ------------------------------------------- */
int PsfOversamplingInfo::GetNRows( )
{
  return nRows_psf;
}


/* ---------------- GetPsfPixels --------------------------------------- */
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



/* END OF FILE: psf_oversampling_info.cpp ------------------------------ */
