/* FILE: parameter_utils.cpp --------------------------------------------- */

// Copyright 2014 by Peter Erwin.
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


// These are utility functions for dealing with parameter vectors (1-D double-precision
// arrays containing parameter values). In particular, these functions help implement
// "fixed" parameters for minimization routines that don't handle them natively
// (e.g., the levmar function dlevmar_bc_dif and its cousins).


#include <stdio.h>
#include "param_struct.h"   // for mp_par structure


// Given an mp_par array parameterLimits (containing nTot mp_par structs), this function
// counts and returns the number of parameters which are *fixed*.
int CountFixedParams( mp_par *parameterLimits, int nTot )
{
  int  nFixed = 0;
  
  for (int i = 0; i < nTot; i++)
    if (parameterLimits[i].fixed == 1)
      nFixed++;
  
  return nFixed;
}


// Given a full set of parameters in params[], return a "condensed" copy of params[] with
// all *fixed* parameters removed, stored in freeParams[]. If there are no fixed parameters,
// then freeParams[] will be identical to params[].
// The nParamsTot-length bool array fixedPars[] specifies whether each parameter is 
// fixed (= true) or free (= false)
void CondenseParamVector( double params[], double freeParams[], int nParamsTot, 
							int nParamsFree,  bool fixedPars[] )
{
  int  j = 0;
  for (int i = 0; i < nParamsTot; i++) {
    if (! fixedPars[i]) {
      freeParams[j] = params[i];
      j++;
    }
  }
  if (j != nParamsFree)
    fprintf(stderr, "ERROR in CondenseParamVector: j (%d) != nParamsFree (%d)",
    		j, nParamsFree);
}


// Same as CondenseParamVector, but working on two vectors of parameter limits
// (lower and upper limits). The input vectors inputMinParamLims and 
// inputMaxParamLims contain lower and upper limits for all nParamsTot
// parameters (though some of these limits may be undefined). The output is the
// two nFree-length vectors lowerLims and upperLims, which have the lower and upper
// limits for only the free parameters. If there are no fixed parameters,
// then lowerLims and upperLims will be identical to the corresponding input
// vectors.
void CondenseParamLimits( double inputMinParamLims[], double inputMaxParamLims[], 
							double lowerLims[], double upperLims[], int nParamsTot, 
							int nParamsFree, bool fixedPars[] )
{
  int  j = 0;
  for (int i = 0; i < nParamsTot; i++) {
    if (! fixedPars[i]) {
      lowerLims[j] = inputMinParamLims[i];
      upperLims[j] = inputMaxParamLims[i];
      j++;
    }
  }
  if (j != nParamsFree)
    fprintf(stderr, "ERROR in CondenseParamLimits: j (%d) != nParamsFree (%d)",
    		j, nParamsFree);
}



// Given a full-sized vector of original parameter values (originalParams; nParamsTot
// in length) and an input "condensed" vector condensedNewParams (nParamsFree in
// length), this functions populates outputParams with the values from originalParams
// for those parameters which are fixed, and the values from condensedNewParams
// for the free parameters.  If there are not fixed parameters, then outputParams
// will be identical to condensedNewParams.
// The nParamsTot-length bool array fixedPars[] specifies whether each parameter is 
// fixed (= true) or free (= false)
void ExpandParamVector( double originalParams[], double condensedNewParams[],
							double outputParams[], int nParamsTot, int nParamsFree,
							bool fixedPars[] )
{
  int  k;
  
  k = 0;
  for (int i = 0; i < nParamsTot; i++) {
    if (fixedPars[i] == true) {
      outputParams[i] = originalParams[i];
    }
    else {
      outputParams[i] = condensedNewParams[k];
      k++;
    }
  }
}

/* END OF FILE: parameter_utils.cpp -------------------------------------- */

