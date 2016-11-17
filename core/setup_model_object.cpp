// Code for setting up ModelObject instances with data

#include <vector>
#include <stdlib.h>

#include "setup_model_object.h"
#include "options_base.h"
#include "model_object.h"

using namespace std;


ModelObject* SetupModelObject( OptionsBase *options, vector<int> nColumnsRowsVector, 
					double *dataPixels, double *psfPixels, double *maskPixels, 
					double *errorPixels, double *psfOversampledPixels, 
					vector<int> xyOversamplePos )
{
  ModelObject *newModelObj;
  int  status;
  
  newModelObj = new ModelObject();

  int  nColumns = nColumnsRowsVector[0];
  int  nRows = nColumnsRowsVector[1];
  status = newModelObj->SetupModelImage(nColumns, nRows);
  if (status < 0) {
    fprintf(stderr, "*** ERROR: Failure in ModelObject::SetupModelImage!\n\n");
    exit(-1);
  }

  return newModelObj;
}
