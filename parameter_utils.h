// Header for for utility functions dealing with parameter vectors

#ifndef _PARAMETER_UTILS_H_
#define _PARAMETER_UTILS_H_


#include "param_struct.h"   // for mp_par structure


int CountFixedParams( mp_par *parameterLimits, int nTot );

void CondenseParamVector( double params[], double freeParams[], int nParamsTot, 
							int nParamsFree,  bool fixedPars[] );

void CondenseParamLimits( double inputMinParamLims[], double inputMaxParamLims[], 
							double lowerLims[], double upperLims[], int nParamsTot, 
							int nParamsFree, bool fixedPars[] );

void ExpandParamVector( double originalParams[], double condensedNewParams[],
							double outputParams[], int nParamsTot, int nParamsFree,
							bool fixedPars[] );


#endif /* _PARAMETER_UTILS_H_ */
