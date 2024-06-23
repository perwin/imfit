// Code for parsing and re-arranging parameter vectors in ModelObjectMultImage

#ifndef _PARAMVECTOR_PROCESSING_H_
#define _PARAMVECTOR_PROCESSING_H_

#include <vector>
#include <tuple>

using namespace std;

int ExtractImageParams( double inputParamsVector[], int imageNumber, 
						double& pixScale, double& rotation, double& intensityScale,
    					double& X0, double& Y0 );

std::tuple<double, double, int> CalculateOffset_X0Y0( double X0_0_ref, double Y0_0_ref, 
								double X0_n_ref, double Y0_n_ref, double X0_0_im, 
								double Y0_0_im, double pixScale_im, double rotation_im );

void AssembleParamsForImage( double externalInputParamsVect[], int nInputParamsTot,
							int nImagesTot, int imageNumber, int nFunctions, 
							vector<int>& paramSizes, bool fblockStartFlags[], int nFuncBlocks, 
							double outputModelParams[] );


#endif  // _PARAMVECTOR_PROCESSING_H_
