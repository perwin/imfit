/* Header file for functions which read in data */

#ifndef _READ_SIMPLE_PARAMS_H_
#define _READ_SIMPLE_PARAMS_H_

#include <string>

long CountDataLines( std::string& fileName );

long ReadSimpleParameterFile( std::string& fileName, long startDataRow, long endDataRow, 
			                 double *xVals );


#endif /* _READ_SIMPLE_PARAMS_H_ */
