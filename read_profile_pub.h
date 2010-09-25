/* Header file for functions which read in data */

#ifndef _READ_PROFILE_H_
#define _READ_PROFILE_H_

#include <string>

long CountDataLines( std::string& fileName );

long ReadDataFile( std::string& fileName, long startDataRow, long endDataRow, 
                        double *xVals, double *yVals, double *yErrs, double *maskVals );


#endif /* _READ_PROFILE_H_ */
