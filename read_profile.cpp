/*    Some utility functions to read in luminosity or other profiles ---
 * i.e., data in two or three columns (1st column = independent variable,
 * 2nd column = data [dependent variable], optional 3rd column = errors)
 *  --- from text files.
 *    Assumes data in two columns or three, separated by spaces and/or tabs.
 *    Comment lines (beginning with "#") are skipped.
 *    
 *    HOW TO USE:
 *    Use CountDataLines() to find out how many data points there are, then
 * allocate vectors for the two columns, then use ReadDataFile() to put the
 * data into the vectors.
 *    Example:
 *    
 *    int  nDataVals;
 *    double  *xVals, yVals;
 *    
 *    nDataVals = CountDataLines(dataFileName);
 *    xVals = (double *)calloc( (size_t)nDataVals, sizeof(double) );
 *    yVals = (double *)calloc( (size_t)nDataVals, sizeof(double) );
 *    if ( (xVals == NULL) || (yVals == NULL) ) {
 *      fprintf(stderr, "\nFailure to allocate memory for input data!\n");
 *      exit(-1);
 *    }
 *    ReadDataFile(dataFileName, xVals, yVals, NULL);
 *
 *    We use NULL for yErrs because we know the file only has two columns,
 * or else we want to ignore the 3rd column.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "read_profile_pub.h"

#define  MAXLINE   1024
#define  FILE_OPEN_ERR_STRING "\n   Couldn't open file \"%s\"\n\n"
#define  NO_DATA_ERR_STRING1 "\n   Missing x data on line %ld of %s!\n\n"
#define  NO_DATA_ERR_STRING2 "\n   Missing y data on line %ld of %s!\n\n"
#define  NO_DATA_ERR_STRING3 "\n   Missing y_err data on line %ld of %s!\n\n"


long CountDataLines( std::string& fileName )
{
  FILE *file_ptr;
  char  buf[MAXLINE];
  long  counter = 0;
  
  if ((file_ptr = fopen(fileName.c_str(), "r")) == NULL) {
    fprintf(stderr, FILE_OPEN_ERR_STRING, fileName.c_str());
    exit(-1);
  }

  // First, count how many data lines there are:
  while ( !feof(file_ptr) ) {
    if ( fgets(buf, MAXLINE - 1, file_ptr) != NULL ) {
      // Count line if it does *not* starts with "#":
      if ( strncmp(buf, "#", 1) != 0 )
        counter++;
    }
  }
  fclose(file_ptr);
  return(counter);
}



long ReadDataFile( std::string& fileName, long startDataRow, long endDataRow, 
                 double *xVals, double *yVals, double *yErrs )
/*    Reads data from a file into 2 or 3 vectors, one for each column.
 * Lines in the file starting with "#" are assumed to be comments and are skipped.
 * Defaults to three-column input (e.g., x  y  y_err); if file has only two
 * columns --- or 3rd column isn't needed --- call the function with
 * yErrs = NULL.  Extra columns beyond the first three (or first two if
 * yErrs is set to NULL) are ignored.
 *    Only data lines with indices >= startDataRow and <= endDataRow are
 * stored (C-style indexing, so first data row = 0).  To read in all the
 * data, 
 * 
 */
{
  FILE *file_ptr;
  char  buf[MAXLINE];
  char  *xStr, *yStr, *eStr;
  long  nRow = 0;
  long  i = -1;   // counts which (input) *data* row we are on
  long  j = -1;   // counts which data row we are storing
  
  if ((file_ptr = fopen(fileName.c_str(), "r")) == NULL) {
    fprintf(stderr, FILE_OPEN_ERR_STRING, fileName.c_str());
    exit(-1);
  }

  while ( !feof(file_ptr) ) {
    nRow++;
    if ( fgets(buf, MAXLINE - 1, file_ptr) != NULL ) {
      // Process line if it does *not* starts with "#":
      if ( strncmp(buf, "#", 1) != 0 ) {
        // this is a valid data line...
        i++;
        // Process line if it falls within user-specified bounds:
        if ( (i >= startDataRow) && (i <= endDataRow) ) {
          // ...and it's within the user-specified ranges
          j++;
          xStr = strtok(buf, " \t");
          if (xStr == NULL) {
            fprintf(stderr, NO_DATA_ERR_STRING1, nRow, fileName.c_str());
            exit(-1);
          }
          xVals[j] = strtod(xStr, (char **)NULL);
          
          yStr = strtok(NULL, " \t");
          if (yStr == NULL) {
            fprintf(stderr, NO_DATA_ERR_STRING2, nRow, fileName.c_str());
            exit(-1);
          }
          if (yErrs == NULL) {
            yStr[strlen(yStr) - 1] = '\0';
            yVals[j] = strtod(yStr, (char **)NULL);
          } else {
            yVals[j] = strtod(yStr, (char **)NULL);
            eStr = strtok(NULL, " \t");
            if (eStr == NULL) {
              fprintf(stderr, NO_DATA_ERR_STRING3, nRow, fileName.c_str());
              exit(-1);
            }
            eStr[strlen(eStr) - 1] = '\0';
            yErrs[j] = strtod(eStr, (char **)NULL);
          }
          // xVals[i] = atof(xStr);
          // yVals[i] = atof(yStr);
          // printf("%f\t%f\n", xVals[i], yVals[i]);
        }    
      }
    }
  }
  
  
  fclose(file_ptr);
  
  return(j);
}
