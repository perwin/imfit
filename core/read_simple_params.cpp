/*    Some utility functions to read in luminosity or other profiles ---
 * i.e., data in two or three columns (1st column = independent variable,
 * 2nd column = data [dependent variable], optional 3rd column = errors,
 * optional 4th column = mask) --- from text files.
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
 *    ReadDataFile(dataFileName, xVals, yVals, NULL, NULL);
 *
 *    We use NULL for yErrs & maskVals because we know the file only has two columns,
 * or else we want to ignore the 3rd & 4th columns.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#ifdef LINUX
#include <string.h>
#endif

#include "read_simple_params.h"
#include "utilities_pub.h"

using namespace std;

#define  MAXLINE   1024
#define  FILE_OPEN_ERR_STRING "\n   Couldn't open file \"%s\"\n\n"
#define  NO_DATA_ERR_STRING1 "\n   Missing data on line %ld of %s!\n\n"


long CountDataLines( string& fileName )
{
  FILE *file_ptr;
  char  buf[MAXLINE];
  long  counter = 0;
  
  if ((file_ptr = fopen(fileName.c_str(), "r")) == nullptr) {
    fprintf(stderr, FILE_OPEN_ERR_STRING, fileName.c_str());
    exit(-1);
  }

  // First, count how many data lines there are:
  while ( !feof(file_ptr) ) {
    if ( fgets(buf, MAXLINE - 1, file_ptr) != nullptr ) {
      // Count line if it does *not* starts with "#":
      if ( strncmp(buf, "#", 1) != 0 )
        counter++;
    }
  }
  fclose(file_ptr);
  return(counter);
}



long ReadSimpleParameterFile( string& fileName, long startDataRow, long endDataRow, 
                 double *xVals )
/*    Reads data from a file into 2--4 vectors, one for each column.
 * Lines in the file starting with "#" are assumed to be comments and are skipped.
 * Defaults to three-column input (e.g., x  y  y_err); if file has only two
 * columns --- or 3rd column isn't needed --- call the function with
 * yErrs = NULL.  Optionally, the 4th column can be read as "mask" values, but
 * only if the 3rd column is *also* read.
 *    Only data lines with indices >= startDataRow and <= endDataRow are
 * stored (C-style indexing, so first data row = 0).
 *
 *    Returns the number of saved rows.
 * 
 */
{
  FILE *file_ptr;
  char  buf[MAXLINE];
  string  currentLine;
  vector<string>  stringPieces;
  int  nPieces;
  long  nRow = 0;
  long  i = -1;   // counts which (input) *data* row we are on
  long  j = -1;   // counts which data row we are storing
  
  if ((file_ptr = fopen(fileName.c_str(), "r")) == nullptr) {
    fprintf(stderr, FILE_OPEN_ERR_STRING, fileName.c_str());
    exit(-1);
  }

  while ( !feof(file_ptr) ) {
    nRow++;
    if ( fgets(buf, MAXLINE - 1, file_ptr) != nullptr ) {
      currentLine = buf;
      // Process line if it does *not* starts with "#":
      if ( currentLine[0] != '#' ) {
        // this is a valid data line...
        i++;
        // Process line if it falls within user-specified bounds:
        if ( (i >= startDataRow) && (i <= endDataRow) ) {
          // ...and it's within the user-specified ranges
          j++;
          SplitString(currentLine, stringPieces);
          nPieces = stringPieces.size();
          
          // check for missing values
          if (nPieces < 1) {
            fprintf(stderr, NO_DATA_ERR_STRING1, nRow, fileName.c_str());
            exit(-1);
          }

          xVals[j] = strtod(stringPieces[0].c_str(), (char **)nullptr);
        }    
      }
    }
  }
  
  
  fclose(file_ptr);
  
  return(j + 1);
}
