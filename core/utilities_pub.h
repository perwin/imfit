/*! \file
   \brief  Utility functions for string processing, testing file existence, etc.
 */

#ifndef _UTILITIES_PUB_H_
#define _UTILITIES_PUB_H_

#include <string>
#include <vector>

#include "mpfit.h"
#include "model_object.h"

using namespace std;

/* constants for use parameter "restriction" when calling NotANumber(): */
const int kAnyInt     =    0;
const int kNonzeroInt =    1;
const int kPosInt     =    2;
const int kAnyReal    =    3;
const int kPosReal    =    4;

/* constants for use with PrepareImageComments: */
const int HDR_MAKEIMAGE = 0;
const int HDR_MODELIMAGE = 1;
const int HDR_RESIDUALIMAGE = 2;
const int HDR_WEIGHTIMAGE = 3;


// String-processing functions

/// like fprintf, but returns a string instead of writing to a file
string PrintToString( const char *fmt, ... );


/// \brief Generates vector of strings containing program name, time stamp,
///        and command-line invocation, intended to be written to start of
///        best-fit output-parameter files, etc.
void MakeOutputHeader( vector<string> *headerLines, const string& programName, 
						const int argc, char *argv[] );

/// \brief Generates vector of strings containing program name, config file name,
///        and PSF image (if any), intended to be written to output FITS image header(s).
void PrepareImageComments( vector<string> *comments, const string &programName, 
                           const string &configFileName, bool psfUsed, 
                           const string &psfFileName, int mode=HDR_MAKEIMAGE,
                           const string &dataFileName="" );

/// \brief Splits a string and returns substrings ("tokens").
///
/// The vector tokens is cleared before adding the substrings.
void SplitString( const string& str, vector<string>& tokens, 
									const string& delimiters = "\n\t " );

/// \brief Same as SplitString, but the pieces of the input string are *added* to the
///        tokens vector, instead of the tokens vector being cleared first
void SplitStringAdd( const string& str, vector<string>& tokens, 
									const string& delimiters = "\t " );

/// Modifies inputString to remove all characters from the first appearance of the delimiter onwards
void ChopComment( string& inputString, char delimiter = '#' );

/// Removes leading and trailing whitespace ("\t ") from a string
void TrimWhitespace( string& stringToModify );


/// Removes leading and trailing square brackets from a string
void StripBrackets( const string& inputFilename, string& strippedFilename );


/// Extracts all coordinates from string of form "[x1:x2,y1:y1]"
void GetAllCoordsFromBracket( const string& bracketString, int *x1, int *x2,
                           int *y1, int *y2 );
void GetAllCoordsFromBracket2( const string& bracketString, vector<int>& x1, vector<int>& x2,
                           vector<int>& y1, vector<int>& y2 );

/// Extracts coordinates x1 and y1 from string of form "[x1:x2,y1:y1]"
void GetStartCoordsFromBracket( const string& bracketString, int *x1, int *y1,
                           const string& fileName );

/// \brief Determines "start" coordinates (lower-left corner) for an image,
///        accounting for any image section provided
void GetPixelStartCoords( const string& inputFilename, int *xStart, int *yStart );


/// \brief Checks to see if image (FITS) file exists (ignores any image-section specification;
///        filenames beginning "ftp:" or "http:" are assumed to exist.
bool ImageFileExists(const char * filename);

/// Checks to see if (local) file exists
bool FileExists(const char * filename);


/// Returns string with current date-and-time (Dow Mon dd hh:mm:ss yyy)
char * TimeStamp( void );


/// Prints "Error in command line:" + errorString, then calls exit(1)
void CommandLineError( char errorString[] );


/// Returns true if string is *not* a number of correct subtype (restriction)
bool NotANumber( const char theString[], int index, int restriction );


/// \brief Returns 1 if entire string is a valid number (positive or negative, including floating-point);
///        returns 0 if not.
int IsNumeric( const char *aString );



#endif /* _UTILITIES_PUB_H_ */
