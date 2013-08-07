/*    Utility functions for testing file existence, reading/parsing command-line 
 * options, etc.
 *
 */

#ifndef _UTILITIES_PUB_H_
#define _UTILITIES_PUB_H_

#include <string>
#include <vector>

#include "mpfit_cpp.h"
#include "model_object.h"

using namespace std;

/* constants for use parameter "restriction" when calling NotANumber(): */
#define kAnyInt          0
#define kNonzeroInt      1
#define kPosInt          2
#define kAnyReal         3
#define kPosReal         4


// String-processing functions

// Splits a string and returns substrings as elements of tokens (tokens is cleared
// before adding the substrings).
void SplitString( const string& str, vector<string>& tokens, 
									const string& delimiters = "\t " );

// Same as SplitString, but the pieces of the input string are *added* to the
// tokens vector, instead of the tokens vector being cleared first
void SplitStringAdd( const string& str, vector<string>& tokens, 
									const string& delimiters = "\t " );

// Removes remainder of string after first occurance of delimiter.
void ChopComment( string& inputString, char delimiter = '#' );

// Removes leading and trailing whitespace ("\t ")from a string
void TrimWhitespace( string& stringToModify );


void StripBrackets( const string& inputFilename, string& strippedFilename );


void GetPixelStartCoords( const string& inputFilename, int *xStart, int *yStart );


bool ImageFileExists(const char * filename);

bool FileExists(const char * filename);


char * TimeStamp( void );


void CommandLineError( char errorString[] );


bool NotANumber( const char theString[], int index, int restriction );



#endif /* _UTILITIES_PUB_H_ */
