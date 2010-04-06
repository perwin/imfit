/* FILE: utilities.cpp --------------------------------------------- */
/*   Several utility routines used by imfit, makeimage, etc.
 */

#include <ctype.h>   /* for isdigit() */
#include <stdio.h>
#include <stdlib.h>  /* for exit() */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//#include "messages_and_defs.h"
#include "utilities_pub.h"
#include "mpfit_cpp.h"
#include "statistics.h"

using namespace std;


/* ------------------- Function Prototypes ----------------------------- */
/* Local Functions: */



/* ---------------- FUNCTION: SplitString() ------------------------ */
// This function tokenizes a string, splitting it into substrings using
// delimiters as the separator (delimiters can be more than one character, in
// which case all of them can serve as delimiters).  The substrings are
// added to the user-supplied vector<string> tokens.
// The default value for delimiter is "\t ", meaning tabs and spaces.
void SplitString( const string& str, vector<string>& tokens, const string& delimiters )
{
  // Skip delimiters at beginning.
  string::size_type  lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type  pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}


/* ---------------- FUNCTION: ChopComment() ------------------------ */
// This function removes the remainder of line after a comment character
// (latter is specified by delimiter, which defaults to '#')
void ChopComment( string& inputString, char delimiter )
{
  string::size_type  loc;
  
  loc = inputString.find(delimiter, 0);
  inputString = inputString.substr(0, loc);
}


/* ---------------- FUNCTION: TrimWhitespace() --------------------- */
// This function removes leading and trailing whitespace from a string; if
// the string is *all* whitespace, then it converts the input string to an
// empty string.  ("Whitespace" = spaces or tabs)
void TrimWhitespace(string& stringToModify)
{
  if (stringToModify.empty())
    return;

  string::size_type  startIndex = stringToModify.find_first_not_of(" \t");
  string::size_type  endIndex = stringToModify.find_last_not_of(" \t");
  if (startIndex == endIndex)
    stringToModify.clear();
  else
    stringToModify = stringToModify.substr(startIndex, (endIndex - startIndex + 1) );
}



/* ---------------- FUNCTION: FileExists() ------------------------- */

bool FileExists(const char * filename)
{
  return ifstream(filename);
}



/* ---------------- FUNCTION: CommandLineError() ------------------- */

void CommandLineError( char errorString[] )
{

  printf("Error in command line:\n   %s\nExiting...\n",
	 errorString);
  exit(1);
}



/* ---------------- FUNCTION: NotANumber() ------------------------- */
// Possible cases:
//    0, 0.0, 0.1, .1
//    -0.1, -.1?
//    -1
bool NotANumber( char theString[], int index, int restriction )
{
  int  theCharacter = theString[index];

  switch (restriction) {
    case kAnyInt:
      if (theCharacter == '-')
        return NotANumber( theString, index + 1, kAnyInt );
      else
        return (bool)( ! isdigit(theCharacter) );
    
    case kPosInt:
      if ( isdigit(theCharacter) && (theCharacter != '0') )
        return false;
      else
        return true;
    
    case kAnyReal:
      switch (theCharacter) {
        case '-':
          return NotANumber( theString, index + 1, kAnyReal );
        case '.':
          return NotANumber( theString, index + 1, kAnyInt );
        default:
          return (bool)( ! isdigit(theCharacter) );
      }  /* end switch (theCharacter) */
    
    case kPosReal:
      // THIS STILL NEEDS WORK!
      switch (theCharacter) {
        case '-':
          return true;
        case '.':
          return NotANumber( theString, index + 1, kAnyInt );
        default:
          return (bool)( ! isdigit(theCharacter) );
      }  /* end switch (theCharacter) */
    
    default:
      return true;
  }  /* end switch (restriction) */
}




/* END OF FILE: utilities.cpp -------------------------------------- */
