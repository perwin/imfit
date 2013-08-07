/* FILE: utilities.cpp --------------------------------------------- */
/*   Several utility routines used by imfit, makeimage, etc.
 */

// Copyright 2010, 2011, 2012, 2013 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.


#include <ctype.h>   /* for isdigit() */
#include <stdio.h>
#include <stdlib.h>  /* for exit() */
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

//#include "messages_and_defs.h"
#include "utilities_pub.h"
#include "mpfit_cpp.h"
#include "statistics.h"

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

  tokens.clear();

  // Skip delimiters at beginning.
  string::size_type  lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type  pos = str.find_first_of(delimiters, lastPos);
  if (pos == string::npos) {
    // oops, no delimiters in this string, so just return the whole string
    tokens.push_back(str);
  }
  else {
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
}


/* ---------------- FUNCTION: SplitStringAdd() --------------------- */
// Same as SplitString, but the pieces of the input string are *added* to the
// tokens vector, instead of the tokens vector being cleared first
void SplitStringAdd( const string& str, vector<string>& tokens, const string& delimiters )
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



/* ---------------- FUNCTION: StripBrackets() ---------------------- */

void StripBrackets( const string& inputFilename, string& strippedFilename )
{
  strippedFilename = inputFilename.c_str();
  ChopComment(strippedFilename, '[');
}



/* ---------------- FUNCTION: GetCoordsFromBracket() ---------------- */
// Given a string of the form "x1:x2,y1:y2", return x1 and y1
// Special case: "*,y1:y2" ==> return 1, y1
// Special case: "x1:x2,*" ==> return x1, 1
void GetCoordsFromBracket( const string& bracketString, int *x1, int *y1,
                           const string& fileName )
{
  vector<string>  sectionPieces, subsectionPieces_x, subsectionPieces_y;
  const string star = string("*");

  // default values indicating errors:
  *x1 = 0;
  *y1 = 0;
  
  SplitString(bracketString, sectionPieces, ",");
  // handle the x part of the section specification
  if (sectionPieces[0] == star)
    *x1 = 1;
  else {
    SplitString(sectionPieces[0], subsectionPieces_x, ":");
    if (subsectionPieces_x.size() != 2) {
      printf("\nWARNING1: Incorrect image section format in \"%s\"!\n",
    					fileName.c_str());
      printf("\t\"%s\"\n", bracketString.c_str());
      return;
    }
    *x1 = atoi(subsectionPieces_x[0].c_str());
  }
  // handle the y part of the section specification
  if (sectionPieces[1] == star)
    *y1 = 1;
  else {
    SplitString(sectionPieces[1], subsectionPieces_y, ":");
    if (subsectionPieces_y.size() != 2) {
      printf("\nWARNING2: Incorrect image section format in \"%s\"!\n",
    					fileName.c_str());
      printf("\t\"%s\"\n", bracketString.c_str());
      return;
    }
    *y1 = atoi(subsectionPieces_y[0].c_str());
  }
}


/* ---------------- FUNCTION: GetPixelStartCoords() ---------------- */

void GetPixelStartCoords( const string& inputFilename, int *xStart, int *yStart )
{
  string::size_type  loc1, loc2, loc3, loc4;
  int  nPieces;
  string  sectionSubstring;
  vector<string>  sectionPieces, subsectionPieces_x, subsectionPieces_y;
  const string star = string("*");
  bool  twoSections = false;
  
  // default values indicating errors:
  *xStart = 0;
  *yStart = 0;
  
  loc1 = inputFilename.find('[', 0);
  if (loc1 == string::npos) {
    // no image section specified, so we're using the entire image
    *xStart = 1;
    *yStart = 1;
    return;
  }
  
  // OK, if we get here, then there's apparently an image section
  loc2 = inputFilename.find(']', loc1);
  if (loc2 == string::npos) {
    printf("\nWARNING: Incorrect image section format in \"%s\"!\n",
    				inputFilename.c_str());
    return;
  }
  // check for possible second bracket group
  loc3 = inputFilename.find("[", loc2);
  if (loc3 != string::npos) {
    // OK, there's more than one set of []
    loc4 = inputFilename.find(']', loc3);
    if (loc4 == string::npos) {
      printf("\nWARNING: Incorrect image section format in \"%s\"!\n",
      				inputFilename.c_str());
      return;
    }
    twoSections = true;
  }
      
  // extract what's inside the first (and possibly only) []
  sectionSubstring = inputFilename.substr(loc1 + 1, loc2 - loc1 - 1);
  SplitString(sectionSubstring, sectionPieces, ",");
  nPieces = sectionPieces.size();
  // two valid possibilites: an image section (nPieces = 2) or an extension number
  // (nPieces = 1)
  if (nPieces == 2) {
    // apparently an image section
    GetCoordsFromBracket(sectionSubstring, xStart, yStart, inputFilename);
    // we found a valid (or invalid) image section; ignore anything else...
    return;
  }
  // if we don't have an image section, we need to have an image extension number
  if ((nPieces != 1) || (sectionSubstring.size() < 1)) {
    printf("\nWARNING: Incorrect image section format in \"%s\"!\n",
  					inputFilename.c_str());
    return;
  }
  if (twoSections == false) {
    // OK, just an image extension and nothing else
    *xStart = 1;
    *yStart = 1;
    return;
  }
      
  if (twoSections == true) {
    // OK, if we get here, there are two bracket groups, and the first was
    // apparently an extension number, so we should expect the second group
    // to be a proper image section
    sectionSubstring = inputFilename.substr(loc3 + 1, loc4 - loc3 - 1);
    SplitString(sectionSubstring, sectionPieces, ",");
    nPieces = sectionPieces.size();
    if (nPieces != 2) {
      printf("\nWARNING: Incorrect image section format in \"%s\"!\n",
    					inputFilename.c_str());
      return;
    }
    // apparently an image section
    GetCoordsFromBracket(sectionSubstring, xStart, yStart, inputFilename);
  }

}



/* ---------------- FUNCTION: ImageFileExists() -------------------- */
// Function which tests for the existence of an image file, with the following
// special cases:
//    1. If filename begins with "ftp:" or "http:", we assume it exists
//    2. Trailing image specifications (e.g. "name.fits[100:200, 100:200]")
// are ignored, since they are not part of the on-disk filename

bool ImageFileExists(const char * filename)
{
  string  ftpString("ftp://");
  string  httpString("http://");
  string  filenameStr(filename);
  string  baseImageFileName;
  
  // Check for possible ftp:// or http://
  if ((filenameStr.find(ftpString) != string::npos) || 
  		(filenameStr.find(httpString) != string::npos)) {
    return true;
  }
  
  StripBrackets(filenameStr, baseImageFileName);
  return FileExists(baseImageFileName.c_str());
}



/* ---------------- FUNCTION: FileExists() ------------------------- */

bool FileExists(const char * filename)
{
  return ifstream(filename);
}



/* ---------------- FUNCTION: TimeStamp() -------------------------- */

char * TimeStamp( void )
{  
  time_t  currentTime;
  char *dateString;
  
  currentTime = time(NULL);
  dateString = ctime(&currentTime);
  // Hack! Shift the null-termination up one character to knock out the
  // '\n' which otherwise ends the string.  This works because the
  // output of ctime() is supposed to be a 26-character string.
  dateString[24] = '\0';
  
  return dateString;
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
bool NotANumber( const char theString[], int index, int restriction )
{
  int  theCharacter = theString[index];

  switch (restriction) {
    case kAnyInt:
      if (theCharacter == '-')
        return NotANumber( theString, index + 1, kAnyInt );
      else
        return (bool)( ! isdigit(theCharacter) );
    
    case kNonzeroInt:
      if (theCharacter == '-')
        return false;
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
