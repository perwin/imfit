/* FILE: commandline_parser.cpp ---------------------------------------- */
/* VERSION 0.5
 *
 *
 * Based loosely on Kishan Thomas' AnyOption class (2006 version).
 *
 * 7--8 Feb 2011: Created.
 *
 */

// Copyright 2011, 2012, 2013 by Peter Erwin.
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


#include <stdio.h>
#include <string>
#include <vector>
#include <map>

#include "commandline_parser.h"
#include "utilities_pub.h"

using namespace std;

//const char  *equalsSign = "=";


/* UTILITY FUNCTIONS */

/* ---------------- FUNCTION: StripLeadingDashes() ---------------------- */
// This function removes leading and trailing whitespace from a string; if
// the string is *all* whitespace, then it converts the input string to an
// empty string.  ("Whitespace" = spaces or tabs)
void StripLeadingDashes( string& stringToModify )
{
  if (stringToModify.empty())
    return;

  string::size_type  startIndex = stringToModify.find_first_not_of("-");
  if (startIndex == string::npos)   // nothing but dashes in this string!
    stringToModify = "";
  else
    stringToModify = stringToModify.substr(startIndex);
}


/* *** DEFINITIONS FOR OPTIONOBJECT CLASS *** */

/* ---------------- CONSTRUCTOR ---------------------------------------- */

OptionObject::OptionObject( )
{
  isFlag = false;
  flagSet = false;
  targetSet = false;
  targetString = "";
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

OptionObject::~OptionObject( )
{
  ;
//  printf("OptionObject being destroyed!\n");
}


/* ---------------- DefineAsFlag --------------------------------------- */
void OptionObject::DefineAsFlag(  )
{
  isFlag = true;
}


/* ---------------- IsFlag --------------------------------------------- */
bool OptionObject::IsFlag(  )
{
  return isFlag;
}


/* ---------------- SetFlag -------------------------------------------- */
void OptionObject::SetFlag(  )
{
  flagSet = true;
}


/* ---------------- FlagSet -------------------------------------------- */
bool OptionObject::FlagSet(  )
{
  return flagSet;
}


/* ---------------- StoreTarget ---------------------------------------- */
void OptionObject::StoreTarget( const char targString[] )
{
  targetString = targString;
  targetSet = true;
}


/* ---------------- TargetSet ------------------------------------------ */
bool OptionObject::TargetSet(  )
{
  return targetSet;
}


/* ---------------- GetTarget ------------------------------------------ */
string& OptionObject::GetTargetString(  )
{
  return targetString;
}





/* *** DEFINITIONS FOR CLINEPARSER CLASS *** */

/* ---------------- CONSTRUCTOR ---------------------------------------- */

CLineParser::CLineParser( )
{
  ignoreUnrecognized = true;
  commandLineEmpty = true;
  verboseLevel = 0;
  errorString1 = "<NULL>";
}


// Note: we use the vector optObjPointers as a way of freeing the memory
// associated with the OptionObjects, rather than using the map optMap, because
// in the latter more than one key can point to the same pointer (so if we
// were to iterate through by key, we could end up trying to free memory
// that's already been freed).  In optObjPointers, each element is a separate,
// discrete pointer to a separate, discrete OptionObject.

/* ---------------- DESTRUCTOR ----------------------------------------- */

CLineParser::~CLineParser( )
{
//  printf("CLineParser object being destroyed!\n");
  
  // free the memory occupied by the OptionObject instances by calling delete
  // on the pointers to each object
  for (int i = 0; i < (int)optObjPointers.size(); i++)
    delete optObjPointers[i];
}



void CLineParser::PrintUsage( )
{
  int  nStrings = (int)usageStrings.size();
  for (int i = 0; i < nStrings; i++)
    printf("%s\n", usageStrings[i].c_str());
}


// Call this method to specify that the parser should treat unrecognized
// options/flags as errors (aborting the parsing process and returning -1)
// instead of just issuing a warning.
void CLineParser::UnrecognizedAreErrors( )
{
  ignoreUnrecognized = false;
}


// Note that the pointers to OptionObjects created in the following methods will be
// stored in the built-in optMap data member; we use the optObjPointers vector to
// keep track of them for purposes of freeing up the allocated memory when this
// CLineParser object is destroyed.

void CLineParser::AddFlag( string shortFlagString )
{
  OptionObject *newOptionObj = new OptionObject;
  optObjPointers.push_back(newOptionObj);
  newOptionObj->DefineAsFlag();
  optMap[shortFlagString] = newOptionObj;
}


void CLineParser::AddFlag( string shortFlagString, string longFlagString )
{
  OptionObject *newOptionObj = new OptionObject;
  optObjPointers.push_back(newOptionObj);
  newOptionObj->DefineAsFlag();
  optMap[shortFlagString] = newOptionObj;
  optMap[longFlagString] = newOptionObj;
}


void CLineParser::AddOption( string shortOptString )
{
  OptionObject *newOptionObj = new OptionObject;
  optObjPointers.push_back(newOptionObj);
  optMap[shortOptString] = newOptionObj;
}


void CLineParser::AddOption( string shortOptString, string longOptString )
{
  OptionObject *newOptionObj = new OptionObject;
  optObjPointers.push_back(newOptionObj);
  optMap[shortOptString] = newOptionObj;
  optMap[longOptString] = newOptionObj;
}


void CLineParser::AddUsageLine( string usageLine )
{
  usageStrings.push_back(usageLine);
}


int CLineParser::ParseCommandLine( int argc, char *argv[] )
{
  string  currentString;
  vector<string>  stringPieces;
  OptionObject  *currentOpt;
  int  i;
  
  i = 0;   // program name is first value (i = 0); we'll skip it
  while (i < (argc - 1)) {
    i++;
    if (verboseLevel > 1)
      printf("i = %d: \"%s\"\n", i, argv[i]);
    if (argv[i][0] == '-') {
      // flag or option
      commandLineEmpty = false;
      currentString = argv[i];
      StripLeadingDashes(currentString);
      if (currentString.size() < 1) {
        fprintf(stderr, "WARNING: isolated \"-\" or \"--\" found on command line!\n");
        return -1;
      }
      if (verboseLevel > 1)
        printf("\tflag or option: %s\n", currentString.c_str());
      // chop the string up if there's an "=" in the middle
      SplitString(currentString, stringPieces, "=");
      if (optMap.count(stringPieces[0]) > 0) {
        // OK, this is a valid flag or option
        currentOpt = optMap[stringPieces[0]];
        if (currentOpt->IsFlag()) {
          // It's a flag, so set it
          currentOpt->SetFlag();
          continue;
        }
        else {
          // It's an option with a target
          // check for "=" format
          if (stringPieces.size() > 1) {
            currentOpt->StoreTarget(stringPieces[1].c_str());
            continue;
          }
          else {   // no equals sign, so should be "-opt target" format
            if ((i + 1) < argc) {
              i++;
              // OPTION: check if target starts with "-" and warn user
              currentOpt->StoreTarget(argv[i]);
              if (verboseLevel > 1)
                printf("\tstoring target \"%s\"...\n", argv[i]);
            }
            else {
              // we need a target for this option, but we've run out of command-line!
              fprintf(stderr, "WARNING: option \"%s\" expects a following argument!\n",
                      argv[i]);
              return -1;
            }
          }
        }
      }
      else {
        fprintf(stderr, "WARNING: Unrecognized command-line option/flag \"%s\"!\n", argv[i]);
        if (! ignoreUnrecognized)
          return -1;
      }
    }
    else {
      // it's an argument, so store it in argStrings vector
      commandLineEmpty = false;
      if (verboseLevel > 1)
        printf("argStrings.size() = %d\n", (int)argStrings.size());
      argStrings.push_back(argv[i]);
      if (verboseLevel > 1)
        printf("... after storing \"%s\", argStrings.size() = %d\n", argv[i], (int)argStrings.size());
    }
  }
  
  return 0;
}


// Checks if associated flag was set
bool CLineParser::CommandLineEmpty( )
{
  return commandLineEmpty;
}




// Checks if associated flag was set
bool CLineParser::FlagSet( string flagName )
{
  if (optMap.count(flagName) > 0)
    return optMap[flagName]->FlagSet();
  else {
    fprintf(stderr, "\nERROR: \"%s\" is not an assigned flag!\n", flagName.c_str());
    return false;
  }
}


// Checks if target was supplied for the associated option
bool CLineParser::OptionSet( string optName )
{
  if (optMap.count(optName) > 0)
    return optMap[optName]->TargetSet();
  else {
    fprintf(stderr, "\nERROR: \"%s\" is not an assigned option!\n", optName.c_str());
    return false;
  }
}


// Returns the stored target string for the associated option
string& CLineParser::GetTargetString( string optName )
{
  if (optMap.count(optName) > 0)
    return optMap[optName]->GetTargetString();
  else {
    fprintf(stderr, "\nERROR: \"%s\" is not an assigned option!\n", optName.c_str());
    return errorString1;
  }
}


int CLineParser::nArguments( )
{
  return (int)argStrings.size();
}


string& CLineParser::GetArgument( int n )
{
  return argStrings[n];
}




/* END OF FILE: commandline_parser.cpp --------------------------------- */
