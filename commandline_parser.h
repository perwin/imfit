/*   Public interfaces for command-line parser.
 */

#ifndef _COMMANDLINE_PARSER_H_
#define _COMMANDLINE_PARSER_H_

#include <string>
#include <vector>
#include <map>

using namespace std;


// Utility functions
void StripLeadingDashes( string& stringToModify );



// Option class
class OptionObject
{
  public:
    // Constructors and Destructors:
    OptionObject( );
    ~OptionObject( );
    
    // Public member functions:
    void DefineAsFlag( );   // specify that this is a "flag" options
    bool IsFlag( );    // is this a flag option?
    void SetFlag( );   // set flag value to true
    bool FlagSet( );   // has flag been set?
    void StoreTarget( const char targString[] );   // e.g., store filename pointed to by option
    bool TargetSet( );   // has target been set (if option is not flag)?
    string& GetTargetString( );

  private:
    // Private member functions:
    
    // Data members:
    bool  isFlag, flagSet, targetSet;
    string  targetString;
};





// Parser class
class CLineParser
{
  public:
    // Constructors and Destructors:
    CLineParser( );
    ~CLineParser( );
    
    // Public member functions:
    void PrintUsage( );
    void UnrecognizedAreErrors( );   // interpret unrecognized flags/options as errors
    void AddFlag( string shortFlagString );
    void AddFlag( string shortFlagString, string longFlagString );
    void AddOption( string shortOptString );
    void AddOption( string shortOptString, string longOptString );
    void AddUsageLine( string usageLine );
    int ParseCommandLine( int argc, char *argv[] );
    bool CommandLineEmpty( );
    bool FlagSet( string flagName );
    bool OptionSet( string optName );
    string& GetTargetString( string optName );
    int nArguments( );
    string& GetArgument( int n );

  private:
  // Private member functions:

  // Data members:
  map<string, OptionObject *>  optMap;   // 	data: map<string, *optObject> of option-objects
  vector<OptionObject *> optObjPointers;   // keep track of all discrete OptionObject pointers
  vector<string>  usageStrings;   // 	data: map<string, *optObject> of option-objects
  vector<string>  argStrings;     // 	data: vector<string> of argument strings
  int  verboseLevel;
  bool  commandLineEmpty;        // true if user supplied *no* options/flags *or* arguments
  bool  ignoreUnrecognized;      // if we encounter an unrecognized option/flag, do we
                                 // ignore it (= just print a warning & continue processing)?
  string  errorString1;
};


#endif  // _COMMANDLINE_PARSER_H_
