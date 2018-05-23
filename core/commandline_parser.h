/** @file
    \brief Public interfaces for command-line parser, including class declarations
           for OptionObject and CLineParser.

 */

#ifndef _COMMANDLINE_PARSER_H_
#define _COMMANDLINE_PARSER_H_

#include <string>
#include <vector>
#include <map>
#include <queue>

using namespace std;


/// \brief Utility function: removes leading dashes from a string
void StripLeadingDashes( string& stringToModify );



/// \brief Class holding info about an individual command-line option/flag
class OptionObject
{
  public:
    // Constructors and Destructors:
    OptionObject( );
    virtual ~OptionObject( );
    
    // Public member functions:
    void DefineAsFlag( );   // specify that this is a "flag" options
    bool IsFlag( );    // is this a flag option?
    void SetFlag( );   // set flag value to true
    bool FlagSet( );   // has flag been set?
    void StoreTarget( const char targString[] );   // e.g., store filename pointed to by option
    bool TargetSet( );   // has target been set (if option is not flag)?
    string& GetTargetString( int n=0 );
    void EnableQueue( );
    bool IsQueue( );
    virtual int NTargetsStored( );

  protected:
    // Data members:
    bool  isFlag, flagSet, targetSet, isQueue;
    vector<string>  targetStrings;
    int  nStrings;
};



/// \brief Class for parsing command line
class CLineParser
{
  public:
    // Constructors and Destructors:
    CLineParser( );
    ~CLineParser( );
    
    // Public member functions:
    void PrintUsage( );
    void UnrecognizedAreErrors( );   //!< interpret unrecognized flags/options as errors
    void AddFlag( const string shortFlagString );
    void AddFlag( const string shortFlagString, const string longFlagString );
    void AddOption( const string shortOptString );
    void AddOption( const string shortOptString, const string longOptString );
    void AddQueueOption( const string shortOptString );
    void AddQueueOption( const string shortOptString, const string longOptString );
    void AddUsageLine( const string usageLine );
    int ParseCommandLine( int argc, char *argv[] );
    bool CommandLineEmpty( );
    bool FlagSet( const string flagName );
    bool OptionSet( const string optName );
    int GetNTargets( const string optName );
    string& GetTargetString( const string optName, int n=0 );
    int nArguments( );
    string& GetArgument( const int n );

  private:
    // Data members:
    map<string, OptionObject *>  optMap;   //!< 	data: map<string, *optObject> of option-objects
    vector<OptionObject *> optObjPointers;   //!< keep track of all discrete OptionObject pointers
    vector<string>  usageStrings;   //!< 	data: map<string, *optObject> of option-objects
    vector<string>  argStrings;     //!< 	data: vector<string> of argument strings
    int  verboseLevel;
    bool  commandLineEmpty;        //!< true if user supplied *no* options/flags *or* arguments
    bool  ignoreUnrecognized;      //!< if we encounter an unrecognized option/flag, do we
                                   //!< ignore it (= just print a warning & continue processing)?
    string  errorString1;
};


#endif  // _COMMANDLINE_PARSER_H_
