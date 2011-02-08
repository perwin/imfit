// Unit tests for code in commandline_parser.cpp

// $ cxxtestgen.py --error-printer -o test_runner.cpp unittest_commandline_parser.h
// $ g++ -Wno-write-strings -o test_runner test_runner.cpp commandline_parser.cpp -I/usr/local/include
//
// [the "-Wno-write-strings" is to suppress warnings when we create and use the
// argv c-string arrays]

// Things to test:
// [x]StripLeadingDashes() function
//
// [x]OptionObject class
//
// [x]CLineParser class


#include <cxxtest/TestSuite.h>

#include <string>
using namespace std;
#include "commandline_parser.h"


// global variables to test command-line inputs
const char *argv1[] = {"progName", "-b"};
const int argc1 = 2;



class NewTestSuite : public CxxTest::TestSuite 
{
public:

  // Tests for StripLeadingDashes function
  
  void testStripLeadingDashes( void )
  {
    string inputStr1, correctStr1;
    
    // correctly strip a single leading dash
    inputStr1 = "-b";
    correctStr1 = "b";
    StripLeadingDashes(inputStr1);
    TS_ASSERT( inputStr1 == correctStr1 );

    // correctly strip double leading dash
    inputStr1 = "--bob";
    correctStr1 = "bob";
    StripLeadingDashes(inputStr1);
    TS_ASSERT( inputStr1 == correctStr1 );

    // correctly handle a non-leading-dash string
    inputStr1 = "alpha-one";
    correctStr1 = "alpha-one";
    StripLeadingDashes(inputStr1);
    TS_ASSERT( inputStr1 == correctStr1 );
    
    // correctly handle strings consisting only of dashes
    inputStr1 = "-";
    correctStr1 = "";
    StripLeadingDashes(inputStr1);
    TS_ASSERT( inputStr1 == correctStr1 );
    inputStr1 = "--";
    correctStr1 = "";
    StripLeadingDashes(inputStr1);
    TS_ASSERT( inputStr1 == correctStr1 );

    // correctly handle empty string
    inputStr1 = "";
    correctStr1 = "";
    StripLeadingDashes(inputStr1);
    TS_ASSERT( inputStr1 == correctStr1 );
  }
  
  

  // Tests for OptionObject class
  
  void testOptionObject_as_flag( void )
  {
    OptionObject *testOptionObj;
    
    testOptionObj = new OptionObject;
    
    // default object is *not* a flag to start with
    TS_ASSERT( testOptionObj->IsFlag() == false );
    
    testOptionObj->DefineAsFlag();
    // should be a flag now
    TS_ASSERT( testOptionObj->IsFlag() == true );

    // flag should not initially be set
    TS_ASSERT( testOptionObj->FlagSet() == false );
    
    // set the flag and test that it has been set
    testOptionObj->SetFlag();
    TS_ASSERT( testOptionObj->FlagSet() == true );
  }


  void testOptionObject_as_option( void )
  {
    OptionObject *testOptionObj;
    string  testTargetString = "test_target";
    string  targetRecipient;
    
    testOptionObj = new OptionObject;
    
    // default object is *not* a flag
    TS_ASSERT( testOptionObj->IsFlag() == false );
    
    // target should initiallly *not* be set
    TS_ASSERT( testOptionObj->TargetSet() == false );

    // set the target string, then test for it
    testOptionObj->StoreTarget( testTargetString.c_str() );
    TS_ASSERT( testOptionObj->TargetSet() == true );
    
    // test to see if we stored the target string correctly
    TS_ASSERT( testOptionObj->TargetSet() == true );
    targetRecipient = testOptionObj->GetTargetString();
    TS_ASSERT( targetRecipient == testTargetString );
  }




  // Tests for CLineParser class
  
  void testCLineParser_setup( void )
  {
    CLineParser *testParser;
    
    testParser = new CLineParser();

    TS_ASSERT( testParser->IsCommandLineEmpty() == true );
  }


  // Test that we create single and double-string flags correctly
  void testCLineParser_CreateFlags( void )
  {
    CLineParser *testParser;
    OptionObject *optionObj;
    
    testParser = new CLineParser();

    testParser->AddFlag("b");
    TS_ASSERT( testParser->IsFlagSet("b") == false );
    testParser->AddFlag("z", "zeta");
    TS_ASSERT( testParser->IsFlagSet("z") == false );
    TS_ASSERT( testParser->IsFlagSet("zeta") == false );
  }


  // Test that we create single and double-string options correctly
  void testCLineParser_CreateOptions( void )
  {
    CLineParser *testParser;
    OptionObject *optionObj;
    
    testParser = new CLineParser();

    testParser->AddOption("x");
    TS_ASSERT( testParser->IsOptionSet("x") == false );
    testParser->AddOption("c", "config");
    TS_ASSERT( testParser->IsOptionSet("c") == false );
    TS_ASSERT( testParser->IsOptionSet("config") == false );
  }


  // Test that we correctly process a command line with no options/flags
  void testCLineParser_ParseEmptyLine( void )
  {
    int  argc = 1;
    char  *argv[] = {"progName"};
    int  status;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("b");
    testParser->AddOption("x");

    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT( status == 0 );
    TS_ASSERT( testParser->IsCommandLineEmpty() == true );
  }


  // Test that we correctly process a command line with mistakes
  void testCLineParser_ParseBad( void )
  {
    // this is a command line with an isolated double-dash
    int  argc1 = 2;
    char  *argv1[] = {"progName", "--"};
    // here, we specify a command line where an option is *not* followed by a target
    int  argc2 = 3;
    char  *argv2[] = {"progName", "-b", "-x"};
    int  status;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("b");
    testParser->AddOption("x");   // this means that "-x" requires a target

    status = testParser->ParseCommandLine(argc1, argv1);
    TS_ASSERT( status == -1 );

    status = testParser->ParseCommandLine(argc2, argv2);
    TS_ASSERT( status == -1 );
  }


  // Test that we correctly process a command line with one flag
  void testCLineParser_ParseSimpleFlag( void )
  {
    int  argc = 2;
    char  *argv[] = {"progName", "-b"};
    int  status;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddOption("x");

    TS_ASSERT( testParser->IsFlagSet("a") == false );
    TS_ASSERT( testParser->IsFlagSet("b") == false );
    TS_ASSERT( testParser->IsOptionSet("x") == false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT( status == 0 );
    TS_ASSERT( testParser->IsCommandLineEmpty() == false );
    TS_ASSERT( testParser->IsFlagSet("a") == false );
    TS_ASSERT( testParser->IsFlagSet("b") == true );
    TS_ASSERT( testParser->IsOptionSet("x") == false );
  }


  // Test that we correctly process a command line with one option
  void testCLineParser_ParseSimpleOption( void )
  {
    int  argc = 3;
    char  *argv[] = {"progName", "-x", "target_for_x"};
    int  status;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddOption("x");

    TS_ASSERT( testParser->IsFlagSet("a") == false );
    TS_ASSERT( testParser->IsFlagSet("b") == false );
    TS_ASSERT( testParser->IsOptionSet("x") == false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT( status == 0 );
    TS_ASSERT( testParser->IsCommandLineEmpty() == false );
    TS_ASSERT( testParser->IsFlagSet("a") == false );
    TS_ASSERT( testParser->IsFlagSet("b") == false );
    TS_ASSERT( testParser->IsOptionSet("x") == true );
    // check that we correctly stored the target
    targetString = testParser->GetTargetString("x");
    TS_ASSERT( targetString == "target_for_x" );
  }


  // Test that we correctly process a more complex command line
  void testCLineParser_ParseComplexLine( void )
  {
    int  argc = 4;
    char  *argv[] = {"progName", "-x", "target_for_x", "-b"};
    int  status;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddOption("x");

    TS_ASSERT( testParser->IsFlagSet("a") == false );
    TS_ASSERT( testParser->IsFlagSet("b") == false );
    TS_ASSERT( testParser->IsOptionSet("x") == false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT( status == 0 );
    TS_ASSERT( testParser->IsCommandLineEmpty() == false );
    TS_ASSERT( testParser->IsFlagSet("a") == false );
    TS_ASSERT( testParser->IsFlagSet("b") == true );
    TS_ASSERT( testParser->IsOptionSet("x") == true );
    // check that we correctly stored the target
    targetString = testParser->GetTargetString("x");
    TS_ASSERT( targetString == "target_for_x" );
  }


  // Test that we correctly extract arguments
  void testCLineParser_ParseArgument( void )
  {
    int  argc1 = 4;
    char  *argv1[] = {"progName", "-x", "target_for_x", "-b"};
    int  argc2 = 3;
    char  *argv2[] = {"progName", "alpha", "beta"};
    int  status;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddOption("x");

    // parse a commmand line with no arguments
    status = testParser->ParseCommandLine(argc1, argv1);
    TS_ASSERT( status == 0 );
    TS_ASSERT( testParser->nArguments() == 0 );
    
    // now check to see if we extract arguments
    status = testParser->ParseCommandLine(argc2, argv2);
    TS_ASSERT( status == 0 );
    // check that arguments were caught
    TS_ASSERT( testParser->nArguments() == 2 );
    TS_ASSERT( testParser->GetArgument(0) == "alpha" );
    TS_ASSERT( testParser->GetArgument(1) == "beta" );
  }


  // Test that we correctly extract flags, options, and arguments
  void testCLineParser_ParseComplexLine_with_args( void )
  {
    int  argc = 6;
    char  *argv[] = {"progName", "alpha", "-x", "target_for_x", "-b", "beta"};
    int  status;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddOption("x");

    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT( status == 0 );
    // check that flags and options were caught
    TS_ASSERT( testParser->IsFlagSet("a") == false );
    TS_ASSERT( testParser->IsFlagSet("b") == true );
    TS_ASSERT( testParser->IsOptionSet("x") == true );
    // check that we correctly stored the target
    targetString = testParser->GetTargetString("x");
    TS_ASSERT( targetString == "target_for_x" );
    // check that arguments were caught
    TS_ASSERT( testParser->nArguments() == 2 );
    TS_ASSERT( testParser->GetArgument(0) == "alpha" );
    TS_ASSERT( testParser->GetArgument(1) == "beta" );
  }

};
