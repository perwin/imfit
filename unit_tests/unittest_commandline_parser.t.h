// Unit tests for code in commandline_parser.cpp

// See run_unittest_cmlineparser.sh for how to compile and run these tests.

// older compilation notes:
// $ cxxtestgen.py --error-printer -o test_runner.cpp unittest_commandline_parser.h
// $ g++ -Wno-write-strings -o test_runner test_runner.cpp commandline_parser.cpp utilities.cpp -I/usr/local/include
//
// [the "-Wno-write-strings" is to suppress warnings when we create and use the
// argv c-string arrays]



#include <cxxtest/TestSuite.h>

#include <string>
using namespace std;
#include "commandline_parser.h"


// global variables to test command-line inputs
const char *argv1[] = {"progName", "-b"};
const int argc1 = 2;



class TestStripLeadingDashes : public CxxTest::TestSuite 
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
    TS_ASSERT_EQUALS( inputStr1, correctStr1 );

    // correctly strip double leading dash
    inputStr1 = "--bob";
    correctStr1 = "bob";
    StripLeadingDashes(inputStr1);
    TS_ASSERT_EQUALS( inputStr1, correctStr1 );

    // correctly handle a non-leading-dash string
    inputStr1 = "alpha-one";
    correctStr1 = "alpha-one";
    StripLeadingDashes(inputStr1);
    TS_ASSERT_EQUALS( inputStr1, correctStr1 );
    
    // correctly handle strings consisting only of dashes
    inputStr1 = "-";
    correctStr1 = "";
    StripLeadingDashes(inputStr1);
    TS_ASSERT_EQUALS( inputStr1, correctStr1 );
    inputStr1 = "--";
    correctStr1 = "";
    StripLeadingDashes(inputStr1);
    TS_ASSERT_EQUALS( inputStr1, correctStr1 );

    // correctly handle empty string
    inputStr1 = "";
    correctStr1 = "";
    StripLeadingDashes(inputStr1);
    TS_ASSERT_EQUALS( inputStr1, correctStr1 );
  }
};
  


class TestOptionObject : public CxxTest::TestSuite 
{
public:

  // Tests for OptionObject class
  
  void testOptionObject_as_flag( void )
  {
    OptionObject *testOptionObj;
    
    testOptionObj = new OptionObject();
    
    // default object is *not* a flag to start with
    TS_ASSERT_EQUALS( testOptionObj->IsFlag(), false );
    
    testOptionObj->DefineAsFlag();
    // should be a flag now
    TS_ASSERT_EQUALS( testOptionObj->IsFlag(), true );

    // flag should not initially be set
    TS_ASSERT_EQUALS( testOptionObj->FlagSet(), false );
    
    // set the flag and test that it has been set
    testOptionObj->SetFlag();
    TS_ASSERT_EQUALS( testOptionObj->FlagSet(), true );
  }


  void testOptionObject_as_option( void )
  {
    OptionObject *testOptionObj;
    string  testTargetString = "test_target";
    string  targetRecipient;
    
    testOptionObj = new OptionObject();
    
    // default object is *not* a flag
    TS_ASSERT_EQUALS( testOptionObj->IsFlag(), false );
    
    // target should initiallly *not* be set
    TS_ASSERT_EQUALS( testOptionObj->TargetSet(), false );

    // set the target string, then test for it
    testOptionObj->StoreTarget( testTargetString.c_str() );
    TS_ASSERT_EQUALS( testOptionObj->TargetSet(), true );
    
    // test to see if we stored the target string correctly
    TS_ASSERT_EQUALS( testOptionObj->TargetSet(), true );
    targetRecipient = testOptionObj->GetTargetString();
    TS_ASSERT_EQUALS( targetRecipient, testTargetString );
  }
};




class TestQeueOptionObject : public CxxTest::TestSuite 
{
public:

  void testQueuOptionObject( void )
  {
    OptionObject *testOptionObj;
    string  testTargetString1 = "test_target1";
    string  testTargetString2 = "test_target2";
    string  targetRecipient;
    
    testOptionObj = new OptionObject();
    testOptionObj->EnableQueue();

    TS_ASSERT_EQUALS( testOptionObj->IsQueue(), true );

    // default object is *not* a flag
    TS_ASSERT_EQUALS( testOptionObj->IsFlag(), false );
    
    // target should initiallly *not* be set
    TS_ASSERT_EQUALS( testOptionObj->TargetSet(), false );
    TS_ASSERT_EQUALS( testOptionObj->NTargetsStored(), 0 );

    // set the target string, then test for it
    testOptionObj->StoreTarget( testTargetString1.c_str() );
    TS_ASSERT_EQUALS( testOptionObj->TargetSet(), true );
    TS_ASSERT_EQUALS( testOptionObj->NTargetsStored(), 1 );
    
    // store another instance
    testOptionObj->StoreTarget( testTargetString2.c_str() );
    TS_ASSERT_EQUALS( testOptionObj->NTargetsStored(), 2 );

    // test to see if we stored the target strings correctly, and if
    // retrieving them decrements the count
    targetRecipient = testOptionObj->GetTargetString(0);
    TS_ASSERT_EQUALS( targetRecipient, testTargetString1 );
    targetRecipient = testOptionObj->GetTargetString(1);
    TS_ASSERT_EQUALS( targetRecipient, testTargetString2 );
  }
};




class TestCLineParser : public CxxTest::TestSuite 
{
public:

  // Tests for CLineParser class
  
  void testCLineParser_setup( void )
  {
    CLineParser *testParser;
    
    testParser = new CLineParser();

    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), true );
  }


  // Test that we create single and double-string flags correctly
  void testCLineParser_CreateFlags( void )
  {
    CLineParser *testParser;
    
    testParser = new CLineParser();

    testParser->AddFlag("b");
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    testParser->AddFlag("z", "zeta");
    TS_ASSERT_EQUALS( testParser->FlagSet("z"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("zeta"), false );
  }


  // Test that we create single and double-string options correctly
  void testCLineParser_CreateOptions( void )
  {
    CLineParser *testParser;
    
    testParser = new CLineParser();

    testParser->AddOption("x");
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    testParser->AddOption("c", "config");
    TS_ASSERT_EQUALS( testParser->OptionSet("c"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("config"), false );
  }


  // Test that we create single and double-string multi-instance ("queue") options correctly
  void testCLineParser_CreateQueueOptions( void )
  {
    CLineParser *testParser;
    
    testParser = new CLineParser();

    testParser->AddQueueOption("x");
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    testParser->AddQueueOption("c", "config");
    TS_ASSERT_EQUALS( testParser->OptionSet("c"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("config"), false );
  }


  // Test that we catch erroneous flag names
  void testCLineParser_CheckBadNames( void )
  {
    CLineParser *testParser;
    
    testParser = new CLineParser();

    testParser->AddFlag("x");
    TS_ASSERT_EQUALS( testParser->FlagSet("x"), false );
    // OK, now try it with an incorrect name
    TS_ASSERT_EQUALS( testParser->FlagSet("qq"), false );
    

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
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), true );
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
    TS_ASSERT_EQUALS( status, -1 );

    status = testParser->ParseCommandLine(argc2, argv2);
    TS_ASSERT_EQUALS( status, -1 );
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

    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
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

    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), true );
    // check that we correctly stored the target
    targetString = testParser->GetTargetString("x");
    TS_ASSERT_EQUALS( targetString, "target_for_x" );
  }


  // Test that we correctly process a command line with one option, using "="
  void testCLineParser_ParseSimpleOption_with_equals( void )
  {
    int  argc = 2;
    char  *argv[] = {"progName", "-x=target_for_x"};
    int  status;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddOption("x");

    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), true );
    // check that we correctly stored the target
    targetString = testParser->GetTargetString("x");
    TS_ASSERT_EQUALS( targetString, "target_for_x" );
  }


  // Test that we correctly process a command line with one queue option
  // and multiple invocations thereof
  void testCLineParser_ParseQueueOption( void )
  {
    int  argc = 5;
    char  *argv[] = {"progName", "-x", "target_for_x_1", "-x", "target_for_x_2"};
    int  status, nTargets;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddQueueOption("x");

    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), true );
    
    // check that we correctly stored the targets
    nTargets = testParser->GetNTargets("x");
    TS_ASSERT_EQUALS( nTargets, 2 );
    targetString = testParser->GetTargetString("x", 0);
    TS_ASSERT_EQUALS( targetString, "target_for_x_1" );
    targetString = testParser->GetTargetString("x", 1);
    TS_ASSERT_EQUALS( targetString, "target_for_x_2" );
  }


  // Test that we correctly process a command line with one queue option
  // and multiple invocations thereof
  void testCLineParser_ParseQueueOption_complex( void )
  {
    int  argc = 15;
    char  *argv[] = {"progName", "-x", "target_for_x_1", "-x", "target_for_x_2",
    				"-y", "target_for_y_1", "-y", "target_for_y_2",
    				"-z", "target_for_z_1", "-z", "target_for_z_2", "-z", "target_for_z_3"};
    int  status, nTargets;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddFlag("a");
    testParser->AddFlag("b");
    testParser->AddQueueOption("x");
    testParser->AddQueueOption("y");
    testParser->AddQueueOption("z");

    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("y"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("z"), false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("y"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("z"), true );
    
    // check that we correctly stored the targets
    nTargets = testParser->GetNTargets("x");
    TS_ASSERT_EQUALS( nTargets, 2 );
    targetString = testParser->GetTargetString("x", 0);
    TS_ASSERT_EQUALS( targetString, "target_for_x_1" );
    targetString = testParser->GetTargetString("x", 1);
    TS_ASSERT_EQUALS( targetString, "target_for_x_2" );

    nTargets = testParser->GetNTargets("y");
    TS_ASSERT_EQUALS( nTargets, 2 );
    targetString = testParser->GetTargetString("y", 0);
    TS_ASSERT_EQUALS( targetString, "target_for_y_1" );
    targetString = testParser->GetTargetString("y", 1);
    TS_ASSERT_EQUALS( targetString, "target_for_y_2" );

    nTargets = testParser->GetNTargets("z");
    TS_ASSERT_EQUALS( nTargets, 3 );
    targetString = testParser->GetTargetString("z", 0);
    TS_ASSERT_EQUALS( targetString, "target_for_z_1" );
    targetString = testParser->GetTargetString("z", 1);
    TS_ASSERT_EQUALS( targetString, "target_for_z_2" );
    targetString = testParser->GetTargetString("z", 2);
    TS_ASSERT_EQUALS( targetString, "target_for_z_3" );
  }


  // Test that we correctly process a command line with one queue option
  // and multiple invocations thereof



  void testCLineParser_ParseQueueOption_complex2( void )
  {
    int  argc = 14;
    char  *argv[] = {"./makeimage", "tests/config_makeimage_gauss-oversample.dat",
    				 "-o", "temptest/oversampled.fits", "--psf", "tests/psf_standard.fits",
    				 "--overpsf", "tests/psf_oversamp.fits", "--overpsf_scale", "3",
    				"--overpsf_region", "100:110,100:110", "--overpsf_region", "50:60,80:90"};
    int  status;
    string  targetString;
    CLineParser  *testParser;
    
    testParser = new CLineParser();
    testParser->AddOption("o");
    testParser->AddOption("psf");
    testParser->AddQueueOption("overpsf");
    testParser->AddQueueOption("overpsf_scale");
    testParser->AddQueueOption("overpsf_region");

    TS_ASSERT_EQUALS( testParser->OptionSet("o"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("psf"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("overpsf"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("overpsf_scale"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("overpsf_region"), false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("o"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("psf"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("overpsf"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("overpsf_scale"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("overpsf_region"), true );
    
    // check that we correctly stored the targets
    targetString = testParser->GetTargetString("o");
    TS_ASSERT_EQUALS( targetString, "temptest/oversampled.fits" );

    targetString = testParser->GetTargetString("psf");
    TS_ASSERT_EQUALS( targetString, "tests/psf_standard.fits" );

    targetString = testParser->GetTargetString("overpsf");
    int nTargets = testParser->GetNTargets("overpsf");
    TS_ASSERT_EQUALS( nTargets, 1 );
    TS_ASSERT_EQUALS( targetString, "tests/psf_oversamp.fits" );

    targetString = testParser->GetTargetString("overpsf", 0);
    TS_ASSERT_EQUALS( targetString, "tests/psf_oversamp.fits" );

    targetString = testParser->GetTargetString("overpsf_scale");
    TS_ASSERT_EQUALS( targetString, "3" );

    nTargets = testParser->GetNTargets("overpsf_region");
    TS_ASSERT_EQUALS( nTargets, 2 );
    targetString = testParser->GetTargetString("overpsf_region", 0);
    TS_ASSERT_EQUALS( targetString, "100:110,100:110" );
    targetString = testParser->GetTargetString("overpsf_region", 1);
    TS_ASSERT_EQUALS( targetString, "50:60,80:90" );
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

    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), false );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), false );
    status = testParser->ParseCommandLine(argc, argv);
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->CommandLineEmpty(), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), true );
    // check that we correctly stored the target
    targetString = testParser->GetTargetString("x");
    TS_ASSERT_EQUALS( targetString, "target_for_x" );
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
    TS_ASSERT_EQUALS( status, 0 );
    TS_ASSERT_EQUALS( testParser->nArguments(), 0 );
    
    // now check to see if we extract arguments
    status = testParser->ParseCommandLine(argc2, argv2);
    TS_ASSERT_EQUALS( status, 0 );
    // check that arguments were caught
    TS_ASSERT_EQUALS( testParser->nArguments(), 2 );
    TS_ASSERT_EQUALS( testParser->GetArgument(0), "alpha" );
    TS_ASSERT_EQUALS( testParser->GetArgument(1), "beta" );
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
    TS_ASSERT_EQUALS( status, 0 );
    // check that flags and options were caught
    TS_ASSERT_EQUALS( testParser->FlagSet("a"), false );
    TS_ASSERT_EQUALS( testParser->FlagSet("b"), true );
    TS_ASSERT_EQUALS( testParser->OptionSet("x"), true );
    // check that we correctly stored the target
    targetString = testParser->GetTargetString("x");
    TS_ASSERT_EQUALS( targetString, "target_for_x" );
    // check that arguments were caught
    TS_ASSERT_EQUALS( testParser->nArguments(), 2 );
    TS_ASSERT_EQUALS( testParser->GetArgument(0), "alpha" );
    TS_ASSERT_EQUALS( testParser->GetArgument(1), "beta" );
  }

};
