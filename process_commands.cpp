/* 
 *
 *    Based on sample code (demo.cpp) provided with Kishan Thomas' anyoptions
 * library.
 */

#include "anyoption.h"


void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  AnyOption *opt = new AnyOption();
  opt->setVerbose(); /* print warnings about unknown options */
  opt->autoUsagePrint(true); /* print usage for bad options */

  /* SET THE USAGE/HELP   */
  opt->addUsage( "" );
  opt->addUsage( "Usage: " );
  opt->addUsage( "   imfit [options] imagefile.fits setupfile.dat" );
  opt->addUsage( " -h  --help  		Prints this help " );
  opt->addUsage( "     --noise noisemap.fits  Noise image " );
  opt->addUsage( "" );


  /* by default all options are checked on the command line and from option/resource file */
  opt->setFlag("help", 'h');
  opt->setOption("noise");      /* an option (takes an argument), supporting only long form */

  /* go through the command line and get the options  */
  opt->processCommandArgs( argc, argv );

  /* Process the results: options */
  if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) {
    opt->printUsage();
    exit(1);
  }
  if( opt->getValue( "title" ) != NULL )
    cout << "noise_image = " << opt->getValue("noise") << endl ;

  /* Process actual arguments, if any: */
  for (int i = 0 ; i < opt->getArgc() ; i++) {
    cout << "arg = " <<  opt->getArgv(i) << endl ;
  }


  delete opt;

}
