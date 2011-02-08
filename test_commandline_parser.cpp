/* Simple code for testing new command-line parsing code
 *
 */

#include <string>
#include <vector>
#include <stdio.h>
#include "utilities_pub.h"
#include "anyoption.h"   // Kishan Thomas' class for command-line option parsing
#include "commandline_parser.h"

using namespace std;


/* ---------------- Definitions & Constants ----------------------------- */
#define MPFIT_SOLVER        1
#define DIFF_EVOLN_SOLVER   2

typedef struct {
  std::string  configFileName;
  std::string  imageFileName;
  bool  noImage;
  std::string  psfFileName;
  bool  psfImagePresent;
  std::string  noiseFileName;
  bool  noiseImagePresent;
  int  errorType;
  std::string  maskFileName;
  bool  maskImagePresent;
  int  maskFormat;
  bool  subsamplingFlag;
  std::string  outputModelFileName;
  bool  saveModel;
  bool  saveBestFitParams;
  std::string  outputParameterFileName;
  bool  useImageHeader;
  double  gain;
  bool  gainSet;
  double  readNoise;
  bool  readNoiseSet;
  int  nCombined;
  bool  nCombinedSet;
  double  originalSky;
  bool  originalSkySet;
  char  modelName[MAXLINE];
  bool  noModel;
  char  paramString[MAXLINE];
  bool  newParameters;
  double  magZeroPoint;
  bool  noParamLimits;
  bool  printImages;
  bool printChiSquaredOnly;
  int  solver;
} commandOptions;



/* ------------------- Function Prototypes ----------------------------- */
/* External functions: */

/* Local Functions: */
void PrintOptionStruct( commandOptions *theOptions );
void ProcessInputOld( int argc, char *argv[], commandOptions *theOptions );
void ProcessInputNew( int argc, char *argv[], commandOptions *theOptions );



int main(int argc, char *argv[])
{
  int  status;
  commandOptions  options;

  // Default settings for commandOptions struct:
  options.configFileName = "default_config_file";
  options.imageFileName = "default_image_file";
  options.psfFileName = "default_psf_file";
  options.noiseFileName = "default_noise_file";
  options.maskFileName = "default_mask_file";
  options.outputModelFileName = "default_output-image_file";
  options.outputParameterFileName = "default_output-params_file";
  options.noImage = true;
  options.psfImagePresent = false;
  options.noiseImagePresent = false;
  options.maskImagePresent = false;
  options.subsamplingFlag = true;
//   std::string  configFileName;
//   std::string  imageFileName;
//   bool  noImage;
//   std::string  psfFileName;
//   bool  psfImagePresent;
//   std::string  noiseFileName;
//   bool  noiseImagePresent;
//   int  errorType;
//   std::string  maskFileName;
//   bool  maskImagePresent;
//   int  maskFormat;
//   bool  subsamplingFlag;
//   std::string  outputModelFileName;
//   bool  saveModel;
//   bool  saveBestFitParams;
//   std::string  outputParameterFileName;
//   bool  useImageHeader;
//   double  gain;
//   bool  gainSet;
//   double  readNoise;
//   bool  readNoiseSet;
//   int  nCombined;
//   bool  nCombinedSet;
//   double  originalSky;
//   bool  originalSkySet;
//   char  modelName[MAXLINE];
//   bool  noModel;
//   char  paramString[MAXLINE];
//   bool  newParameters;
//   double  magZeroPoint;
//   bool  noParamLimits;
//   bool  printImages;
//   bool printChiSquaredOnly;
//   int  solver;
// } commandOptions;

  printf("\nStarting ...\n\n");
  
//  printf("Processing command-line using Kishan Thomas's anyoption.cpp...\n");
//  ProcessInputOld(argc, argv, &options);

  printf("\n\nProcessing command-line using new code [NOTHING YET!]...\n");
  ProcessInputNew(argc, argv, &options);


  PrintOptionStruct(&options);
  
  printf("\nDone.\n");
}



void PrintOptionStruct( commandOptions *theOptions )
{
  printf("\nHere are the current settings of the commandOptions structure:\n");
  printf("\toptions->configFileName = %s\n", theOptions->configFileName.c_str());
  printf("\toptions->imageFileName = %s\n", theOptions->imageFileName.c_str());
  printf("\toptions->psfFileName = %s\n", theOptions->psfFileName.c_str());
  printf("\toptions->noiseFileName = %s\n", theOptions->noiseFileName.c_str());
  printf("\toptions->maskFileName = %s\n", theOptions->maskFileName.c_str());
  printf("\toptions->outputModelFileName = %s\n", theOptions->outputModelFileName.c_str());
  printf("\toptions->outputParameterFileName = %s\n", theOptions->outputParameterFileName.c_str());
  
  printf("\toptions->noImage = %d\n", theOptions->noImage);
  printf("\toptions->psfImagePresent = %d\n", theOptions->psfImagePresent);
  printf("\toptions->noiseImagePresent = %d\n", theOptions->noiseImagePresent);
  printf("\toptions->maskImagePresent = %d\n", theOptions->maskImagePresent);
  printf("\toptions->subsamplingFlag = %d\n", theOptions->subsamplingFlag);
}



void ProcessInputNew( int argc, char *argv[], commandOptions *theOptions )
{
  CLineParser  *optParser = new CLineParser();
  
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   imfit [options] imagefile.fits");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
  optParser->AddUsageLine("     --list-functions         Prints list of available functions (components)");
  optParser->AddUsageLine("     --chisquare-only         Print chi^2 of input model and quit");
  optParser->AddUsageLine(" -c  --config <config-file>   configuration file");
  optParser->AddUsageLine("     --de                     Use differential evolution solver instead of L-M");
  optParser->AddUsageLine("     --noise <noisemap.fits>  Noise image");
  optParser->AddUsageLine("     --mask <mask.fits>       Mask image");
  optParser->AddUsageLine("     --psf <psf.fits>         PSF image");
  optParser->AddUsageLine("     --nosubsampling          Do *not* do pixel subsampling near centers");
  optParser->AddUsageLine("     --save-params <output-file>          Save best-fit parameters in config-file format");
  optParser->AddUsageLine("     --save-model <outputname.fits>       Save best-fit model image");
  optParser->AddUsageLine("     --use-headers            Use image header values for gain, readnoise [NOT IMPLEMENTED YET]");
  optParser->AddUsageLine("     --sky <sky-level>        Original sky background (ADUs) which was subtracted from image");
  optParser->AddUsageLine("     --gain <value>           Image gain (e-/ADU)");
  optParser->AddUsageLine("     --readnoise <value>      Image read noise (e-)");
  optParser->AddUsageLine("     --ncombined <value>      Number of images averaged to make final image (if counts are average or median)");
  optParser->AddUsageLine("     --errors-are-variances   Indicates that values in noise image = variances (instead of sigmas)");
  optParser->AddUsageLine("     --errors-are-weights     Indicates that values in noise image = weights (instead of sigmas)");
  optParser->AddUsageLine("     --mask-zero-is-bad       Indicates that zero values in mask = *bad* pixels");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --printimage             Print out images (for debugging)");
  optParser->AddUsageLine("");

  optParser->AddFlag("h", "help");
  optParser->AddFlag("de");
  optParser->AddFlag("use-headers");
  optParser->AddFlag("nosubsampling");
  optParser->AddOption("c");
  optParser->AddOption("mask");
//  optParser->PrintUsage();

  printf("Is \"h\" set? %d\n", optParser->IsFlagSet("h"));
  printf("Is \"help\" set? %d\n", optParser->IsFlagSet("help"));
  printf("Is \"de\" set? %d\n\n", optParser->IsFlagSet("de"));
  
  int status = optParser->ParseCommandLine(argc, argv);
  if (status < 0) {
    fprintf(stderr, "\nFailure in command-line processing; quitting...\n\n");
    return;
  }

  printf("Is \"h\" set? %d\n", optParser->IsFlagSet("h"));
  printf("Is \"help\" set? %d\n", optParser->IsFlagSet("help"));
  printf("Is \"de\" set? %d\n", optParser->IsFlagSet("de"));
  printf("Does \"mask\" have a target? %d\n", optParser->IsOptionSet("mask"));
  if (optParser->IsOptionSet("mask"))
    printf("\ttarget for \"mask\" = \"%s\"\n", optParser->GetTargetString("mask").c_str());
  printf("\n");
  
  int nArgs = optParser->nArguments();
  printf("%d arguments found:\n", nArgs);
  for (int i = 0; i < nArgs; i++) {
    printf("   %s", optParser->GetArgument(i).c_str());
  }
  
  
  printf("\n");
  
  delete optParser;
}



void ProcessInputOld( int argc, char *argv[], commandOptions *theOptions )
{

  AnyOption *opt = new AnyOption();
  opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* SET THE USAGE/HELP   */
  opt->addUsage("Usage: ");
  opt->addUsage("   imfit [options] imagefile.fits");
  opt->addUsage(" -h  --help                   Prints this help");
  opt->addUsage("     --list-functions         Prints list of available functions (components)");
  opt->addUsage("     --chisquare-only         Print chi^2 of input model and quit");
  opt->addUsage(" -c  --config <config-file>   configuration file");
  opt->addUsage("     --de                     Use differential evolution solver instead of L-M");
  opt->addUsage("     --noise <noisemap.fits>  Noise image");
  opt->addUsage("     --mask <mask.fits>       Mask image");
  opt->addUsage("     --psf <psf.fits>         PSF image");
  opt->addUsage("     --nosubsampling          Do *not* do pixel subsampling near centers");
  opt->addUsage("     --save-params <output-file>          Save best-fit parameters in config-file format");
  opt->addUsage("     --save-model <outputname.fits>       Save best-fit model image");
  opt->addUsage("     --use-headers            Use image header values for gain, readnoise [NOT IMPLEMENTED YET]");
  opt->addUsage("     --sky <sky-level>        Original sky background (ADUs) which was subtracted from image");
  opt->addUsage("     --gain <value>           Image gain (e-/ADU)");
  opt->addUsage("     --readnoise <value>      Image read noise (e-)");
  opt->addUsage("     --ncombined <value>      Number of images averaged to make final image (if counts are average or median)");
  opt->addUsage("     --errors-are-variances   Indicates that values in noise image = variances (instead of sigmas)");
  opt->addUsage("     --errors-are-weights     Indicates that values in noise image = weights (instead of sigmas)");
  opt->addUsage("     --mask-zero-is-bad       Indicates that zero values in mask = *bad* pixels");
  opt->addUsage("");
  opt->addUsage("     --printimage             Print out images (for debugging)");
  opt->addUsage("");


  /* by default all options are checked on the command line and from option/resource file */
  opt->setFlag('q');
  opt->setFlag("help", 'h');
  opt->setFlag("list-functions");
  opt->setFlag("printimage");
  opt->setFlag("chisquare-only");
  opt->setFlag("de");
  opt->setFlag("use-headers");
  opt->setFlag("errors-are-variances");
  opt->setFlag("errors-are-weights");
  opt->setFlag("mask-zero-is-bad");
  opt->setFlag("nosubsampling");
  opt->setOption("noise");      /* an option (takes an argument), supporting only long form */
  opt->setOption("mask");      /* an option (takes an argument), supporting only long form */
  opt->setOption("psf");      /* an option (takes an argument), supporting only long form */
  opt->setOption("save-params");      /* an option (takes an argument), supporting only long form */
  opt->setOption("save-model");      /* an option (takes an argument), supporting only long form */
  opt->setOption("sky");        /* an option (takes an argument), supporting only long form */
  opt->setOption("gain");        /* an option (takes an argument), supporting only long form */
  opt->setOption("readnoise");        /* an option (takes an argument), supporting only long form */
  opt->setOption("ncombined");        /* an option (takes an argument), supporting only long form */
  opt->setOption("config", 'c');

  /* parse the command line:  */
  opt->processCommandArgs( argc, argv );


  /* Process the results: actual arguments, if any: */
  if (opt->getArgc() > 0) {
    theOptions->imageFileName = opt->getArgv(0);
    theOptions->noImage = false;
    printf("\tImage file = %s\n", theOptions->imageFileName.c_str());
  }

  /* Process the results: options */
  // First two are options which print useful info and then exit the program
  if ( opt->getFlag('q') ) {
  	printf("\n\nimfit_main: q option found!\n\n");
    delete opt;
    exit(1);
  }
  if ( opt->getFlag("help") || opt->getFlag('h') || (! opt->hasOptions()) ) {
    opt->printUsage();
    delete opt;
    exit(1);
  }
  if (opt->getFlag("list-functions")) {
    printf("\nHere we print the available functions ...\n");
    delete opt;
    exit(1);
  }

  if (opt->getFlag("printimage")) {
    theOptions->printImages = true;
  }
  if (opt->getFlag("chisquare-only")) {
    printf("\t No fitting will be done!\n");
    theOptions->printChiSquaredOnly = true;
  }
  if (opt->getFlag("de")) {
  	printf("\t Differential Evolution selected!\n");
  	theOptions->solver = DIFF_EVOLN_SOLVER;
  }
  if (opt->getFlag("nosubsampling")) {
    theOptions->subsamplingFlag = false;
  }
  if (opt->getFlag("use-header")) {
    theOptions->useImageHeader = true;
  }
  if (opt->getFlag("errors-are-variances")) {
    theOptions->errorType = WEIGHTS_ARE_VARIANCES;
  }
  if (opt->getFlag("errors-are-weights")) {
    theOptions->errorType = WEIGHTS_ARE_WEIGHTS;
  }
  if (opt->getFlag("mask-zero-is-bad")) {
    theOptions->maskFormat = MASK_ZERO_IS_BAD;
  }
  if (opt->getValue("config") != NULL) {
    theOptions->configFileName = opt->getValue("config");
    printf("\tconfiguration file = %s\n", theOptions->configFileName.c_str());
  }
  if (opt->getValue("noise") != NULL) {
    theOptions->noiseFileName = opt->getValue("noise");
    theOptions->noiseImagePresent = true;
    printf("\tnoise image = %s\n", theOptions->noiseFileName.c_str());
  }
  if (opt->getValue("psf") != NULL) {
    theOptions->psfFileName = opt->getValue("psf");
    theOptions->psfImagePresent = true;
    printf("\tPSF image = %s\n", theOptions->psfFileName.c_str());
  }
  if (opt->getValue("mask") != NULL) {
    theOptions->maskFileName = opt->getValue("mask");
    theOptions->maskImagePresent = true;
    printf("\tmask image = %s\n", theOptions->maskFileName.c_str());
  }
  if (opt->getValue("save-model") != NULL) {
    theOptions->outputModelFileName = opt->getValue("save-model");
    theOptions->saveModel = true;
    printf("\toutput best-fit model image = %s\n", theOptions->outputModelFileName.c_str());
  }
  if (opt->getValue("save-params") != NULL) {
    theOptions->outputParameterFileName = opt->getValue("save-params");
    theOptions->saveBestFitParams = true;
    printf("\toutput best-fit parameter file = %s\n", theOptions->outputParameterFileName.c_str());
  }
  if (opt->getValue("sky") != NULL) {
    if (NotANumber(opt->getValue("sky"), 0, kAnyReal)) {
      fprintf(stderr, "*** WARNING: sky should be a real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->originalSky = atof(opt->getValue("sky"));
    theOptions->originalSkySet = true;
    printf("\toriginal sky level = %g ADU\n", theOptions->originalSky);
  }
  if (opt->getValue("gain") != NULL) {
    if (NotANumber(opt->getValue("gain"), 0, kPosReal)) {
      fprintf(stderr, "*** WARNING: gain should be a positive real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->gain = atof(opt->getValue("gain"));
    theOptions->gainSet = true;
    printf("\tgain = %g e-/ADU\n", theOptions->gain);
  }
  if (opt->getValue("readnoise") != NULL) {
    if (NotANumber(opt->getValue("readnoise"), 0, kPosReal)) {
      fprintf(stderr, "*** WARNING: read noise should be a non-negative real number!\n");
      delete opt;
      exit(1);
    }
    theOptions->readNoise = atof(opt->getValue("readnoise"));
    theOptions->readNoiseSet = true;
    printf("\tread noise = %g e-\n", theOptions->readNoise);
  }
  if (opt->getValue("ncombined") != NULL) {
    if (NotANumber(opt->getValue("ncombined"), 0, kPosInt)) {
      fprintf(stderr, "*** WARNING: ncombined should be a positive integer!\n");
      delete opt;
      exit(1);
    }
    theOptions->nCombined = atoi(opt->getValue("ncombined"));
    theOptions->nCombinedSet = true;
    printf("\tn_combined = %d\n", theOptions->nCombined);
  }

  delete opt;

}

