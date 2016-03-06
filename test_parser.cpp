/* Simple code for testing config_file_parser.cpp
 *
 */

#include <string>
#include <vector>
#include <stdio.h>
#include "mpfit.h"
#include "utilities_pub.h"
#include "config_file_parser.h"

using namespace std;


string  configFileName("test_config_sets.dat");
string  configFileName2("test_config_sets_and_extra.dat");


int main(int argc, char *argv[])
{
  vector<string>  functionList;
  vector<double>  parameterList;
  vector<mp_par>  paramLimits;
  vector<int>  setStarts;   // which function number marks start of new set
  bool  paramLimitsExist = false;
  bool mode2D_flag = true;   // we're testing the reading of 2D config files (non-profilefit)
  int  nSets, nParamsTot;
  int  status;
  configOptions  userConfigOptions, userConfigOptions2;
  
  printf("\nStarting ...\n\n");
  

  /* Read configuration file */
  if (! FileExists(configFileName.c_str())) {
    printf("\n*** WARNING: Unable to find or open configuration file \"%s\"!\n\n", 
           configFileName.c_str());
    return -1;
  }
  
  
  // Read file in limited mode (ignore parameter limits)
  printf("Starting simple (makeimage) mode:\n");
  status = ReadConfigFile(configFileName2, mode2D_flag, functionList, parameterList, 
                              setStarts, userConfigOptions);

  printf("User config options:\n");
  for (int j = 0; j < userConfigOptions.nOptions; j++) {
    printf("%s = %s\n", userConfigOptions.optionNames[j].c_str(),
    				userConfigOptions.optionValues[j].c_str());
  }
  printf("\n");

  printf("Number of sets = %d, number of functions = %d, number of parameters = %d\n", 
            (int)setStarts.size(), (int)functionList.size(), (int)parameterList.size());
  printf("\tSet starts (function numbers):");
  for (int i = 0; i < (int)setStarts.size(); i++)
    printf(" %d", setStarts[i]);
  printf("\n");
  nParamsTot = parameterList.size();
  printf("%d parameters found:\n", nParamsTot);
  for (int i = 0; i < nParamsTot; i++) {
    printf("\t#%d = %f\n", i, parameterList[i]);
  }
  printf("Done with simple (makeimage) mode...\n\n");

  
  printf("Starting full (imfit) mode...\n");
  functionList.clear();
  parameterList.clear();
  paramLimits.clear();
  setStarts.clear();
  // Read file in full mode (recognize parameter limits)
  status = ReadConfigFile(configFileName2, mode2D_flag, functionList, parameterList, paramLimits, 
                setStarts, paramLimitsExist, userConfigOptions2);
  
  printf("User config options:\n");
  for (int j = 0; j < userConfigOptions2.nOptions; j++) {
    printf("%s = %s\n", userConfigOptions2.optionNames[j].c_str(),
    				userConfigOptions2.optionValues[j].c_str());
  }
  printf("\n");
  printf("Number of functions = %d, number of parameters = %d\n", (int)functionList.size(),
            (int)parameterList.size());
  nSets = setStarts.size();
  printf("There are %d sets, with sets starting at\n", nSets);
  for (int i = 0; i < nSets; i++) {
    printf("\tfunction %d\n", setStarts[i]);
  }
  
  nParamsTot = parameterList.size();
  printf("%d parameters found:\n", nParamsTot);
  for (int i = 0; i < nParamsTot; i++) {
    printf("\t#%d = %f\n", i, parameterList[i]);
  }
  
  if (paramLimitsExist) {
    printf("\nParameter limits vector size = %d\n", (int)paramLimits.size());
    for (int i = 0; i < (int)paramLimits.size(); i++) {
      printf("#%d: fixed = %d, limited = %d,%d, limits = %f,%f\n", i, paramLimits[i].fixed,
                paramLimits[i].limited[0], paramLimits[i].limited[1],
                paramLimits[i].limits[0], paramLimits[i].limits[1]);
    }
  }
  printf("\nDone.\n");
}

