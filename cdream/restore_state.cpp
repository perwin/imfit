#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>

using namespace std;

#include "dream.h"
#include "array.h"

// [ ]TAG: lik
// int dream_restore_state( const dream_pars* p, Array3D<double>& state, vector<double *>& lik2,
//     					vector<double>& pCR, int& inBurnIn )
int dream_restore_state( const dream_pars* p, Array3D<double>& state, Array2D<double>& lik,
    					vector<double>& pCR, int& inBurnIn )
{
  int prevLines = 0;
  if (p->appendFile && (p->outputRootname != "" || p->outputRootname != "-")) {
    if (p->verboseLevel > 0)
      printf("Restoring previous state... ");
    prevLines = p->maxEvals;
    for (int i = 0; i < p->numChains; ++i) {
      if (p->verboseLevel > 0)
        printf("%d ", i);
      int line = 0;
      ostringstream chainFilename("");
      chainFilename << p->outputRootname << "." << i + 1 << ".txt";
      ifstream ifile(chainFilename.str().c_str());
      if (! ifile) {
        prevLines = -1;
        fprintf(stderr, "   pre-existing MCMC output file \"%s\" not found!", chainFilename.str().c_str());
        break;
      } else {
        string input;
        while (! getline(ifile,input).eof()) {
          if (line >= p->maxEvals) 
            break;
          if (input.length() == 0) 
            continue;
          if (input[0] == '#') 
            continue;
          istringstream istr(input);
          for (int j = 0; j < p->nvar; ++j) 
            istr >> state(line,i,j); 
          // [ ]TAG: lik
          // istr >> lik2[i][line] >> inBurnIn;
          istr >> lik(line,i) >> inBurnIn;
          for (int j(0); j < p->nCR; ++j) 
            istr >> pCR[j];
          ++line;
        }
        ifile.close();
        if (prevLines > line) 
          prevLines = line - 1;
      }
    }
  }
  return prevLines;
}


