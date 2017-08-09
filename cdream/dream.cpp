// ===========================================================================
//
// DREAM 
// -----
//
// Gabriel E Leventhal
// Institute of Integrative Biology
// ETH Zurich
// Universitätstrasse 16
// 8092 Zürich
// Switzerland
//
// gabriel.leventhal@env.ethz.ch
// http://www.leventhal.ch
//
// DREAM algorithm:
//
// Vrugt, J. A., ter Braak, C. J. F., Diks, C. G. H., Robinson, B. A., Hyman, 
// J. M., Higdon, D., 2009. Accelerating Markov chain Monte Carlo simulation 
// by differential evolution with self-adaptive randomized subspace sampling. 
// International Journal of Nonlinear Sciences and Numerical Simulation 
// 10 (3), 273-290. DOI: 10.1515/IJNSNS.2009.10.3.273
//
// ===========================================================================

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <streambuf>
#include <vector>
#include <algorithm>
#include <map>
#include <stdio.h>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>

#include <rng/RngStream.h>

#include "dream.h"
#include "dream_params.h"

#include "model_object.h"
#include "utilities_pub.h"


int dream( const dream_pars* p, rng::RngStream* rng )
{
  int inBurnIn = (p->burnIn > 0);
  int burnInStart = 0;
  int genNumber = 0;
  int nLikelihoodEvals = 0;
  bool converged = false;
  
  ModelObject * theModel = (ModelObject *)p->extraData;

  // opening message
  if (p->verboseLevel > 0) {
    printf("Beginning DREAM:\n");
    for (int i = 0; i < p->nvar; ++i) {
      printf("%16s = [%8.1f, %8.1f]", p->parameterNames[i].c_str(), p->varLo[i], p->varHi[i]);
      if (p->varLock[i]) 
        printf(" -- FIXED");
      printf("\n");
    }

    printf("Collapse outliers: %d\n", p->collapseOutliers);
  }
  
  // MCMC chains
  Array3D<double> state(p->maxEvals, p->numChains, p->nvar);
  // [x]TAG: lik
  vector<double *> lik2(p->numChains);
  for (int i = 0; i < p->numChains; i++) {
    double * tempArrayPtr = new double[p->maxEvals + 1];
    lik2.push_back(tempArrayPtr);
  }
  Array2D<double> lik(p->maxEvals + 1, p->numChains);

  Array2D<double> proposal(p->numChains, p->nvar);
  Array2D<double> proposal_two(p->numChains, p->nvar);
  proposal.set_all(0.0);
  proposal_two.set_all(0.0);

  vector<double> scaleReduction(p->nvar, 0.0);
  vector<double> pCR(p->nCR, 1.0/p->nCR);
  
  // PE: needed for easier printing of individual steps
  double *tempParams = (double *)malloc(p->nvar * sizeof(double));

  // =========================================================================
  // read previous state

  // [ ]TAG: lik
  int prevLines = dream_restore_state(p, state, lik, pCR, inBurnIn);
  if (prevLines < 0) {
    fprintf(stderr, "DREAM: previous MCMC output files not found!\n");
    return DREAM_EXIT_NO_OUTPUT_FILES;
  }
  // =========================================================================
  // open output file

  ostringstream chainFilename;
  vector<ostream*> oout;
  ios_base::openmode fmode = (p->appendFile) ? (ios_base::out | ios_base::app) : ios_base::out;

  // PE: changed output chain-file names so they start with 1, not 0
  oout.resize(p->numChains, &cout);
  if (p->outputRootname != "" && p->outputRootname != "-") {
    for (int i = 0; i < p->numChains; ++i) {
      chainFilename.str("");
      chainFilename << p->outputRootname << "." << i + 1 << ".txt";
      oout[i] = new ofstream(chainFilename.str().c_str(), fmode);
      oout[i]->setf(ios::scientific, ios::floatfield);
      oout[i]->precision(12);
      for (int n = 0; n < p->outputHeaderLines.size(); n++)
        *oout[i] << p->outputHeaderLines[n];
      if (p->appendFile)
        *oout[i] << "# --- Resuming DREAM ---" << endl;
    }
  }

  // =========================================================================
  // Initialize with latin hypercube sampling if not a resumed run

  int do_calc(1);

  if (! p->appendFile) {
    Array2DView<double> initVar(state.n_y(), state.n_z(), state.pt(0,0,0));
    // [ ]TAG: lik
    ArrayView<double> initLik(lik.n_y(), lik.pt(0,0));
    dream_initialize(p, rng, initVar, initLik);

    // save initial state of each chain
    for (int i = 0; i < p->numChains; ++i) {
      for (int j = 0; j < p->nvar; ++j) 
        tempParams[j] = state(0,i,j);
      string paramString = theModel->PrintModelParamsHorizontalString(tempParams);
      *oout[i] << paramString << " ";
      // [.]TAG: lik
      // *oout[i]  << lik2[0][i] << " " << inBurnIn << " " << " ";
      *oout[i]  << lik(0,i) << " " << inBurnIn << " " << " ";
      for (int j = 0; j < p->nCR; ++j) 
        *oout[i] << pCR[j] << " ";
      *oout[i] << 1 << endl;
    }
  }

  // =========================================================================

  double gamma;
  int r1, r2;
  double crossRate = 1.0;
  int curRun = 0;
  int gammaGeneration = 0;
  int delta;
  int numAccepted = 0;
  int ireport = 0;
  double drand = 0.0;

  vector<int> acceptStep(p->numChains, 0);

  vector<double> pairDiff(p->nvar);
  vector<double> e(p->nvar);
  vector<double> epsilon(p->nvar);
  vector<bool> updatePar(p->nvar, false);
  vector<int>  updateDim(p->numChains, p->nfree);
  vector<double> step(p->nvar);

  vector<double> bstar(p->nvar, p->bstar_zero);
//  bstar[4] = 1;   // PE: WTF?
 
  // ======================================================================
  // RUN MCMC

  vector<unsigned> L(p->nCR, 0);      // candidates for crossover
  vector<unsigned> totalSteps(p->nCR, 0);
  Array2D<int> CRm(p->numChains, p->loopSteps);
  vector<double> delta_tot(p->nCR, 1.0);
  vector<double> sd(p->nvar, 0.0);
  vector<double> delta_normX(p->numChains, 0.0);

  double delta_sum(0.0);
  double pCR_sum(0.0);

  for (int t = prevLines + 1; t < p->maxEvals; ++t) {   // loop over time/generations

    // beginning of loop, generate crossover probabilities
    if (genNumber == 0) { 
      gen_CR(rng, pCR, CRm, L); 
      numAccepted = 0;
    }

    for (int i = 0; i < p->numChains; ++i) {
      // loop over individual chains to generate new proposals
      if (p->deltaMax > 1) 
        rng->uniform_int(1, &delta, 1, p->deltaMax + 1);
      else 
        delta = 1;

      // generate proposal
      // PE: i = chain index, j = variable index
      for (int j = 0; j < p->nvar; ++j) 
        proposal(i,j) = state(t - 1, i, j);
      // pick pairs
      vector<int> r1(delta, 0);
      vector<int> r2(delta, 0);
      for (int a(0); a < delta; ++a) {
        rng->uniform_int(1, &r1[a], 0, p->numChains - 1);
        rng->uniform_int(1, &r2[a], 0, p->numChains - 1);
        if (r1[a] >= i) 
          ++r1[a];
        if (p->numChains > 2) {
          while (r2[a] == r1[a] || r2[a] == i) 
            ++r2[a] %= p->numChains;
        }
      }
      // PE: updateDim keeps track of how many variables will be replaced;
      // differs for each chain in each generation
      // (nFree - # variables which crossover events set equal to current state value)
      updateDim[i] = p->nfree;
      for (int j = 0; j < p->nvar; ++j) {
        if (! p->varLock[j]) {
          updatePar[j] = true;
          // calculate random values
          rng->uniform(1, &drand);
          e[j] = p->noise*(2.0*drand - 1.0);
          rng->gaussian(1, &epsilon[j], 0.0, bstar[j]);
          // compute pairwise comparisons
          pairDiff[j] = 0.0;
          for (int a(0); a < delta; ++a) {
            if (r1[a] != r2[a]) 
              pairDiff[j] += state(t - 1, r1[a], j) - state(t - 1, r2[a], j);
          }
          // check for crossover events
          crossRate = 1.0*CRm(i, genNumber) / p->nCR;
          if (crossRate < 1.0) {
            rng->uniform(1, &drand);
            if (drand < 1.0 - crossRate) {
              step[j] = 0.0;
              updatePar[j] = false;
              --updateDim[i];
            }
          }
        }
      }
      if (updateDim[i] > 0) {
        // every 5 generations, set gamma = 1.0 to promote long jumps outside local mode
        if ((t > 1) && ((t % 5) == 0))
          gamma = 1.0;
        else
          gamma = 2.38/sqrt(2.0*updateDim[i]*delta);
        for (int j = 0; j < p->nvar; ++j) {
          if (updatePar[j]) {
            // calculate step for this dimension
            step[j] = (1 + e[j])*gamma*pairDiff[j] + epsilon[j];
            // update proposal
            proposal(i,j) = state(t - 1,i,j) + step[j];
          } else {
            proposal(i,j) = state(t - 1,i,j);
          }
        }
      } else {   // PE: entire proposal is actually identical to current state
        for (int j = 0; j < p->nvar; ++j) proposal(i,j) = state(t - 1,i,j);
        }
    }  // end of loop(i) over chains


    for (int i = 0; i < p->numChains; ++i) {
      // loop over individual chains to calculate likelihoods of proposals
      // and determine acceptances
      if (updateDim[i] > 0) {
        do_calc = 1;
        for (int j = 0; j < p->nvar; ++j) {
          if (! p->varLock[j]) {
            if (proposal(i,j) < p->varLo[j] || proposal(i,j) > p->varHi[j]) {
              do_calc = 0;
              break;
            }
          }
        }
        if (p->recalcLik + inBurnIn > 0) {
          // [.]TAG: lik
          // lik2[i][t - 1] = p->fun(i, t - 1, state.pt(t - 1, i), p->extraData, true);
          lik(t - 1,i) = p->fun(i, t - 1, state.pt(t - 1, i), p->extraData, true);
          nLikelihoodEvals++;
        }
        if (do_calc) {
          // [.]TAG: lik
          // lik2[i][t] = p->fun(i, t, proposal(i), p->extraData, false);
          lik(t,i) = p->fun(i, t, proposal(i), p->extraData, false);
          nLikelihoodEvals++;
          // if (p->vflag) cout << ". Likelihood = " << lik(t,i) << endl;
        } else   // [.] TAG: lik
          // lik2[i][t] = -INFINITY;
          lik(t,i) = -INFINITY;
      } else {
        for (int j = 0; j < p->nvar; ++j) 
          proposal(i,j) = state(t - 1,i,j);
        if (p->recalcLik + inBurnIn > 0) {
          // [.]TAG: lik
          // lik2[i][t] = p->fun(i, t, proposal(i), p->extraData, true);
          lik(t,i) = p->fun(i, t, proposal(i), p->extraData, true);
          nLikelihoodEvals++;
        } else {  // [.]TAG: lik
          // lik2[i][t] = lik2(t - 1,i);
          lik(t,i) = lik(t - 1,i);
        }
      }
    }

    for (int i = 0; i < p->numChains; ++i) {
      // [.]TAG: lik
      // double newLikelihood = lik2[i][t];
      // double prevLikelihood = lik2[i][t - 1];
      double newLikelihood = lik(t,i);
      double prevLikelihood = lik(t - 1,i);
      if (newLikelihood == -INFINITY) 
        acceptStep[i] = 0;
      else if (newLikelihood >= prevLikelihood) 
        acceptStep[i] = 1;
      else {
        rng->uniform(1, &drand);
        if (log(drand) < newLikelihood - prevLikelihood) 
          acceptStep[i] = 1;
        else 
          acceptStep[i] = 0;
      }
      if (p->verboseLevel > 1) {
        cout << t << " " << prevLikelihood 
             << " <- " << acceptStep[i] << "|" << updateDim[i]
             << " -> " << newLikelihood << endl;
      }
      if (acceptStep[i]) {
        ++numAccepted;
        for (int j = 0; j < p->nvar; ++j) 
          state(t,i,j) = proposal(i,j);
      } else {
        for (int j = 0; j < p->nvar; ++j)
          state(t,i,j) = state(t - 1,i,j);
        // [.]TAG: lik
        // lik2[i][t] = lik2[i][t - 1];
        lik(t,i) = lik(t - 1,i);
      }
    }  // end of loop(i) over individual chains


    // ---------------------------------------------------------------------
    // update pCR if in burn-in phase
    if (inBurnIn && p->pCR_update) {
      // get standard deviations between the chains
      for (int j = 0; j < p->nvar; ++j) {
        sd[j] = gsl_stats_sd(state.pt(t,0,j), p->nvar, p->numChains);
//          if (! p->varLock[j] && sd[j] == 0.0) {
//            bstar[j] *= 2;
//            cerr << "Variable " << j << " has collapsed. Increasing stochasticity to " << bstar[j] << "." << endl;
//          }
      }
      for (int i = 0; i < p->numChains; ++i) {
        if (acceptStep[i]) {
          // get Euclidian distance
          delta_normX[i] = 0.0;
          for (int j = 0; j < p->nvar; ++j) {
            if (! p->varLock[j] && sd[j] > 0.0) {
              delta_normX[i] += gsl_pow_2((state(t,i,j) - state(t - 1,i,j))/sd[j]);
            }
          }
          delta_tot[CRm(i, genNumber) - 1] += delta_normX[i];
        }
      }
    }


    if (++genNumber >= p->loopSteps) {
      // every p->loopSteps [default = 10] generations, do the following
      //    1. If in burn-in phase, update CR delta
      //    2. Check for outlier chains
      //    3. If *not* in burn-in phase, check for Gelman-Rubin convergence
      //       *if* it's been p->gelmanEvals loopSteps since last G-R check
      //    4. Otherwise, check to see if we've passed the burn-in phase generation limit
      genNumber = 0;

      if (inBurnIn && p->pCR_update) {
        // get total delta
        delta_sum = 0.0;
        for (int m = 0; m < p->nCR; ++m) {
          delta_sum += delta_tot[m];
          totalSteps[m] += L[m];
        }
        if (delta_sum > 0.0) {
          pCR_sum = 0.0;
          for (int m = 0; m < p->nCR; ++m) {
            pCR[m] = (delta_tot[m]/delta_sum) / totalSteps[m];  // relative jump size per step
            pCR_sum += pCR[m];
          }
          for (int m = 0; m < p->nCR; ++m) {
            pCR[m] /= pCR_sum;
          }
        }
      }

      if (p->collapseOutliers && t < p->reenterBurnin*p->maxEvals) {
        // remove outlier chains
        vector<double> meanlik(p->numChains, -INFINITY);
        vector<bool> outliers(p->numChains, false);
        // [ ]TAG: lik
        check_outliers(t, lik, meanlik, outliers);
        int best_chain = gsl_stats_max_index(meanlik.data(), 1, p->numChains);
        for (int i = 0; i < p->numChains; ++i) {
          if (outliers[i] && i != best_chain) {
            // chain is an outlier
            // [.] TAG: lik
            // lik2[i][t] = lik2[best_chain][t]
            lik(t,i) = lik(t,best_chain);
            for (int j = 0; j < p->nvar; ++j) 
              state(t,i,j) = state(t,best_chain,j);
            // PE: if we're not currently in burn-in (and p->burnIn > 0), re-enter it!
            if (! inBurnIn && p->burnIn > 0) {
              if (p->verboseLevel > 0) {
                printf("[%d] Outlier chain detected [%d] outside of burn-in.", t, i);
                printf(" Moving to best chain [%d] and re-entering burn-in.\n", best_chain);
              }
              burnInStart = t;
              curRun = 0;
              inBurnIn = 1;
            }
          }
        }
      }

      // check if in burn-in period
      if (! inBurnIn) {
        // calculate Gelman-Rubin convergence
        if (++curRun >= p->gelmanEvals) {
          if (p->verboseLevel > 0)
            printf("[%d] performing convergence diagnostics:", t);
          gelman_rubin(state, scaleReduction, p->varLock, burnInStart + p->burnIn, t);
          // estimate variance
          int exitLoop(p->nvar);
          if (p->verboseLevel > 0)
            printf(" GR (it %d): ", t);
          for (int j(0); j < p->nvar; ++j) {
            if (! p->varLock[j]) {
              // scale reduction factor
              if (scaleReduction[j] < p->scaleReductionCrit)
                --exitLoop;
              if (p->verboseLevel > 0)
                printf("[%d]%.5f ", j, scaleReduction[j]);
            } else 
              --exitLoop;
          }
          if (p->verboseLevel > 0)
            printf("\n");
          // check for convergence
          if (exitLoop <= 0) {
            if (p->verboseLevel > 0)
              printf("Convergence detected!\n");
            converged = true;
            break;
          }
          // reset counter
          curRun = 0;
        }
      } else {
        // check if burn-in is finished
        if (t >= burnInStart + p->burnIn) {
          inBurnIn = 0;
          if (p->verboseLevel > 0)
            printf("[%d] exiting burn-in.\n", t);
        }
      }
    }  // end of "if (++genNumber >= p->loopSteps)"


    ++ireport;
    if (ireport >= p->report_interval) {
      ireport = 0;
      for (int i = 0; i < p->numChains; ++i) {
        for (int j = 0; j < p->nvar; j++)
          tempParams[j] = state(t,i,j);
        string paramString = theModel->PrintModelParamsHorizontalString(tempParams);
        *oout[i] << paramString << " ";
        // TAG: lik
        // *oout[i] << lik2[i][t] << " " << (t < burnInStart + p->burnIn) << " " << " ";
        *oout[i] << lik(t,i) << " " << (t < burnInStart + p->burnIn) << " " << " ";
        for (int j(0); j < p->nCR; ++j) 
          *oout[i] << pCR[j] << " ";
        *oout[i] << acceptStep[i] << endl;
      }
    }
  }  // end of loop(i) over generations
  
  

  for (int i(0); i < p->numChains; ++i) {
    // close output files
    if (oout[i] != NULL) 
      delete oout[i];
  }

  // PE: memory cleanup of added stuff
  free(tempParams);
  for (int i = 0; i < p->numChains; i++)
    delete[] lik2[i];
  
  
  if (p->verboseLevel > 0)
    printf("%d likelihood function calls\n", nLikelihoodEvals);
  if (! converged) {
    if (p->verboseLevel > 0)
      printf("Maximum number of iterations reached.\n");
    return DREAM_EXIT_MAX_ITERATIONS;
  }
  else
    return DREAM_EXIT_CONVERGENCE;
}

