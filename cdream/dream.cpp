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

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "dream.h"
#include "dream_params.h"

// MPI variables
int mpi_rank = 0;
int mpi_ntasks = 0;
clock_t calc_time = 0;
clock_t run_time = 0;
clock_t last_time = 0;


int dream( const dream_pars* p, rng::RngStream* rng ) {
  int inBurnIn = (p->burnIn > 0);
  int burnInStart = 0;
  int genNumber = 0;

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
  
  if (mpi_rank == 0) {
    // MCMC chains
    Array3D<double> state(p->maxEvals, p->numChains,p->nvar);
    Array2D<double> lik(p->maxEvals + 1, p->numChains);

    Array2D<double> proposal(p->numChains, p->nvar);
    Array2D<double> proposal_two(p->numChains, p->nvar);
    proposal.set_all(0.0);
    proposal_two.set_all(0.0);

    vector<double> scaleReduction(p->nvar, 0.0);
    vector<double> pCR(p->nCR, 1.0/p->nCR);

    // =========================================================================
    // read previous state

    int prevLines = dream_restore_state(p, state, lik, pCR, inBurnIn);

    // =========================================================================
    // open output file

    ostringstream chainFilename;
    vector<ostream*> oout;
    ios_base::openmode fmode = (p->appendFile) ? (ios_base::out | ios_base::app) : ios_base::out;

    oout.resize(p->numChains, &cout);
    if (p->outputRootname != "" && p->outputRootname != "-") {
      for (int i = 0; i < p->numChains; ++i) {
        chainFilename.str("");
        chainFilename << p->outputRootname << "." << i << ".txt";
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
    // Initialize with latin hypbercube sampling if not a resumed run

    int do_calc(1);

    if (! p->appendFile) {
      Array2DView<double> initVar(state.n_y(), state.n_z(), state.pt(0,0,0));
      ArrayView<double> initLik(lik.n_y(), lik.pt(0,0));
      dream_initialize(p, rng, initVar, initLik);

      // save initial state of each chain
      for (int i = 0; i < p->numChains; ++i) {
        for (int j = 0; j < p->nvar; ++j) 
          *oout[i] << state(0,i,j) << " ";
        *oout[i]  << lik(0,i) << " " << inBurnIn << " " << 0 << " ";
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
//    bstar[4] = 1;   // WTF?
 
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

    for (int t = prevLines + 1; t < p->maxEvals; ++t) {

      // beginning of loop, generate crossover probabilities
      if (genNumber == 0) { 
        gen_CR(rng, pCR, CRm, L); 
        numAccepted = 0;
      }

      for (int i = 0; i < p->numChains; ++i) {
        if (p->deltaMax > 1) 
          rng->uniform_int(1, &delta, 1, p->deltaMax + 1);
        else 
          delta = 1;

        // generate proposal
        // PE: i = chain index, j = variable index
        for (int j = 0; j < p->nvar; ++j) 
          proposal(i,j) = state(t - 1,i,j);
        if (gammaGeneration++ == 5) {
          // PE: every 5 steps, we make a large proposal jump (gamma = 1.0)
          // PE: (here, we seem to be doing a simplified version of Vrugt et al.,
          // which sets gamma = 1 but otherwise does the full proposal generation)
          for (int j = 0; j < p->nvar; ++j) {
            if (p->varLock[j]) 
              continue;
            gammaGeneration = 0;
            do {
              rng->uniform_int(1, &r1, 0, p->numChains - 1);
              rng->uniform_int(1, &r2, 0, p->numChains - 1);
              if (r1 >= i) 
                ++r1;
              if (r2 >= i) 
                ++r2;
            } while (r1 == r2 && p->numChains > 2);
            step[j] = state(t - 1, r1, j) - state(t - 1, r2, j);
            proposal(i, j) = state(t - 1, i, j) + step[j];
          }
          // (TODO) if step == 0 ==> use Cholskey decomposition!!
        } else {
          // PE: Standard proposal generation
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
          // PE: updateDim keeps track of how many variables will be replaced
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
        }
      }

      // calculate likelihood and acceptance probability
      for (int i = 0; i < p->numChains; ++i) {
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
            lik(t - 1,i) = p->fun(i, t - 1, state.pt(t - 1, i), p->extraData, true);
          }
          if (do_calc) {
            lik(t,i) = p->fun(i, t, proposal(i), p->extraData, false);
            // if (p->vflag) cout << ". Likelihood = " << lik(t,i) << endl;
          } else 
            lik(t,i) = -INFINITY;
        } else {
          for (int j = 0; j < p->nvar; ++j) 
            proposal(i,j) = state(t - 1,i,j);
          if (p->recalcLik + inBurnIn > 0) {
            lik(t,i) = p->fun(i, t, proposal(i), p->extraData, true);
          } else {
            lik(t,i) = lik(t - 1,i);
          }
        }
      }

      for (int i = 0; i < p->numChains; ++i) {
        if (lik(t,i) == -INFINITY) 
          acceptStep[i] = 0;
        else if (lik(t,i) >= lik(t - 1,i)) 
          acceptStep[i] = 1;
        else {
          rng->uniform(1, &drand);
          if (log(drand) < lik(t,i) - lik(t - 1,i)) 
            acceptStep[i] = 1;
          else 
            acceptStep[i] = 0;
        }
        if (p->verboseLevel > 1) {
          cout << t << " " << lik(t - 1,i) 
               << " <- " << acceptStep[i] << "|" << updateDim[i]
               << " -> " << lik(t,i) << endl;
        }
        if (acceptStep[i]) {
          ++numAccepted;
          for (int j = 0; j < p->nvar; ++j) 
            state(t,i,j) = proposal(i,j);
        } else {
          for (int j = 0; j < p->nvar; ++j)
            state(t,i,j) = state(t - 1,i,j);
          lik(t,i) = lik(t - 1,i);
        }
      }

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
          check_outliers(t, lik, meanlik, outliers);
          int best_chain = gsl_stats_max_index(meanlik.data(), 1, p->numChains);
          for (int i = 0; i < p->numChains; ++i) {
            if (outliers[i] && i != best_chain) {
              // chain is an outlier
              lik(t,i) = lik(t,best_chain);
              for (int j = 0; j < p->nvar; ++j) 
                state(t,i,j) = state(t,best_chain,j);
              if (! inBurnIn && p->burnIn > 0) {
                cerr << "[" << t << "] Outlier chain detected [" << i << "] outside of burn in."
                     << " Moving to chain " << best_chain << " and re-entering burn in." << endl;
                burnInStart = t;
                curRun = 0;
                inBurnIn = 1;
              }
            }
          }
        }

        // check if in burn in period
        if (! inBurnIn) {
          // calculate Gelman-Rubin convergence
          if (curRun++ >= p->gelmanEvals) {
            cout << "[" << t << "] performing convergence diagnostics:";
            gelman_rubin(state, scaleReduction, p->varLock, burnInStart + p->burnIn, t);
            // estimate variance
            int exitLoop(p->nvar);
            cout << " GR (it " << t << "): ";
            for (int j(0); j < p->nvar; ++j) {
              if (! p->varLock[j]) {
                // scale reduction factor
                if (scaleReduction[j] < p->scaleReductionCrit)
                  --exitLoop;
                cout << "[" << j << "]" << scaleReduction[j] << " ";
              } else 
                --exitLoop;
            }
            cout << endl;
            // check for convergence
            if (exitLoop <= 0) 
              break;
            // reset counter
            curRun = 0;
          }
        } else {
          // check if burn-in is finished
          if (t >= burnInStart + p->burnIn) {
            inBurnIn = 0;
            cerr << "[" << t << "] exiting burn in." << endl;
          }
        }
      }

      ++ireport;
      if (ireport >= p->report_interval) {
        ireport = 0;
        for (int i = 0; i < p->numChains; ++i) {
          for (int j = 0; j < p->nvar; ++j) 
            *oout[i] << state(t,i,j) << " ";
          *oout[i] << lik(t,i) << " " << (t < burnInStart + p->burnIn) << " " << genNumber << " ";
          for (int j(0); j < p->nCR; ++j) 
            *oout[i] << pCR[j] << " ";
          *oout[i] << acceptStep[i] << endl;
        }
      }
    }

    for (int i(0); i < p->numChains; ++i) {
      if (oout[i] != NULL) 
        delete oout[i];
    }

  } else {
    // MPI slaves
  }

#ifdef USE_MPI
  last_time = clock();
  fprintf(stdout, "Calcluation efficiency (Slave %d): %f\n",
          mpi_rank, 1.0*calc_time/last_time);
  MPI::Finalize();
#endif

  return EXIT_SUCCESS;
}

