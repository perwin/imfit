#include "dream.h"

void gelman_rubin( Array3D<double>& state,  vector<double>& scaleReduction,
				 	const int* lockVar, int first, int numIter, int adjustDF ) 
{
  int maxEvals = state.n_x();
  int numChains = state.n_y();
  int numPars = state.n_z();
  int t = (numIter <= 0) ? maxEvals : numIter;

  Array2D<double> chainMean(numChains, numPars);
  Array2D<double> chainVar(numChains, numPars);
  Array2D<double> chainMean2(numChains, numPars);
  vector<double> chainBetMean(numPars, 0.0);
  vector<double> chainBetMean2(numPars, 0.0);
  vector<double> chainBetVar(numPars, 0.0);
  vector<double> meanVar(numPars, 0.0);
  vector<double> varVar(numPars, 0.0);
  vector<double> estimatedVar(numPars, 0.0);

  chainMean.set_all(0.0);
  chainVar.set_all(0.0);
  if ((int) scaleReduction.size() < numPars) 
    scaleReduction.resize(numPars, 0.0);

  int i, j;
  int t0 = (t + first)/2;
  int n = t - t0 - 1;

  // get within-chain means and variances
  for (i = 0; i < numChains; ++i) {
    for (j = 0; j < numPars; ++j) {
      if (! lockVar[j]) {
        chainMean(i,j) = gsl_stats_mean(state.pt(t0,i,j), numChains*numPars, n);
        chainMean2(i,j) = gsl_pow_2(chainMean(i,j));
        chainVar(i,j) = gsl_stats_variance_m(state.pt(t0,i,j), numChains*numPars, n, chainMean(i,j));
      }
    } 
  }

  // get between-chain mean and variance
  for (j = 0; j < numPars; ++j) {
    if (! lockVar[j]) {
      chainBetMean[j] = gsl_stats_mean(chainMean.pt(0,j), numPars, numChains);
      chainBetMean2[j] = gsl_stats_mean(chainMean2.pt(0,j), numPars, numChains);
      chainBetVar[j] = gsl_stats_variance_m(chainMean.pt(0,j), numPars, numChains, chainBetMean[j]);
      meanVar[j] = gsl_stats_mean(chainVar.pt(0,j), numPars, numChains);
      varVar[j] = gsl_stats_variance_m(chainVar.pt(0,j), numPars, numChains, meanVar[j]);
    }
  }

  // estimate variance
  for (j = 0; j < numPars; ++j) {
    if (! lockVar[j]) {
      int m = numChains;
      double Bdivn = chainBetVar[j];           // between-chain variance
      double W = meanVar[j];                   // mean sample variance
      double sigma2 = (n - 1.0)/n*W + Bdivn;   // estimated variance
      double Vhat = sigma2 + Bdivn/numChains;  // pooled variance estimate
      double Rhat = Vhat/W;                    // potential scale reduction factor
      double cv1 = gsl_stats_covariance_m(chainVar.pt(0,j), numPars, chainMean2.pt(0,j), numPars,
                                          numChains, meanVar[j], chainBetMean2[j]);
      double cv2 = gsl_stats_covariance_m(chainVar.pt(0,j), numPars, chainMean.pt(0,j), numPars,
                                          numChains, meanVar[j], chainBetMean[j]);
      // var.V <- (+ (1 + 1/Nchain)^2 * (2 * b^2)/(Nchain - 1) + 
      //  2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
      double varVhat = gsl_pow_2((n - 1.0)/n)*varVar[j]/m
                       + gsl_pow_2((m + 1.0)/(m*n))*(2./(m - 1.0))*gsl_pow_2(Bdivn*n)
                       + 2*(m + 1.0)*(n - 1.0)/(m*n*n)*(n/m)*(cv1 - 2.*chainBetMean[j]*cv2);
      double df = 2.0*gsl_pow_2(Vhat)/varVhat;

      if (adjustDF)
        scaleReduction[j] = (W > 0.0) ? sqrt(Rhat*(df + 3.)/(df + 1.)) : 1.0;
      else
        scaleReduction[j] = (W > 0.0) ? sqrt(Rhat) : 1.0;

      /*
      printf("[%d] Between chain variance = %g\n",j,Bdivn);
      printf("[%d] Mean sample variance   = %g\n",j,W);
      printf("[%d] Estimated variance     = %g\n",j,sigma2);
      printf("[%d] Vhat                   = %g\n",j,Vhat);
      printf("[%d] var(Vhat)              = %g\n",j,varVhat);
      printf("[%d] R estimate             = %g\n",j,Rhat);
      printf("[%d] CV1                    = %g\n",j,cv1);
      printf("[%d] CV2                    = %g\n",j,cv2);
      printf("[%d] df                     = %g\n",j,df);
      */
    }
  }
}

