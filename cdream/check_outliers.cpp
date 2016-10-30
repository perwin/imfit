#include "dream.h"

void check_outliers( int t, Array2D<double>& lik, vector<double>& meanlik,
                    vector<bool>& outliers )
{
  int t0 = t/2;
  int numChains = lik.n_y();
  double Q1;
  double Q3;
  double IQR;
  double UR;
  vector<double> liksrt(numChains);
  if ((int) meanlik.size() < numChains) 
    meanlik.resize(numChains, -INFINITY);
  for (int i(0); i < numChains; ++i) {
    meanlik[i] = gsl_stats_mean(lik.pt(t0,i), numChains, t - t0);
    liksrt[i] = meanlik[i];
  }
  sort(liksrt.begin(), liksrt.end());
  Q1 = gsl_stats_quantile_from_sorted_data(liksrt.data() ,1 ,liksrt.size() ,0.25);
  Q3 = gsl_stats_quantile_from_sorted_data(liksrt.data() ,1 ,liksrt.size() ,0.75);
  IQR = Q3 - Q1;
  UR = Q1 - 2*IQR;
  for (int i(0); i < numChains; ++i) {
    outliers[i] = false;
    if (meanlik[i] < UR) 
      outliers[i] = true;
  }
}


