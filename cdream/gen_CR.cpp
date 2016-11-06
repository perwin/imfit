#include <stdio.h>

#include "dream.h"

void gen_CR( rng::RngStream* rng, const vector<double>& pCR, Array2D<int>& CRm, 
			vector<unsigned>& L )
{
  size_t numChains(CRm.n_x());
  size_t loopSteps(CRm.n_y());
  size_t size(numChains*loopSteps);
  if (L.size() < pCR.size()) {
    fprintf(stderr, "genCR: Bad array size.\n");
    return;
  }
  // pick candidates for each crossover value
  rng->multinomial(pCR.size(), size, pCR.data(), L.data());
  size_t j = 0;
  for (size_t m = 0; m < L.size(); ++m) {
    for (size_t i = 0; i < L[m]; ++i) {
      CRm[j] = m + 1;
      ++j;
    }
  }
  rng->shuffle(CRm.pt(0,0), size);
}


