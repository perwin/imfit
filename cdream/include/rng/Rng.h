#ifndef __RNG_H__
#define __RNG_H__

#include "RngStream.h"
#include <vector>

namespace rng {
  class Rng {
    public:
      Rng() : seed(time(NULL)) {}
      virtual ~Rng() {}
      virtual void alloc(size_t n) = 0;
      virtual void free() = 0;
      inline void set_seed(unsigned long s) { seed = s; }

      inline RngStream* operator[](size_t i) { return streams[i]; }

      inline void uniform(size_t k, size_t n, double* r, double a = 0.0, double b = 1.0) {
        if (streams.size() > k) {
          streams[k]->uniform(n,r,a,b);
        }
      }

      inline unsigned long getSeed() const { return seed; }

    protected:
      std::vector<RngStream*> streams;
      unsigned long seed;
  };
}

#endif

