#ifndef __GSLRNG_H__
#define __GSLRNG_H__

#include "Rng.h"
#include "GSLStream.h"

namespace rng {
  class GSLRng : public Rng {
    public:
      GSLRng() {}
      virtual ~GSLRng() { free(); }

      inline void alloc(size_t n) {
        if (streams.size() == 0) {
          streams.resize(n,NULL);
          for (size_t i = 0; i < n; ++i) {
            streams[i] = new GSLStream;
            streams[i]->alloc(seed+i);
          }
        }
      }

      inline void free() {
        while (streams.size() > 0) {
          streams.back()->free();
          streams.pop_back();
        }
      }

      inline void set_type(const gsl_rng_type* T) { type = T; }

    protected:
      const gsl_rng_type* type;
      gsl_rng** rng;
  };
}

#endif

