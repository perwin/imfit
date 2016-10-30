#ifndef __MKLRNG_H__
#define __MKLRNG_H__

#include "Rng.h"
#include "MKLStream.h"

#include <iostream>

namespace rng {
  class MKLRng : public Rng {
    public:
      MKLRng() { /* std::cout << "Using MKL." << std::endl; */ }
      virtual ~MKLRng() { free(); }

      inline void alloc(size_t n) {
        if (streams.size() == 0) {
          streams.resize(n,NULL);
          for (size_t i = 0; i < n; ++i) {
            streams[i] = new MKLStream;
            streams[i]->alloc(seed);
            streams[i]->leapfrog(i,n);
          }
        }
      }

      inline void free() {
        while (streams.size() > 0) {
          streams.back()->free();
          streams.pop_back();
        }
      }
  };
}

#endif

