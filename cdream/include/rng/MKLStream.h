#ifndef __MKLSTREAM_H__
#define __MKLSTREAM_H__

#include <mkl_vsl.h>
#include "RngStream.h"

namespace rng {
  class MKLStream : public RngStream {
    public:
      MKLStream() : brng(VSL_BRNG_MCG31) {}
      virtual ~MKLStream() { free(); }

      inline void alloc(unsigned long seed = time(NULL)) {
        if (! is_alloc) {
          vslNewStream(&stream,brng,seed);
          is_alloc = 1;
        }
      }

      inline void free() {
        if (is_alloc) {
          vslDeleteStream(&stream); 
          is_alloc = 0;
        }
      }

      inline void copy(const MKLStream& s) {
        if (is_alloc) free();
        vslCopyStream(&stream,s.stream);
        is_alloc = 1;
      }

      inline void leapfrog(size_t k, size_t n) {
        vslLeapfrogStream(stream,k,n); 
      }

      inline void uniform(size_t n, double* r, double a = 0.0, double b = 1.0) {
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,n,r,a,b);
      }

      inline void uniform_int(size_t n, int* r, int a = 0, int b = 10) {
        viRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,n,r,a,b);
      }

      inline void poisson(size_t n, int* k, double lambda) {
        viRngPoisson(VSL_RNG_METHOD_POISSON_PTPE,stream,n,k,lambda);
      }

      inline void set_type(const MKL_INT type) { brng = type; }

      inline void shuffle(int* x, size_t n) {
        std::cerr << "Int shuffling not implemented!" << std::endl;
      }

      inline void shuffle(double* x, size_t n) {
        std::cerr << "Double shuffling not implemented!" << std::endl;
      }

      inline void gaussian(size_t n, double* r, double mu = 0.0, double sigma = 1.0) { 
        std::cerr << "Gaussian not implemented!" << std::endl;
      };

      inline void multinomial(size_t k, size_t n, const double* p, unsigned* a) {
        std::cerr << "Multinomial not implemented!" << std::endl;
      }


    protected:
      MKL_INT brng;
      VSLStreamStatePtr stream;
  };
}

#endif

