#ifndef __RNGSTREAM_H__
#define __RNGSTREAM_H__

#include <cstdlib>
#include <iostream>
#include <ctime>

namespace rng {
  inline void make_discrete( size_t n, const double* w, double* x ) {
    x[0] = w[0];
    for (size_t i = 1; i < n; ++i) x[i] = x[i - 1] + w[i];
  }
  
  inline double dprob( size_t i, const double* w ) {
    return w[i] - ((i > 0) ? w[i - 1] : 0);
  }

  class RngStream {
    public:
      RngStream() : is_alloc(0) {}
      virtual ~RngStream() {}
      virtual void alloc( unsigned long seed = time(NULL) ) = 0;
      virtual void free() = 0;

      inline void leapfrog( size_t k, size_t n ) {}

      virtual void uniform( size_t n, double* r, double a = 0.0, double b = 1.0 ) = 0;
      virtual void uniform_int( size_t n, int* r, int a = 0, int b = 100 ) = 0;
      virtual void poisson( size_t n, int* k, double lambda ) = 0;

      int discrete( size_t n, const double* w ) {
        double* w0 = new double[n + 1];
        make_discrete(n, w, w0);
        int j = discrete_x(n, w0);
        delete[] w0;
        return j;
      }

      inline int pick( const double* x, size_t n ) {
        if (x[n - 1] <= 0.0) return -1;
        double r;
        uniform(1, &r, 0, x[n - 1]);
        int i = (int) (r/x[n - 1]*n);
        if (i < 0) {
          std::cerr << "RngStream::pick : "
                    << r << " " << x[n - 1] << " " << n << std::endl;
          return -1;
        }
        if (r > x[i]) {
          while (r > x[i] && i < (int)n) ++i;
        } else {
          while (i > 0) {
            if (r > x[i - 1]) break;
            --i;
          }
        }
        return i;
      }
      inline int discrete_x( size_t n, const double* w) { return pick(w,n ); }

      virtual void multinomial( size_t k, size_t n, const double* p, unsigned* a ) = 0;
      virtual void gaussian( size_t n, double* r, double mu, double sigma ) = 0;
      virtual void shuffle( int* x, size_t n ) = 0;
      virtual void shuffle( double* x, size_t n ) = 0;

//      template<typename T> void shuffle(T* x, size_t n) {
//        std::cerr << "Shuffling not implemented!" << std::endl;
//      }
//      virtual void multinomial(size_t k, size_t n, const double* p, unsigned int* x) {
//        std::cerr << "Multinomial not implemented!" << std::endl;
//      }

    protected:
      int is_alloc;
  };
}

#endif // __RNGSTREAM_H__

