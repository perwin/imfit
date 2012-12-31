/* FILE: integrator.cpp ------------------------------------------------ */
/* VERSION 0.1
 *
 * Code for performing a line-of-sight integration using GSL QAGS integration.
 * Everything inside a single function, as local variables, to ensure thread
 * safety (e.g., for use with OpenMP).
 *
 * NOTE: Trial use of gsl_integration_qagi (integrating from -infty to +infty
 * sometimes worked, but sometimes failed (e.g., for 3D exponential disk when
 * inclination >~ 55 deg, though i = 90 worked).
 * Conclusion: safer to stick with gsl_integration_qags and do finite integration
 * to +/- large_number (s1 and s2).
 *
 * NOTE: Trial change of LIMIT_SIZE from 1000 to 10000, or RELATIVE_TOL from 1.0e-6 to
 * 1.0e-5, had no effect on integration failures for edges of moderately inclined 
 * Ferrers bar, so not much reason to change them.
 */
 
//#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include "integrator.h"

#define LIMIT_SIZE   1000
#define RELATIVE_TOL  1.0e-6

double  Integrate( gsl_function F, double s1, double s2 )
{
  double  result, error;
  int  status;
  gsl_integration_workspace * workspace;
  
  // allocate and free the workspace object here (referencing it with a local
  // variable) to ensure thread safety
  workspace = gsl_integration_workspace_alloc(LIMIT_SIZE);
//  gsl_integration_qagi(&F, 0, 1e-6, 1000, workspace, &result, &error);
  status = gsl_integration_qags(&F, s1, s2, 0, RELATIVE_TOL, LIMIT_SIZE, workspace, &result, &error);
  gsl_integration_workspace_free(workspace);
  
  return result;
}

  //  gsl_integration_qags(&F, -integLimit, integLimit, 0, 1e-7, 1000, 
//  						workspace, &totalIntensity, &error); 

//      double f (double x, void * params) {
//        double alpha = *(double *) params;
//        double f = log(alpha*x) / sqrt(x);
//        return f;
//      }
//      
//      int
//      main (void)
//      {
//        gsl_integration_workspace * w 
//          = gsl_integration_workspace_alloc (1000);
//        
//        double result, error;
//        double expected = -4.0;
//        double alpha = 1.0;
//      
//        gsl_function F;
//        F.function = &f;
//        F.params = &alpha;
//      
//        gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
//                              w, &result, &error); 
//      
//        printf ("result          = % .18f\n", result);
//        printf ("exact result    = % .18f\n", expected);
//        printf ("estimated error = % .18f\n", error);
//        printf ("actual error    = % .18f\n", result - expected);
//        printf ("intervals =  %d\n", w->size);
//      
//        gsl_integration_workspace_free (w);
//      
//        return 0;
//      }
