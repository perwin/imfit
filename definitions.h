/*! \file
    \brief Generally useful definitions (debugging levels, fit statistics, solvers, etc.) 

    Definitions of constants referring to debugging levels, which fit statisitic
    is being used, which minimizer/solver is being used, format of error/weight
    image, definition of good/bad pixels in mask, and max buffer sizes for
    lines of text and for filenames.
 */

#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_


#define  MAXLINE   1024
#define MAX_FILENAME_LENGTH  512



/* DEBUGGING LEVELS: */
const int  DEBUG_NONE  =             0;
const int  DEBUG_BASIC =             1;
const int  DEBUG_2     =             2;
const int  DEBUG_3     =             3;
const int  DEBUG_ALL   =            10;


/* OPTIONS FOR FIT STATISTICS: */
#define FITSTAT_CHISQUARE       1   //! standard chi^2
#define FITSTAT_CASH            2   //! standard (minimal) Cash statistic
#define FITSTAT_POISSON_MLR     3   //! Poisson Maximum Likelihood Ratio statistic


/* SOLVER OPTIONS: */
#define NO_FITTING             0
#define MPFIT_SOLVER           1
#define DIFF_EVOLN_SOLVER      2
#define NMSIMPLEX_SOLVER       3
#define ALT_SOLVER             4
#define GENERIC_NLOPT_SOLVER   5

/* TYPE OF INPUT ERROR/WEIGHT IMAGE */
#define  WEIGHTS_ARE_SIGMAS     100  //! "weight image" pixel value = sigma
#define  WEIGHTS_ARE_VARIANCES  110  //! "weight image" pixel value = variance (sigma^2)
#define  WEIGHTS_ARE_WEIGHTS    120  //! "weight image" pixel value = weight

#define  MASK_ZERO_IS_GOOD        10  //! "standard" input mask format (good pixels = 0)
#define  MASK_ZERO_IS_BAD         20  //! alternate input mask format (good pixels = 1)


#endif /* _DEFINITIONS_H_ */
