/* Generally useful definitions (sizes of different kind of buffers, etc.) */

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


/* TYPE OF INPUT ERROR/WEIGHT IMAGE */
#define  WEIGHTS_ARE_SIGMAS     100  // "weight image" pixel value = sigma
#define  WEIGHTS_ARE_VARIANCES  110  // "weight image" pixel value = variance (sigma^2)
#define  WEIGHTS_ARE_WEIGHTS    120  // "weight image" pixel value = weight

#define  MASK_ZERO_IS_GOOD        10  // "standard" input mask format (good pixels = 0)
#define  MASK_ZERO_IS_BAD         20  // alternate input mask format (good pixels = 1)


#endif /* _DEFINITIONS_H_ */
