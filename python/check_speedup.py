#!/usr/bin/env python
#
# A Python script to compare the elapsed time of two different imfit fits.
# The assumed use is to copare a standard, OpenMP-enabled fit (the first file)
# against the same fit with OpenMP turned off (--max-threads=1 option).
# For simplicity, we assume the OpenMP-enabled fit should be at least two times
# faster; if it isn't an error message is printed and 1 is returned (signaling
# failure).
#
# We assume each output file has, as its third-to-last line, something like
# "(Elapsed time: 0.291234 sec for fit, 0.294756 sec total)"
#
# The script uses sys.exit() to return either 0 for success or 1 for some kind
# of failure; this is for use with shell scripts for regression tests, etc.

from __future__ import print_function

import sys, os, optparse, re


# look for and capture elapsed time in output line from imfit
findTime = re.compile(r"Elapsed time:\s(?P<time>[0-9]+\.[0-9]+)\ssec")

# required ratio
MINIMUM_TIME_RATIO = 2.0

# predefine some ANSI color codes
RED  = '\033[31m' # red
NC = '\033[0m' # No Color


def CompareTimes( line1, line2, fileName1, fileName2 ):
    """
    Simple comparison of two different imfit output files, with the assumption
    that the elapsed time for second (line2 from fileName2) should be at least
    twice that for the first.
    """
    s1 = findTime.search(line1)
    s2 = findTime.search(line2)
    if (s1 is None) or (s2 is None):
        txt = "\n\t" + RED + ">>> WARNING:" + NC
        print(txt + " unable to find \"elapsed time\" output in %s and/or %s!\n" % (fileName1, fileName2))
        sys.exit(1)
    t1 = float(s1.group('time'))
    t2 = float(s2.group('time'))
    ratio = t2/t1
    #print(t1, t2, ratio)
    if (ratio < MINIMUM_TIME_RATIO):
        return False, ratio
    else:
        return True, ratio



def main(argv=None):

    usageString = "%prog output_file_1 output_file_2\n"
    parser = optparse.OptionParser(usage=usageString, version="%prog ")

    (options, args) = parser.parse_args(argv)
    
    if (len(args)) < 3:
        msg = "you must supply two text-file names!\n"
        print(RED + "ERROR: " + msg + NC)
        sys.exit(1)
    regularFile = args[1]
    noThreadingFile = args[2]
    if not os.path.exists(regularFile):
        msg = "unable to find file %s!\n" % regularFile
        print(RED + "ERROR: " + msg + NC)
        sys.exit(1)
    if not os.path.exists(noThreadingFile):
        msg = "unable to find file %s!\n" % noThreadingFile
        print(RED + "ERROR: " + msg + NC)
        sys.exit(1)
    
    # get output lines with elapsed times
    elapsed_time1 = open(regularFile).readlines()[-3]
    elapsed_time2 = open(noThreadingFile).readlines()[-3]
    result, ratio = CompareTimes(elapsed_time1, elapsed_time2, regularFile, noThreadingFile)
    if (result is False):
        txt = "\n\t" + RED + ">>> WARNING:" + NC
        msg = txt + " OpenMP-enabled fit was slower than expected!"
        msg += " (ratio of max-thread=1 fit to standard fit = %g)!\n" % ratio
        print(msg)
        sys.exit(1)
    else:
        print(" OK.")
        sys.exit(0)


if __name__ == '__main__':
    
    main(sys.argv)
