/* FILE: count_cpu_cores.cpp --------------------------------------- */

// Code for determining the number of *physical* CPU cores on a computer;
// useful for determining how many threads to use for OpenMP and (especially)
// FFTW.
// This has a single function in two versions, one for Linux (parses the
// /proc/cpuinfo file) and one for macOS (uses sysctlbyname function).

// Copyright 2022 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <stdio.h>

#include "utilities_pub.h"
#include "count_cpu_cores.h"

using namespace std;


#ifdef LINUX
// Linux version of GetPhysicalCoreCount() -- parses /proc/cpuinfo "file"

#include <stdlib.h>

int GetPhysicalCoreCount( )
{
  string  fileName = "/proc/cpuinfo";
  ifstream  inputFileStream;
  string  inputLine;
  vector<string>  inputLines;
  vector<string>  tokens;
  set<string>  physIDStrings;   // use std::set so we only store unique strings
  // search strings
  string  physIDStr = "physical id";
  string  coreStr = "cpu cores";
  int  nSockets = 0;
  int  nCoresPerSocket = 0;

  inputFileStream.open(fileName.c_str());
  if( ! inputFileStream ) {
     cerr << "Error opening input stream for file " << fileName.c_str() << endl;
  }
  while ( getline(inputFileStream, inputLine) ) {
    if ((inputLine.size() > 0) && (inputLine[0] != '#')) {
      inputLines.push_back(inputLine);
    }
  }
  inputFileStream.close();

  for (int i = 0; i < inputLines.size(); i++) {
    // collect unique instances of "physical id" string
    if (inputLines[i].find(physIDStr) != string::npos) {
      physIDStrings.insert(inputLines[i]);
    }
    else if (inputLines[i].find(coreStr) != string::npos) {
      // get nCoresPerSocket (assume all instances are the same and just use the last)
      SplitString(inputLines[i], tokens);
      nCoresPerSocket = atol(tokens[tokens.size() - 1].c_str());
    }
  }
  
  nSockets = physIDStrings.size();
  return nSockets * nCoresPerSocket;
}

#else
// macOS version -- much simpler!

#include <sys/types.h>
#include <sys/sysctl.h>

int GetPhysicalCoreCount( )
{
  int  nCores;
  size_t  size = sizeof(nCores);
  
  // we use hw.perflevel0.physicalcpu instead of hw.physicalcpu to ensure we only count
  // high-performance cores on Apple Silicon CPUs
  sysctlbyname("hw.perflevel0.physicalcpu", &nCores, &size, NULL, 0);
  return nCores;
}

#endif



/* END OF FILE: count_cpu_cores.cpp ------------------------------------ */
