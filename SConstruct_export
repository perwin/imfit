# Scons "makefile" for imfit and related programs
# 
# To use this, type the following on the command line:
#    $ scons <target-name>
# where <target-name> is, e.g., "imft", "imfit_db", etc.
# "scons" by itself will build *all* targets
# 
# To clean up files associated with a given target:
#    $ scons -c <target-name>
# To clean up all targets (including programs):
#    $ scons -c


# *** SPECIAL STUFF ***
# To build a version with extra, experimental functions
#    $ scons --extra-funcs <target-name>
#
# To build a version *without* OpenMP enabled
#    $ scons --no-openmp <target-name>
#
# To build a version *without* FFTW threading:
#    $ scons --no-threading <target-name>
#
# To build a version with full debugging printouts:
#    $ scons define=DEBUG <target-name>
#
# To build export version ("fat" binaries, all libraries statically linked):
#    $ scons --fat --static <target-name>
#

# *** EXPORT CONFIGURATIONS ***
# MacOS X fat binaries
# $ scons --static --fat

# To add one or more directories to the header or library search paths:
#    $ scons --header-path=/path/to/header/dir
# OR $ scons --header-path=/path/to/header/dir:/alt/path:/another/path
#    $ scons --lib-path=/path/to/lib/dir
# etc.


# Copyright 2010, 2011, 2012, 2013 by Peter Erwin.
# 
# This file is part of Imfit.
# 
# Imfit is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.
# 
# Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
# 
# You should have received a copy of the GNU General Public License along
# with Imfit.  If not, see <http://www.gnu.org/licenses/>.


# Operating-system determination via os.uname:
# First element of tuple is basic OS; 3rd element is version number;
# 5th element is processor architecture (e.g., "i386", "sun4u", "i686")
#    os.uname()[0] = "Darwin" --> Mac OS X
#    os.uname()[0] = "SunOS" --> Solaris
#    os.uname()[0] = "Linux" --> Linux

import os

# Version definition for imfit + makeimage
PACKAGE_VERSION = "PACKAGE_VERSION=\"1.0b2\""

# the following is for when we want to force static linking to the GSL library
# (Change these if the locations are different on your system)
STATIC_GSL_LIBRARY_FILE_MACOSX = File("/usr/local/lib/libgsl.a")
STATIC_GSL_LIBRARY_FILE1_LINUX = File("/usr/lib/libgsl.a")
STATIC_GSL_LIBRARY_FILE2_LINUX = File("/usr/lib/libgslcblas.a")

# the following is for when we want to force static linking to the NLopt library
# (Change these if the locations are different on your system)
STATIC_NLOPT_LIBRARY_FILE_MACOSX = File("/usr/local/lib/libnlopt.a")
STATIC_NLOPT_LIBRARY_FILE1_LINUX = File("/usr/local/lib/libnlopt.a")


FUNCTION_SUBDIR = "function_objects/"

os_type = os.uname()[0]

cflags_opt = ["-O2", "-g0"]
cflags_db = ["-Wall", "-g3"]

base_defines = ["ANSI", "USING_SCONS", PACKAGE_VERSION]

# libraries needed for imfit, makeimage, psfconvolve, & other 2D programs
lib_list = ["fftw3", "cfitsio", "m"]
# libraries needed for profilefit and psfconvolve1d compilation
lib_list_1d = ["fftw3", "m"]


include_path = ["/usr/local/include", FUNCTION_SUBDIR]
lib_path = ["/usr/local/lib"]
link_flags = []

# system-specific setup
if (os_type == "Darwin"):   # OK, we're compiling on Mac OS X
	# Note: if for some reason you need to compile to 32-bit -- e.g., because
	# your machine is 32-bit only, or because the fftw3 and cfitsio libraries
	# are 32-bit, use the following
#	cflags_opt.append("-m32")
#	link_flags = ["-m32"]
	cflags_db = ["-Wall", "-Wshadow", "-Wredundant-decls", "-Wpointer-arith", "-g3"]
if (os_type == "Linux"):
	# change the following path definitions as needed
	include_path.append("/usr/include")
#	lib_list.append("pthread")
	if os.getlogin() == "erwin":
		include_path.append("/home/erwin/include")
		lib_path.append("/home/erwin/lib")
	# When compiled under Linux, -O3 causes mysterious "invalid pointer" error at end of run
	cflags_opt = ["-O3", "-g0"]
	cflags_db = ["-Wall", "-Wshadow", "-Wredundant-decls", "-Wpointer-arith", "-g3"]
	# silly Linux doesn't have OpenBSD string routines built in, so we'll have to include them
	base_defines = base_defines + ["LINUX"]
#	lib_list.append("gslcblas")
defines_opt = base_defines
#defines_db = base_defines + ["DEBUG"]
defines_db = base_defines

extra_defines = []


# Default settings for compilation
useGSL = True
useNLopt = True
useFFTWThreading = True
useOpenMP = True
useExtraFuncs = False
useStaticLibs = False
buildFatBinary = False


# Define some user options
AddOption("--lib-path", dest="libraryPath", type="string", action="store", default=None,
	help="colon-separated list of additional paths to search for libraries")
AddOption("--header-path", dest="headerPath", type="string", action="store", default=None,
	help="colon-separated list of additional paths to search for header files")
AddOption("--no-threading", dest="fftwThreading", action="store_false", 
	default=True, help="compile programs *without* FFTW threading")
AddOption("--no-gsl", dest="useGSL", action="store_false", 
	default=True, help="do *not* use GNU Scientific Library")
AddOption("--no-nlopt", dest="useNLopt", action="store_false", 
	default=True, help="do *not* use NLopt library")
AddOption("--no-openmp", dest="useOpenMP", action="store_false", 
	default=False, help="compile *without* OpenMP support")
AddOption("--extra-funcs", dest="useExtraFuncs", action="store_true", 
	default=False, help="compile additional FunctionObject classes for testing")

# Define some more arcane options (e.g., for making binaries for distribution)
AddOption("--static", dest="useStaticLibs", action="store_true", 
	default=False, help="force static library linking")
AddOption("--fat", dest="makeFatBinaries", action="store_true", 
	default=False, help="generate a \"fat\" (32-bit + 64-bit Intel) binary for Mac OS X")


if GetOption("headerPath") is not None:
	extraPaths = GetOption("headerPath").split(":")
	print "extra header search paths: ", extraPaths
	include_path += extraPaths
if GetOption("libraryPath") is not None:
	extraPaths = GetOption("libraryPath").split(":")
	print "extra library search paths: ", extraPaths
	lib_path += extraPaths
if GetOption("fftwThreading") is False:
	useFFTWThreading = False
if GetOption("useGSL") is False:
	useGSL = False
if GetOption("useNLopt") is False:
	useNLopt = False
if GetOption("useOpenMP") is True:
	useOpenMP = True
if GetOption("useExtraFuncs") is True:
	useExtraFuncs = True
if GetOption("useStaticLibs") is True:
	useStaticLibs = True
if GetOption("makeFatBinaries") is True:
	buildFatBinary = True


if useFFTWThreading:   # default is to do this
	lib_list.insert(0, "fftw3_threads")
	lib_list_1d.insert(0, "fftw3_threads")
	if (os_type == "Linux"):
		lib_list.append("pthread")
		lib_list_1d.append("pthread")
	extra_defines.append("FFTW_THREADING")

if useGSL:   # default is to do this
	if useStaticLibs:
		if (os_type == "Darwin"):
			lib_list.append(STATIC_GSL_LIBRARY_FILE_MACOSX)
		else:
			# assuming we're on a Linux system
			lib_list.append(STATIC_GSL_LIBRARY_FILE1_LINUX)
			lib_list.append(STATIC_GSL_LIBRARY_FILE2_LINUX)
	else:
		lib_list.append("gsl")
	
	# and stuff for 1D programs:
	lib_list_1d.append("gsl")
	if (os_type == "Linux"):
		lib_list.append("gslcblas")
		lib_list_1d.append("gslcblas")
else:
	extra_defines.append("NO_GSL")

if useNLopt:   # default is to do this
	if useStaticLibs:
		if (os_type == "Darwin"):
			lib_list.append(STATIC_NLOPT_LIBRARY_FILE_MACOSX)
			lib_list_1d.append(STATIC_NLOPT_LIBRARY_FILE_MACOSX)
		else:
			# assuming we're on a Linux system
			lib_list.append(STATIC_NLOPT_LIBRARY_FILE1_LINUX)
			lib_list_1d.append(STATIC_NLOPT_LIBRARY_FILE1_LINUX)
	else:
		lib_list.append("nlopt")	
		lib_list_1d.append("nlopt")	
else:
	extra_defines.append("NO_NLOPT")

if useOpenMP:   # default is to do this (turn this off with "--no-openmp")
	cflags_opt.append("-fopenmp")
	cflags_db.append("-fopenmp")
	link_flags.append("-fopenmp")
	extra_defines.append("USE_OPENMP")

if useExtraFuncs:   # default is to *not* do this; user must specify with "--openmp"
	extra_defines.append("USE_EXTRA_FUNCS")


if buildFatBinary and (os_type == "Darwin"):
	# note that we have to specify "-arch xxx" as "-arch", "xxx", otherwise SCons
	# passes "-arch xxx" wrapped in quotation marks, which gcc/g++ chokes on.
	cflags_opt += ["-arch", "i686", "-arch", "x86_64"]
	cflags_db += ["-arch", "i686", "-arch", "x86_64"]
	link_flags += ["-arch", "i686", "-arch", "x86_64"]

# 	if key == 'mode':
# 		if value == "export":   # "scons mode=export"  [for compiling "export" versions]
# 			# build a fat Intel (32-bit/64-bit) binary (works on Mac OS X)
# 			cflags_opt.append("-arch i386 -arch x86_64")
# 			cflags_db.append("-arch i386 -arch x86_64")
# 			link_flags.append("-arch i386 -arch x86_64")
# 			useStaticLibs = True


defines_db = defines_db + extra_defines
defines_opt = defines_opt + extra_defines




# "Environments" for compilation:
# "env_debug" is environment with debugging options turned on
# "env_opt" is an environment for optimized compiling

env_debug = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_db, LINKFLAGS=link_flags, CPPDEFINES=defines_db )
env_opt = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_opt, LINKFLAGS=link_flags, CPPDEFINES=defines_opt )


# Checks for libraries and headers -- if we're not doing scons -c:
if not env_opt.GetOption('clean'):
	conf_opt = Configure(env_opt)
	cfitsioFound = conf_opt.CheckLibWithHeader('cfitsio', 'fitsio.h', 'c')
	fftwFound = conf_opt.CheckLibWithHeader('fftw3', 'fftw3.h', 'c')
	fftwThreadsFound = conf_opt.CheckLib('fftw3_threads')
	nloptFound = conf_opt.CheckLibWithHeader('nlopt', 'nlopt.h', 'c')
	gslFound = conf_opt.CheckLib('gsl')
	libsOK = False
	if cfitsioFound and fftwFound:
		libsOK = True
	else:
		print("ERROR: Failed to find one or more required libraries and/or header files (cfitsio and/or fftw3)!")
		print("\tMake sure they are installed; if necessary, include correct path to library with --lib-path option")
		print("\tand correct path to header with --header-path option")
		exit(1)
	if useFFTWThreading and not fftwThreadsFound:
		print("ERROR: Failed to find fftw3_threading library!")
		print("\tSuggestion: include correct path to library with --lib-path option")
		print("\tOR run SCons with --no-threading option")
		exit(1)
	if useGSL and not gslFound:
		print("ERROR: Failed to find gsl library!")
		print("\tSuggestion: include correct path to library with --lib-path option")
		print("\tOR run SCons with --no-gsl option")
		exit(1)
	if useNLopt and not nloptFound:
		print("ERROR: Failed to find nlopt library!")
		print("\tSuggestion: include correct path to library with --lib-path option")
		print("\tOR run SCons with --no-nlopt option")
		exit(1)
	env_opt = conf_opt.Finish()



# We have separate lists of object names (what we want the .o files to be called) and
# source names (.cpp, .c) so that we can specify separate debugging and optimized compilations.

# Pure C code
c_obj_string = """mp_enorm statistics mersenne_twister"""
c_objs = c_obj_string.split()
c_sources = [name + ".c" for name in c_objs]


# C++ code

# ModelObject and related classes:
modelobject_obj_string = """model_object convolver"""
modelobject_objs = modelobject_obj_string.split()
modelobject_sources = [name + ".cpp" for name in modelobject_objs]

# Function objects:
functionobject_obj_string = """function_object func_gaussian func_exp func_gen-exp  
		func_sersic func_gen-sersic func_core-sersic func_broken-exp
		func_broken-exp2d func_moffat func_flatsky func_gaussian-ring 
		func_gaussian-ring2side func_edge-on-disk_n4762 func_edge-on-disk_n4762v2 
		func_edge-on-ring func_edge-on-ring2side"""
if useGSL:
	# the following modules require GSL be present
	functionobject_obj_string += " func_edge-on-disk"
	functionobject_obj_string += " integrator"
	functionobject_obj_string += " func_expdisk3d"  # requires integrator
	functionobject_obj_string += " func_brokenexpdisk3d"  # requires integrator
	functionobject_obj_string += " func_gaussianring3d"  # requires integrator
if useExtraFuncs:
	# experimental extra functions for personal testing
	functionobject_obj_string += " func_broken-exp-bar"
	if useGSL:
		functionobject_obj_string += " func_brokenexpbar3d"
		functionobject_obj_string += " func_boxytest3d"

functionobject_objs = [ FUNCTION_SUBDIR + name for name in functionobject_obj_string.split() ]
functionobject_sources = [name + ".cpp" for name in functionobject_objs]


# Base files for imfit:
imfit_base_obj_string = """commandline_parser utilities image_io levmar_fit mpfit 
		diff_evoln_fit DESolver config_file_parser add_functions print_results 
		bootstrap_errors imfit_main"""
if useNLopt:
	imfit_base_obj_string += " nmsimplex_fit"
imfit_base_objs = imfit_base_obj_string.split()
imfit_base_sources = [name + ".cpp" for name in imfit_base_objs]

# Base files for makeimage:
makeimage_base_obj_string = """commandline_parser utilities image_io config_file_parser 
			add_functions makeimage_main"""
makeimage_base_objs = makeimage_base_obj_string.split()
makeimage_base_sources = [name + ".cpp" for name in makeimage_base_objs]


# imfit: put all the object and source-code lists together
imfit_objs = imfit_base_objs + modelobject_objs + functionobject_objs + c_objs
imfit_sources = imfit_base_sources + modelobject_sources + functionobject_sources + c_sources

# makeimage: put all the object and source-code lists together
makeimage_objs = makeimage_base_objs + modelobject_objs + functionobject_objs + c_objs
makeimage_sources = makeimage_base_sources + modelobject_sources + functionobject_sources + c_sources



# Finally, define the actual targets
# specify ".do" as the suffix for "full-debug" object code
imfit_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env_debug.Program("imfit_db", imfit_dbg_objlist)
imfit_opt_objlist = [ env_opt.Object(obj, src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env_opt.Program("imfit", imfit_opt_objlist)

makeimage_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(makeimage_objs, makeimage_sources) ]
env_debug.Program("makeimage_db", makeimage_dbg_objlist)
env_opt.Program("makeimage", makeimage_sources)




