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
# To build a version with OpenMP enabled
#    $ scons --openmp <target-name>
#
# To build a version *without* FFTW threading:
#    $ scons --no-threading <target-name>
#
# To build a version with full debugging printouts:
#    $ scons define=DEBUG <target-name>
#
# To build export version (with OpenMP, "fat" binaries, all libraries statically linked):
#    $ scons --openmp --fat --static <target-name>
#

# *** EXPORT CONFIGURATIONS ***
# MacOS X fat binaries
# $ scons --openmp --static --fat

# To add one or more directories to the header or library search paths:
#    $ scons --header-path=/path/to/header/dir
# OR $ scons --header-path=/path/to/header/dir:/alt/path:/another/path
#    $ scons --lib-path=/path/to/lib/dir
# etc.


# Operating-system determination via os.uname:
# First element of tuple is basic OS; 3rd element is version number;
# 5th element is processor architecture (e.g., "i386", "sun4u", "i686")
#    os.uname()[0] = "Darwin" --> Mac OS X
#    os.uname()[0] = "SunOS" --> Solaris
#    os.uname()[0] = "Linux" --> Linux

import os

# the following is for when we want to force static linking to the GSL library
# (Change these if the locations are different on your system)
STATIC_GSL_LIBRARY_FILE_MACOSX = File("/usr/local/lib/libgsl.a")
STATIC_GSL_LIBRARY_FILE1_LINUX = File("/usr/lib/libgsl.a")
STATIC_GSL_LIBRARY_FILE2_LINUX = File("/usr/lib/libgslcblas.a")


FUNCTION_SUBDIR = "function_objects/"

os_type = os.uname()[0]

cflags_opt = ["-O2", "-g0"]
cflags_db = ["-Wall", "-g3"]

base_defines = ["ANSI"]

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
useFFTWThreading = True
useOpenMP = False
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
AddOption("--openmp", dest="useOpenMP", action="store_true", 
	default=False, help="compile with OpenMP support [LIMITED AND EXPERIMENTAL!]")

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
if GetOption("useOpenMP") is True:
	useOpenMP = True
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
		# KLUDGE (adding Pnl library)...
		#lib_list.append("pnl")
	
	# and stuff for 1D programs:
	lib_list_1d.append("gsl")
	if (os_type == "Linux"):
		lib_list.append("gslcblas")
		lib_list_1d.append("gslcblas")
else:
	extra_defines.append("NO_GSL")

if useOpenMP:   # default is to *not* do this; user must specify with "--openmp"
	cflags_opt.append("-fopenmp")
	cflags_db.append("-fopenmp")
	link_flags.append("-fopenmp")


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


# Allow use to specify extra definitions via command line
# (e.g., "scons define=NO_FFTW_THREADING mode=openmp"):
# extra_defines = ["FFTW_THREADING"]
# for key, value in ARGLIST:
# 	if key == 'define':
# 		extra_defines.append(value)
# 		if value == "NO_FFTW_THREADING":   # "scons define=NO_FFTW_THREADING"
# 			# OK, user has requested FFTW threading *not* be used
# 			lib_list.remove("fftw3_threads")
# 			if (os_type == "Linux"):
# 				lib_list.remove("pthread")
# 			extra_defines.remove("FFTW_THREADING")
# 		if value == "NO_GSL":   # "scons define=NO_GSL"
# 			# OK, user does *not* have GSL library installed
# 			defines_opt.append("NO_GSL")
# 			defines_db.append("NO_GSL")
# 			lib_list.remove("gsl")
# 			lib_list.remove("gslcblas")
# 			useGSL = False
# 		if value == "STATIC":   # "scons define=STATIC"
# 			# force compilation with static linking to libraries
# 			link_flags.append("-static")
# 	if key == 'mode':
# 		if value == "export":   # "scons mode=export"  [for compiling "export" versions]
# 			# build a fat Intel (32-bit/64-bit) binary (works on Mac OS X)
# 			cflags_opt.append("-arch i386 -arch x86_64")
# 			cflags_db.append("-arch i386 -arch x86_64")
# 			link_flags.append("-arch i386 -arch x86_64")
# 			useStaticLibs = True
# 		if value == "openmp":   # "scons mode=openmp"  [for compiling with OpenMP enabled]
# 			cflags_opt.append("-fopenmp")
# 			cflags_db.append("-fopenmp")
# 			link_flags.append("-fopenmp")
# 		if value == "profile":   # "scons mode=profile"  [for profiling the code]
# 			if (os_type == "Linux"):
# 				cflags_opt.append("-pg")
# 				cflags_db.append("-pg")
# 				link_flags.append("-pg")

defines_db = defines_db + extra_defines
defines_opt = defines_opt + extra_defines



# "Environments" for compilation:
# "env_debug" is environment with debugging options turned on
# "env_opt" is an environment for optimized compiling

env_debug = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_db, LINKFLAGS=link_flags, CPPDEFINES=defines_db )
env_opt = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_opt, LINKFLAGS=link_flags, CPPDEFINES=defines_opt )


# We have separate lists of object names (what we want the .o files to be called) and
# source names (.cpp, .c) so that we can specify separate debugging and optimized compilations.

# Pure C code
c_obj_string = """mp_enorm statistics"""
c_objs = c_obj_string.split()
c_sources = [name + ".c" for name in c_objs]


# C++ code

# ModelObject and related classes:
modelobject_obj_string = """model_object convolver"""
modelobject_objs = modelobject_obj_string.split()
modelobject_sources = [name + ".cpp" for name in modelobject_objs]

# Function objects:
functionobject_obj_string = """function_object func_gaussian func_exp func_gen-exp  
		func_sersic func_gen-sersic func_flat-exp func_broken-exp func_broken-exp2d
		func_moffat func_flatsky func_gaussian-ring func_gaussian-ring2side
		func_edge-on-disk_n4762 func_edge-on-disk_n4762v2 
		func_edge-on-ring func_edge-on-ring2side"""
if useGSL:
	# the following modules require GSL be present
	functionobject_obj_string += " func_edge-on-disk"
	functionobject_obj_string += " func_expdisk3d func_gaussianring3d integrator"
functionobject_objs = [ FUNCTION_SUBDIR + name for name in functionobject_obj_string.split() ]
functionobject_sources = [name + ".cpp" for name in functionobject_objs]


# Base files for imfit:
imfit_base_obj_string = """commandline_parser utilities image_io mpfit diff_evoln_fit DESolver
		config_file_parser add_functions print_results imfit_main"""
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




# *** Other programs (nonlinfit, psfconvolve, older stuff)

env_1d = Environment( CPPPATH=include_path, LIBS=lib_list_1d, LIBPATH=lib_path,
						CCFLAGS=cflags_db, LINKFLAGS=link_flags, CPPDEFINES=defines_db )

# mersenne_twister code is only used by profilefit
c_mersenne_obj_string = """mersenne_twister"""
c_mersenne_objs = c_mersenne_obj_string.split()
c_mersenne_sources = [name + ".c" for name in c_mersenne_objs]

# ModelObject1d and related classes:
modelobject1d_obj_string = """model_object model_object_1d"""
modelobject1d_objs = modelobject1d_obj_string.split()
modelobject1d_sources = [name + ".cpp" for name in modelobject1d_objs]

# 1D FunctionObject classes:
functionobject1d_obj_string = """function_object func1d_gaussian func1d_exp func1d_sersic 
		func1d_core-sersic func1d_broken-exp func1d_moffat func1d_delta func1d_sech 
		func1d_sech2 func1d_vdksech func1d_gaussian2side 
		func1d_n1543majmin_circbulge func1d_n1543majmin func1d_n1543majmin2"""
functionobject1d_objs = [ FUNCTION_SUBDIR + name for name in functionobject1d_obj_string.split() ]
functionobject1d_sources = [name + ".cpp" for name in functionobject1d_objs]

# Base files for profilefit:
profilefit_base_obj_string = """utilities commandline_parser mpfit diff_evoln_fit DESolver
		read_profile config_file_parser add_functions_1d print_results 
		convolver convolver1d bootstrap_errors_1d profilefit_main"""
profilefit_base_objs = profilefit_base_obj_string.split()
profilefit_base_sources = [name + ".cpp" for name in profilefit_base_objs]

# profilefit: put all the object and source-code lists together
profilefit_objs = profilefit_base_objs + modelobject1d_objs + functionobject1d_objs + c_objs + c_mersenne_objs
profilefit_sources = profilefit_base_sources + modelobject1d_sources + functionobject1d_sources + c_sources + c_mersenne_sources

# psfconvolve1d: put all the object and source-code lists together
psfconvolve1d_objs = ["psfconvolve1d_main", "commandline_parser", "utilities",
					"read_profile", "convolver1d"]
psfconvolve1d_sources = [name + ".cpp" for name in psfconvolve1d_objs]



# test_commandline: put all the object and source-code lists together
test_commandline_objs = ["test_commandline_parser", "anyoption", 
			"commandline_parser", "utilities"]
test_commandline_sources = [name + ".cpp" for name in test_commandline_objs]


# source+obj lists for older or less-used programs:
# readimage: put all the object and source-code lists together
readimage_sources = ["readimage_main.cpp", "image_io.cpp"]

# psfconvolve: put all the object and source-code lists together
psfconvolve_objs = ["psfconvolve_main", "commandline_parser", "utilities",
					"image_io", "convolver"]
psfconvolve_sources = [name + ".cpp" for name in psfconvolve_objs]

# test_parser: put all the object and source-code lists together
testparser_objs = ["test_parser", "config_file_parser", "utilities"]
testparser_sources = [name + ".cpp" for name in testparser_objs]




# profile fit is fast, so we don't really need an "optimized" version
profilefit_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(profilefit_objs, profilefit_sources) ]
env_1d.Program("profilefit", profilefit_dbg_objlist)

psfconvolve_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(psfconvolve_objs, psfconvolve_sources) ]
env_debug.Program("psfconvolve", psfconvolve_dbg_objlist)

psfconvolve1d_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(psfconvolve1d_objs, psfconvolve1d_sources) ]
env_1d.Program("psfconvolve1d", psfconvolve1d_dbg_objlist)


# timing: variation on makeimage designed to time image-generation and convolution
# Base files for timing:
timing_base_obj_string = """commandline_parser utilities image_io config_file_parser 
			add_functions timing_main"""
timing_base_objs = timing_base_obj_string.split()
timing_base_sources = [name + ".cpp" for name in timing_base_objs]

timing_sources = timing_base_sources + modelobject_sources + functionobject_sources + c_sources

env_opt.Program("timing", timing_sources)


# test harnesses, etc.:
test_commandline_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(test_commandline_objs, test_commandline_sources) ]
env_debug.Program("test_commandline", test_commandline_objlist)

# older programs
# env_opt.Program("readimage", readimage_sources)
# testparser_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(testparser_objs, testparser_sources) ]
# env_debug.Program("testparser", testparser_objlist)


# Run the tests (note that these work because we define nonexistent files
# (e.g., "azprofile_test") as targets; since the specified commands never
# create these "files", they are always "out-of-date" and always force 
# execution of the command.
#env_debug.Command("nonlinfit_test", "nonlinfit3_db", "./do_test_suite3")

