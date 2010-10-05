# Scons "makefile" for imfit and related programs
# 
# To use this, type the following on the command line:
#    $ scons <target-name>
# where <target-name> is, e.g., "imft", "imfit_db", etc.
# "scons" by itself will build *all* targets
# 
# To build a version with FFTW3 threading enabled:
#    $ scons define=FFTW_THREADING <target-name>
#
# To build a version with full debugging printouts:
#    $ scons define=DEBUG <target-name>
#
# To clean up files associated with a given target:
#    $ scons -c <target-name>
# To clean up all targets (including programs):
#    $ scons -c


# Operating-system determination via os.uname:
# First element of tuple is basic OS; 3rd element is version number;
# 5th element is processor architecture (e.g., "i386", "sun4u", "i686")
#    os.uname()[0] = "Darwin" --> Mac OS X
#    os.uname()[0] = "SunOS" --> Solaris
#    os.uname()[0] = "Linux" --> Linux

import os

FUNCTION_SUBDIR = "function_objects/"

os_type = os.uname()[0]

cflags_opt = ["-O2", "-g0"]
cflags_db = ["-Wall", "-g3"]

base_defines = ["ANSI"]

lib_list = ["fftw3", "cfitsio", "m", "gsl"]

# system-specific stuff
if (os_type == "Darwin"):   # OK, we're compiling on Mac OS X
	# /sw/include and /sw/lib are for GNU Scientific Library (gsl) stuff
	include_path = ["/usr/local/include", "/sw/include", FUNCTION_SUBDIR]
	lib_path = ["/usr/local/lib", "/sw/lib"]
	# Use "-m32" for both compiling and linking, to make sure we work in 
	# 32-bit mode for Mac OS X (10.6 defaults to 64-bit compilation otherwise,
	# which would cause problems when linking with fftw3 and cfitsio libraries):
	cflags_opt.append("-m32")
	cflags_db = ["-Wall", "-Wshadow", "-Wredundant-decls", "-Wpointer-arith", "-g3", "-m32"]
	link_flags = ["-m32"]
if (os_type == "Linux"):
	# When compiled under Linux, -O3 causes mysterious "invalid pointer" error at end of run
	cflags_opt = ["-O3", "-g0"]
	cflags_db = ["-Wall", "-Wshadow", "-Wredundant-decls", "-Wpointer-arith", "-g3"]
	# silly Linux doesn't have OpenBSD string routines built in, so we'll have to include them
	base_defines = base_defines + ["NO_BSD_STRINGS"]
	include_path = ["/home/erwin/include", "/usr/include", FUNCTION_SUBDIR]
	lib_path = ["/home/erwin/lib"]
	lib_list.append("gslcblas")
	link_flags = None
defines_opt = base_defines
#defines_db = base_defines + ["DEBUG"]
defines_db = base_defines


# Allow use to specify extra definitions via command line
# (e.g., "scons define=FFTW_THREADING"):
extra_defines = []
for key, value in ARGLIST:
	if key == 'define':
		extra_defines.append(value)
		if value == "FFTW_THREADING":
			# OK, user has requested FFTW threading to be turned on
			lib_list.insert(0, "fftw3_threads")
			if (os_type == "Linux"):
				lib_list.append("pthread")

defines_db = defines_db + extra_defines
defines_opt = defines_opt + extra_defines



# "env_debug" is environment with debugging options turned on
# "env_opt" is an environment for optimized compiling [not really used yet]

env_debug = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_db, LINKFLAGS=link_flags, CPPDEFINES=defines_db )
# env_norm = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
# 						CCFLAGS=cflags_opt, LINKFLAGS=link_flags, CPPDEFINES=defines_opt )
env_opt = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_opt, LINKFLAGS=link_flags, CPPDEFINES=defines_opt )

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
		func_sersic func_gen-sersic func_flat-exp func_broken-exp func_broken-exp2d
		func_edge-on-disk func_moffat func_flatsky
		func_edge-on-disk_n4762 func_edge-on-disk_n4762v2 
		func_edge-on-ring func_edge-on-ring2side"""
functionobject_objs = [ FUNCTION_SUBDIR + name for name in functionobject_obj_string.split() ]
functionobject_sources = [name + ".cpp" for name in functionobject_objs]

# ModelObject1d and related classes:
modelobject1d_obj_string = """model_object model_object_1d"""
modelobject1d_objs = modelobject1d_obj_string.split()
modelobject1d_sources = [name + ".cpp" for name in modelobject1d_objs]

# 1D FunctionObject classes:
functionobject1d_obj_string = """function_object func1d_gaussian func1d_exp func1d_sersic 
		func1d_broken-exp func1d_moffat func1d_delta func1d_sech func1d_sech2 func1d_vdksech"""
functionobject1d_objs = [ FUNCTION_SUBDIR + name for name in functionobject1d_obj_string.split() ]
functionobject1d_sources = [name + ".cpp" for name in functionobject1d_objs]


# Base files for imfit:
imfit_base_obj_string = """utilities anyoption image_io mpfit diff_evoln_fit DESolver
		config_file_parser add_functions print_results imfit_main"""
imfit_base_objs = imfit_base_obj_string.split()
imfit_base_sources = [name + ".cpp" for name in imfit_base_objs]

# Base files for profilefit:
profilefit_base_obj_string = """utilities anyoption mpfit diff_evoln_fit DESolver
		read_profile config_file_parser add_functions_1d print_results 
		convolver convolver1d profilefit_main"""
profilefit_base_objs = profilefit_base_obj_string.split()
profilefit_base_sources = [name + ".cpp" for name in profilefit_base_objs]

# Base files for makeimage:
makeimage_base_obj_string = """anyoption utilities image_io config_file_parser 
			add_functions makeimage_main"""
makeimage_base_objs = makeimage_base_obj_string.split()
makeimage_base_sources = [name + ".cpp" for name in makeimage_base_objs]

# Base files for timing:
timing_base_obj_string = """anyoption utilities image_io config_file_parser 
			add_functions timing_main"""
timing_base_objs = timing_base_obj_string.split()
timing_base_sources = [name + ".cpp" for name in timing_base_objs]


# imfit: put all the object and source-code lists together
imfit_objs = imfit_base_objs + modelobject_objs + functionobject_objs + c_objs
imfit_sources = imfit_base_sources + modelobject_sources + functionobject_sources + c_sources

# profilefit: put all the object and source-code lists together
profilefit_objs = profilefit_base_objs + modelobject1d_objs + functionobject1d_objs + c_objs
profilefit_sources = profilefit_base_sources + modelobject1d_sources + functionobject1d_sources + c_sources

# makeimage: put all the object and source-code lists together
makeimage_objs = makeimage_base_objs + modelobject_objs + functionobject_objs + c_objs
makeimage_sources = makeimage_base_sources + modelobject_sources + functionobject_sources + c_sources
# timing: variation on makeimage designed to time image-generation and convolution
timing_sources = timing_base_sources + modelobject_sources + functionobject_sources + c_sources


# readimage: put all the object and source-code lists together
#readimage_sources = ["readimage_main.cpp", "readimage.cpp"]
readimage_sources = ["readimage_main.cpp", "image_io.cpp"]

# psfconvolve: put all the object and source-code lists together
#psfconvolve_sources_old = ["psfconvolve_main_old.cpp", "anyoption.cpp", "image_io.cpp"]
psfconvolve_objs = ["psfconvolve_main", "anyoption", "image_io", "convolver"]
psfconvolve_sources = [name + ".cpp" for name in psfconvolve_objs]

# psfconvolve1d: put all the object and source-code lists together
psfconvolve1d_objs = ["psfconvolve1d_main", "anyoption", "read_profile", "convolver1d"]
psfconvolve1d_sources = [name + ".cpp" for name in psfconvolve1d_objs]

# test_parser: put all the object and source-code lists together
testparser_objs = ["test_parser", "config_file_parser", "utilities"]
testparser_sources = [name + ".cpp" for name in testparser_objs]

# test_commandline: put all the object and source-code lists together
test_commandline_objs = ["test_commandline_parser", "anyoption", "utilities"]
test_commandline_sources = [name + ".cpp" for name in test_commandline_objs]


# Finally, define the actual targets
# specify ".do" as the suffix for "full-debug" object code
imfit_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env_debug.Program("imfit_db", imfit_dbg_objlist)
imfit_opt_objlist = [ env_opt.Object(obj, src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env_opt.Program("imfit", imfit_opt_objlist)

# profile fit is fast, so we don't really need an "optimized" version
profilefit_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(profilefit_objs, profilefit_sources) ]
env_debug.Program("profilefit", profilefit_dbg_objlist)

makeimage_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(makeimage_objs, makeimage_sources) ]
env_debug.Program("makeimage_db", makeimage_dbg_objlist)
env_opt.Program("makeimage", makeimage_sources)

psfconvolve_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(psfconvolve_objs, psfconvolve_sources) ]
env_debug.Program("psfconvolve", psfconvolve_dbg_objlist)

psfconvolve1d_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(psfconvolve1d_objs, psfconvolve1d_sources) ]
env_debug.Program("psfconvolve1d", psfconvolve1d_dbg_objlist)

env_opt.Program("timing", timing_sources)

# test harnesses, etc.:
test_commandline_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(test_commandline_objs, test_commandline_sources) ]
env_debug.Program("test_commandline", test_commandline_objlist)

# older programs
env_opt.Program("readimage", readimage_sources)
testparser_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(testparser_objs, testparser_sources) ]
env_debug.Program("testparser", testparser_objlist)


# Run the tests (note that these work because we define nonexistent files
# (e.g., "azprofile_test") as targets; since the specified commands never
# create these "files", they are always "out-of-date" and always force 
# execution of the command.
#env_debug.Command("nonlinfit_test", "nonlinfit3_db", "./do_test_suite3")

