# Scons "makefile" for imfit program
# 
# To use this, type the following on the command line:
#    $ scons <target-name>
# where <target-name> is, e.g., "imft", "imfit_db", etc.
# "scons" by itself will build *all* targets
# 
# To clean up everything (including programs):
#    $ scons -c


# Operating-system determination via os.uname:
# First element of tuple is basic OS; 3rd element is version number;
# 5th element is processor architecture (e.g., "i386", "sun4u", "i686")
#    os.uname()[0] = "Darwin" --> Mac OS X
#    os.uname()[0] = "SunOS" --> Solaris
#    os.uname()[0] = "Linux" --> Linux

import os

os_type = os.uname()[0]

cflags_opt = ["-O3", "-g0"]
cflags_db = ["-Wall", "-g3"]

base_defines = ["ANSI"]
# slightly experimental system-specific stuff
if (os_type == "Darwin"):
	# OK, we're compiling on Mac OS X
	include_path = ["/usr/local/include"]
	lib_path = ["/usr/local/lib"]
	# Use "-m32" for both compiling and linking, to make sure we work in 
	# 32-bit mode for Mac OS X (10.6 defaults to 64-bit compilation otherwise,
	# which would cause problems when linking with fftw3 library):
	cflags_opt = ["-O3", "-g0", "-m32"]
	cflags_db = ["-Wall", "-g3", "-m32"]
	link_flags = ["-m32"]
if (os_type == "Linux"):
	# silly Linux doesn't have OpenBSD string routines built in, so we'll
	# have to include them
	base_defines = base_defines + ["NO_BSD_STRINGS"]
	include_path = ["/usr/include"]
	lib_path = ["/usr/lib"]
	link_flags = None
defines_opt = base_defines
defines_db = base_defines + ["DEBUG"]


lib_list = ["fftw3", "cfitsio", "m"]


# "env_norm" is currently a semi-debugging environment (most debugging options
# 			turned on, but not -DDEBUG
# "env_debug" is identical, except that -DDEBUG is also specified
# "env_opt" is an environment for optimized compiling [not really used yet]

env_debug = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_db, LINKFLAGS=link_flags, CPPDEFINES=defines_db )
env_norm = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_opt, LINKFLAGS=link_flags, CPPDEFINES=defines_opt )
env_opt = Environment( CPPPATH=include_path, LIBS=lib_list, LIBPATH=lib_path,
						CCFLAGS=cflags_opt, LINKFLAGS=link_flags, CPPDEFINES=defines_opt )

# We have separate lists of object names (what we want the .o files to be 
# called) and source names (.cpp, .c) so that we can specify separate 
# debugging and optimized compilations.

# Pure C code
c_obj_string = """mp_enorm statistics mersenne_twister"""
c_objs = c_obj_string.split()
c_sources = [name + ".c" for name in c_objs]

# C++ code
# ModelObject and various FunctionObject classes:
modelobject_obj_string = """model_object function_object func_gaussian func_exp
		func_sersic func_flat-exp func_broken-exp"""
modelobject_objs = modelobject_obj_string.split()
modelobject_sources = [name + ".cpp" for name in modelobject_objs]

# ModelObject1d and various 1D FunctionObject classes:
modelobject1d_obj_string = """model_object model_object_1d function_object
		func1d_exp func1d_sersic"""
modelobject1d_objs = modelobject1d_obj_string.split()
modelobject1d_sources = [name + ".cpp" for name in modelobject1d_objs]

# Base files for imfit:
imfit_base_obj_string = """utilities anyoption image_io mpfit diff_evoln_fit DESolver
		config_file_parser add_functions print_results imfit_main"""
imfit_base_objs = imfit_base_obj_string.split()
imfit_base_sources = [name + ".cpp" for name in imfit_base_objs]

# Base files for profilefit:
profilefit_base_obj_string = """utilities anyoption mpfit diff_evoln_fit DESolver
		read_profile config_file_parser add_functions_1d print_results profilefit_main"""
profilefit_base_objs = profilefit_base_obj_string.split()
profilefit_base_sources = [name + ".cpp" for name in profilefit_base_objs]

# Base files for makeimage:
makeimage_base_obj_string = """anyoption utilities image_io config_file_parser 
			add_functions makeimage_main"""
makeimage_base_objs = makeimage_base_obj_string.split()
makeimage_base_sources = [name + ".cpp" for name in makeimage_base_objs]


# imfit: put all the object and source-code lists together
imfit_objs = imfit_base_objs + modelobject_objs + c_objs
imfit_sources = imfit_base_sources + modelobject_sources + c_sources

# profilefit: put all the object and source-code lists together
profilefit_objs = profilefit_base_objs + modelobject1d_objs + c_objs
profilefit_sources = profilefit_base_sources + modelobject1d_sources + c_sources

# makeimage: put all the object and source-code lists together
makeimage_objs = makeimage_base_objs + modelobject_objs + c_objs
makeimage_sources = makeimage_base_sources + modelobject_sources + c_sources

# readimage: put all the object and source-code lists together
readimage_sources = ["readimage_main.cpp", "readimage.cpp"]

# psfconvolve: put all the object and source-code lists together
psfconvolve_sources_old = ["psfconvolve_main_old.cpp", "anyoption.cpp", "image_io.cpp"]
psfconvolve_sources = ["psfconvolve_main.cpp", "anyoption.cpp", "image_io.cpp",
	"convolver.cpp"]

# test_parser: put all the object and source-code lists together
testparser_sources = ["test_parser.cpp", "config_file_parser.cpp", "utilities.cpp"]



# Finally, define the actual targets
# specify ".do" as the suffix for "full-debug" object code
full_dbg_objlist = [ env_debug.Object(obj + ".do", src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env_debug.Program("imfit_db", full_dbg_objlist)

base_norm_objlist = [ env_norm.Object(obj, src) for (obj,src) in zip(imfit_objs, imfit_sources) ]
env_norm.Program("imfit", base_norm_objlist)

env_norm.Program("profilefit", profilefit_sources)
env_norm.Program("makeimage", makeimage_sources)
env_norm.Program("readimage_test", readimage_sources)
env_norm.Program("psfconvolve", psfconvolve_sources)
env_norm.Program("psfconvolve_old", psfconvolve_sources_old)
env_norm.Program("testparser", testparser_sources)


# Specify that debug-compilation produces object files with "_debug" 
# at end of root name by making a list of scons Objects which pair
# "<rootname>_debug" with "<rootname.cpp>"
# nonlinfit_debug_Objects = [ env_debug.Object(objname + "_debug", sourcename) 
# 							for objname,sourcename in zip(nonlinfit_objs, nonlinfit_sources) ]
# env_debug.Program("nonlinfit3_db", nonlinfit_debug_Objects)


# Run the tests (note that these work because we define nonexistent files
# (e.g., "azprofile_test") as targets; since the specified commands never
# create these "files", they are always "out-of-date" and always force 
# execution of the command.
#env_debug.Command("nonlinfit_test", "nonlinfit3_db", "./do_test_suite3")

