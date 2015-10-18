# Lists of files for use by make_distributions.py

binary_only_files = """
imfit
makeimage
"""

documentation_files = """
imfit_howto.pdf
imfit_howto.tex
"""

misc_required_files = """
COPYING.txt
DISCLAIMER
README.txt
CHANGELOG.md
"""

# header files in top-level directory
source_header_files_top = """
definitions
mersenne_twister
mp_enorm
option_struct_imfit
option_struct_makeimage
param_struct
sample_configs
statistics
"""

# header files in core/
source_header_files_core = """
add_functions
bootstrap_errors
commandline_parser
config_file_parser
convolver
downsample
estimate_memory 
image_io
model_object
oversampled_region
print_results
utilities_pub
"""

# C files in top-level directory 
source_files_c = """
mersenne_twister
mp_enorm
statistics
"""


# the following are C++ files

source_files_core = """
add_functions
bootstrap_errors
commandline_parser 
config_file_parser
convolver
downsample
estimate_memory 
image_io 
imfit_main
makeimage_main
model_object
oversampled_region
print_results
utilities 
"""

source_files_solvers ="""
levmar_fit
mpfit 
diff_evoln_fit
DESolver
nmsimplex_fit
nlopt_fit
dispatch_solver
solver_results
"""

source_files_funcobj = """
function_object 
func_gaussian 
func_exp 
func_gen-exp  
func_sersic 
func_gen-sersic 
func_core-sersic 
func_broken-exp 
func_broken-exp2d
func_king
func_king2
func_moffat 
func_flatsky 
func_gaussian-ring 
func_gaussian-ring2side
func_edge-on-disk_n4762 
func_edge-on-disk_n4762v2 
func_edge-on-ring 
func_edge-on-ring2side
func_edge-on-disk
func_brokenexpdisk3d 
func_expdisk3d 
func_gaussianring3d 
integrator
"""


example_files = """
config_exponential_ic3478_256.dat
config_sersic_ic3478_256.dat
ic3478rss_256.fits
ic3478rss_256_mask.fits
config_makeimage_moffat_psf.dat
psf_moffat_51.fits
README_examples.txt
"""

testing_scripts = """
do_imfit_tests
do_makeimage_tests
"""

python_files = """
py_startup_test.py
compare_fits_files.py
compare_imfit_printouts.py
diff_printouts.py
"""


test_files = """
config_imfit_expdisk32.dat
imfit_config_ic3478_64x64.dat
imfit_config_ic3478_64x64b.dat
imfit_config_ic3478_64x64c.dat
imfit_config_n3073.dat
config_imfit_pgc35772.dat
config_imfit_poisson.dat
config_imfit_flatsky.dat
config_3x3_flatsky.dat
config_imfit_gauss-oversample-test2.dat
config_makeimage_sersictest512_bad1.dat
config_makeimage_sersictest512_bad2.dat
config_makeimage_sersictest512_bad3.dat
config_makeimage_sersictest512_bad4.dat
config_makeimage_sersictest512_bad5.dat
config_imfit_sersictest512_badlimits1.dat
config_imfit_sersictest512_badlimits2.dat
config_imfit_sersictest512_badlimits3.dat
config_imfit_badparamline.dat
uniform_image32.fits
testimage_expdisk32.fits
testimage_poisson_lowsn20.fits
ic3478rss_64x64.fits
ic3478rss_64x64_sigma.fits
ic3478rss_64x64_variance.fits
n3073rss_small.fits
n3073rss_small_cps.fits
n3073rss_small_mask.fits
pgc35772_continuum.fits
pgc35772_mask.fits
totalmask_64x64.fits
biggertest_orig.fits
gensersictest_orig.fits
sersic+exp_orig.fits
gensersictest612_conv_cutout512.fits
flatsky_128x128.fits
testimage_3x3_nan.fits
testimage_3x3_onezero.fits
mask_for_onezero.fits
oversamp_test4.fits
psf_standard.fits
psf_oversamp.fits
imfit_textout1
imfit_textout2
imfit_textout3
imfit_textout3b
imfit_textout3c_tail
imfit_textout3d
imfit_textout3d2
imfit_textout3e
imfit_textout3e2
imfit_textout4
imfit_textout4b
imfit_textout4c
imfit_textout4d
imfit_textout4e
imfit_textout4e2
imfit_textout5a_tail
imfit_textout5b_tail
imfit_textout5c_tail
imfit_textout5d_tail
imfit_textout6
imfit_textout7a
imfit_textout7b
imfit_textout7c
imfit_textout_bad1
imfit_textout_bad2
imfit_textout_bad3
imfit_textout_bad4
imfit_textout_bad5
imfit_textout_bad6
imfit_textout_bad7
imfit_textout_bad8
imfit_textout_badnloptname
config_biggertest_4c.dat
config_makeimage_gensersic512.dat
config_makeimage_sersic+exp512.dat
config_makeimage_sersic+exp512_nosize.dat
psf_moffat_35.fits
makeimage_textout0
makeimage_textout1
makeimage_textout2
makeimage_textout3
makeimage_textout4
makeimage_textout5
makeimage_textout6
makeimage_textout7
makeimage_textout8
makeimage_textout9
makeimage_textout10
makeimage_textout11
makeimage_textout12
makeimage_textout13
makeimage_textout14
makeimage_textout15
"""
