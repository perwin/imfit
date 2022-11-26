# Lists of files for use by make_distributions.py

binary_only_files = """
imfit
imfit-mcmc
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
cdream/LICENSE_cdream.txt
"""

extras_files = """
imfit_completions.bash
README.md
"""

python_files_for_binary_dist = """
python/imfit.py
python/imfit_funcs.py
"""

# header files in core/
source_header_files_core = """
add_functions
bootstrap_errors
commandline_parser
config_file_parser
convolver
count_cpu_cores
definitions
downsample
estimate_memory 
getimages
image_io
mersenne_twister
model_object
mp_enorm
options_base
options_imfit
options_makeimage
options_mcmc
oversampled_region
param_struct
print_results
psf_oversampling_info
sample_configs
setup_model_object
statistics
utilities_pub
"""

# header files in cdream/
source_header_files_mcmc = """
array
dream
dream_params
include/rng/GSLRng
include/rng/GSLStream
include/rng/MKLRng
include/rng/MKLStream
include/rng/Rng
include/rng/RngStream
"""
                      

# the following are C++ files

source_files_core = """
add_functions
bootstrap_errors
commandline_parser 
config_file_parser
convolver
count_cpu_cores
downsample
estimate_memory 
getimages
image_io 
imfit_main
makeimage_main
mcmc_main
mersenne_twister
model_object
mp_enorm
oversampled_region
print_results
psf_oversampling_info
setup_model_object
statistics
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

source_files_mcmc ="""
check_outliers
dream
dream_initialize
dream_pars
gelman_rubin
gen_CR
restore_state
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
func_tilted-sky-plane
func_gaussian-ring 
func_gaussian-ring2side
func_gaussian-ring-az 
func_edge-on-disk_n4762 
func_edge-on-disk_n4762v2 
func_n4608disk 
func_edge-on-ring 
func_edge-on-ring2side
func_edge-on-disk
func_brokenexpdisk3d 
func_expdisk3d 
func_gaussianring3d 
func_ferrersbar2d
func_ferrersbar3d
func_flatbar
func_pointsource
func_pointsource-rot
helper_funcs
helper_funcs_3d
integrator
psf_interpolators
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
do_mcmc_tests
do_makeimage_tests
"""

python_files = """
py_startup_test.py
check_speedup.py
compare_fits_files.py
compare_imfit_printouts.py
diff_printouts.py
imfit.py
imfit_funcs.py
"""


# for tests/imfit_reference/
test_files_imfit = """
config_imfit_expdisk32.dat
imfit_config_ic3478_64x64.dat
imfit_config_ic3478_64x64b.dat
imfit_config_ic3478_64x64c.dat
imfit_config_n3073.dat
config_imfit_pgc35772.dat
config_imfit_gauss-oversample-test2.dat
config_imfit_2gauss_small.dat
config_imfit_poisson.dat
config_imfit_flatsky.dat
config_imfit_small-gaussian.dat
config_3x3_flatsky.dat
config_imfit_ptsource.dat
config_imfit_sersictest512_badlimits1.dat
config_imfit_sersictest512_badlimits2.dat
config_imfit_sersictest512_badlimits3.dat
config_imfit_badparamline.dat
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
imfit_textout4e3
imfit_textout5a_tail
imfit_textout5b_tail
imfit_textout5c_tail
imfit_textout5d_tail
imfit_textout5e_tail
imfit_textout6
imfit_textout6b
imfit_textout6c
imfit_textout6d
imfit_textout6e
imfit_textout6f
imfit_textout6g
imfit_textout6h
imfit_textout7a
imfit_textout7b
imfit_textout7c
imfit_textout7d
imfit_textout7e
imfit_textout7f
imfit_textout8a
imfit_textout_bad0
imfit_textout_bad1
imfit_textout_bad2
imfit_textout_bad3
imfit_textout_bad4
imfit_textout_bad5
imfit_textout_bad6
imfit_textout_bad7
imfit_textout_bad8
imfit_textout_badnloptname
"""

# for tests/makeimage_reference/
test_files_makeimage = """
config_biggertest_4c.dat
config_makeimage_gensersic512.dat
config_makeimage_sersic+exp512b.dat
config_makeimage_sersic+exp512.dat
config_makeimage_sersic+exp512_nosize.dat
config_makeimage_gauss-oversample.dat
config_makeimage_gauss-oversample2.dat
config_makeimage_pointsource-oversample.dat
config_makeimage_sersictest512_bad1.dat
config_makeimage_sersictest512_bad2.dat
config_makeimage_sersictest512_bad3.dat
config_makeimage_sersictest512_bad4.dat
config_makeimage_sersictest512_bad5.dat
makeimage_textout0
makeimage_textout1
makeimage_textout2
makeimage_textout3
makeimage_textout4
makeimage_textout5
makeimage_textout5b
makeimage_textout6
makeimage_textout7
makeimage_textout8
makeimage_textout9
makeimage_textout10
makeimage_textout11
makeimage_textout12
makeimage_textout13
makeimage_textout13b
makeimage_textout13c
makeimage_textout14
makeimage_textout15
savefluxes_out.dat
"""

# for tests/imfit-mcmc_reference/
test_files_mcmc = """
config_imfit_faintstar.dat
config_imfit_faintstar_nolims.dat
mcmc_textout1
mcmc_ref1.1.txt_skip3
mcmc_ref2.1.txt_skip3
mcmc_test_short_ref.1.txt_last100
mcmc_ref4.1.txt_skip3
"""

# for tests/
test_files = """
bestfit_params_2gauss_small_tail.dat
bootstrap_output_seed10_tail.dat
bootstrap_output_seed10_2_tail.dat
uniform_image32.fits
testimage_expdisk32.fits
testimage_poisson_lowsn20.fits
testimage_3x3_ones.fits
testimage_3x3_allzeros.fits
testimage_3x3_mask-with-nan.fits
ic3478rss_64x64.fits
ic3478rss_64x64_sigma.fits
ic3478rss_64x64_variance.fits
n3073rss_small.fits
n3073rss_small_cps.fits
n3073rss_small_mask.fits
pgc35772_continuum.fits
pgc35772_mask.fits
faintstar.fits
totalmask_64x64.fits
biggertest_orig.fits
gensersictest_orig.fits
sersic+exp_orig.fits
gensersictest612_conv_cutout512.fits
testimage_2gauss_psf.fits
twogaussian_psf+2osamp_noisy.fits
flatsky_128x128.fits
testimage_3x3_nan.fits
testimage_3x3_onezero.fits
testimage_3x3_ones.fits
mask_for_onezero.fits
oversamp_test4.fits
test_emptyhdu.fits
test_multiextension_hdu0empty.fits
test_table.fits
psf_standard.fits
psf_oversamp.fits
psf_moffat_35.fits
psf_moffat_35_n4699z.fits
psf_moffat_35_oversamp3.fits
psf_moffat_fwhm2.fits
psf_moffat_fwhm2_35x35.fits
test_dump_mcmc2a
test_dump_mcmc2b
oversampled_orig.fits
oversampled2both_orig.fits
oversampled_pointsource.fits
mcmc_data/mcmc_test_short.1.txt
mcmc_data/mcmc_test_short.2.txt
mcmc_data/mcmc_test_short.3.txt
mcmc_data/mcmc_test_short.4.txt
mcmc_data/mcmc_test_short.5.txt
mcmc_data/mcmc_test_short.6.txt
"""
