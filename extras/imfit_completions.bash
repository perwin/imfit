# Bash tab auto-completion definitions for imfit, imfit-mcmc, and makeimage

# First, define the function to be called when the user has typed the command
# name and then pressed TAB.
# Then, register the function using "complete -F <function_name> <command_name>"
# (add -f to ensure that filename completion also works)

_imfit() 
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="--help --version --list-functions --list-parameters --sample-config --config
--noise --mask --psf --overpsf --overpsf_scale --overpsf_region
--save-params --save-model --save-residual --save-weights 
--sky --gain --readnoise --exptime --ncombined --exptime 
--errors-are-variances --errors-are-weights --mask-zero-is-bad 
--model-errors --cashstat --poisson-mlr --mlr --ftol
--nm --nlopt --de --bootstrap --save-bootstrap --chisquare-only --fitstat-only
--quiet --silent --loud --max-threads --seed --nosubsampling"

    if [[ ${cur} == -* ]] ; then
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    fi
}
complete -f -F _imfit imfit


_imfit_mcmc() 
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="--help --version --list-functions --list-parameters --sample-config --config
--noise --mask --psf --overpsf --overpsf_scale --overpsf_region
--save-params --save-model --save-residual --save-weights 
--sky --gain --readnoise --exptime --ncombined --exptime 
--errors-are-variances --errors-are-weights --mask-zero-is-bad 
--model-errors --cashstat --poisson-mlr --mlr
--output --append --nchains --max-chain-length --burnin-length --gelman-evals
--uniform-offset --gaussian-offset
--quiet --silent --loud --max-threads --seed --nosubsampling"

    if [[ ${cur} == -* ]] ; then
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    fi
}
complete -f -F _imfit_mcmc imfit-mcmc


_makeimage() 
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="--help --version --list-functions --list-parameters --sample-config --output
--refimage --psf --overpsf --overpsf_scale --overpsf_region --ncols --nrows --nosubsampling
--output-functions --print-fluxes --save-fluxes --estimation-size --zero-point --nosave 
--timing --max-threads --debug"

    if [[ ${cur} == -* ]] ; then
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    fi
}
complete -f -F _makeimage makeimage

