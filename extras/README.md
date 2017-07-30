## extras

This directory currently contains code for tab auto-completions in the
Bash shell.

To enable this, copy the file `imfit_completions.bash` somewhere convenient and
add the following to your `.profile`, `.bash_profile`, or `.bashrc` file:

    source /path/to/imfit_completions.bash

You should then (assuming you restart your shell, open a new terminal
window, or something similar) be able to type "`imfit --`", press the
TAB key a couple of times, and get a listing of the command-line flags
and options; typing part of a flag/option and then pressing TAB should
complete the command. (And similarly for `imfit-mcmc` and `makeimage`.)

If your operating system makes use of something like an `/etc/bash_completion.d/`
directory and your shell startup files automatically source this, you could
copy `imfit_completions.bash` there, of course. (Assuming you have permission
to modify that directory.)
