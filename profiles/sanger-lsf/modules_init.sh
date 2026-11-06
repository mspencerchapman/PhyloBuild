# Defines the farm's `module` command for non-login shells. LSF jobs don't inherit it from the login shell
# (it's an exported bash function), so lsf_submit.sh points BASH_ENV here and every bash in the job sources it.
module() {
    eval "$(/usr/bin/tclsh8.6 /usr/lib/x86_64-linux-gnu/modulecmd.tcl bash "$@")"
}
