#!/bin/sh
#
# This script is part of the Gamma-ray Polarimetry Simulation Toolkit
# ===================================================================
#
# Based on the X-Calibur Simulation Toolkit (XST)
# by Henric Krawczynski and Fabian Kislat
#
# A detailed manual can be found at 
#    https://github.com/ComPair/GPST
#
# For Change Log see the file gpst.C
#

# Determine where we're at
mypath=$(cd `dirname $0` && pwd)
my_so=${mypath}/gpst_C.so

# Default values for parameters
debug="false"
debug_g=""
filename=""
have_filename=0
force=""
clean="false"
large="false"

# Print version info
print_version() {
    echo '.L '${mypath}'/gpst.C+
        gpst_version()' | root -l -b
}

# Print help text
print_help() {
    echo -ne "\x1b[1m"
    print_version
    echo ""
    echo -e "Usage:\x1b[0m"
    echo "    $0 [-h] [-v] [-c] [-l] [-g] [-f|++] input_file"
    echo -e "\x1b[1mParameters:\x1b[0m"
    echo "    -h   Print this help text and exit"
    echo "    -v   Print version number and exit"
    echo "    -c   Clean. Do not print version number and date on the plot"
    echo "    -l   Large labels"
    echo "    -g   Turn on debugging (both debug output from the script"
    echo "         and debug compiler options)"
    echo "    ++   Force recompilation (same as -f)"
    echo "    -f   Force recompilation (same as ++)"
    echo "    input_file"
    echo "         File with model and simulation parameters"
    echo ""
    echo "A detailed manual can be found at"
    echo -e "    \x1b[4mhttps://github.com/ComPair/GPST\x1b[0m"
}


# Check command line options
until [ $have_filename -eq 1 ]; do
    if [ $# -lt 1 ]; then
	echo "GPST ERROR: No input file specified, not enough arguments" >&2
	exit 1
    fi
    
    case $1 in
# Exit after showing help or version info
	-h)
	    print_help
	    exit 0
	    ;;
	-v)
	    print_version
	    exit 0
	    ;;
	-c)
	    clean="true"
	    ;;
	-l)
	    large="true"
	    ;;
	-g)
	    debug="true"
	    debug_g="g"
	    ;;
	-f|\+\+)
	    force="+"
	    ;;
	*)
	    filename="$1"
	    have_filename=1
	    ;;
    esac
    shift
done

# Determine if the library needs to be recompiled
if [ -f ${my_so} ] && [ "x${force}" == "x" ]; then
    objdump --syms ${my_so} | grep "\\.debug" &>/dev/null
    have_debug_syms=$?
    if [ ${debug} == "true" ] && [ ${have_debug_syms} -ne 0 ]; then
	force="+"
	echo -e "\x1b[1mGPST INFO:\x1b[0m Debug option changed. Forcing library rebuild."
    elif [ ${debug} == "false" ] && [ ${have_debug_syms} -eq 0 ]; then
	force="+"
	echo -e "\x1b[1mGPST INFO:\x1b[0m Debug option changed. Forcing library rebuild."
    fi
fi


# Tell GPST where to look for the data
export GPST_DATA_DIR=${mypath}

# Run GPST
if [ ${debug} == "true" ]; then
    echo "GPST DEBUG: Executing root -l ${mypath}/gpst.C+${force}${debug_g}(\"${filename}\", ${clean}, ${large}, ${debug})"
fi
root -l ${mypath}"/gpst.C+${force}${debug_g}(\"${filename}\", ${clean}, ${large}, ${debug})"
