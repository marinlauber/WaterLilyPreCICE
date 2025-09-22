#!/usr/bin/env bash
set -e -u

if [[ $# -gt 0 ]];
then
    # run the case
    ccx_preCICE -i $1 -bg -precice-participant Calculix

    # postprocess
    python post.py -F "$1".frd
fi