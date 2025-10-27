#!/usr/bin/env bash
set -e -u

echo "cleaning up old files..."
rm -f *.frd *.dat *.png *.log *.sta *.cvg *.out *.fdb *.pvd *.12d
rm -rf datp/

if [[ $# -gt 0 ]];
then
    # run the case
    ccx_preCICE -i $1 -bg -precice-participant Calculix

    # postprocess
    python post.py -F "$1".frd
fi