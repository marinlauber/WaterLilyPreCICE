#!/usr/bin/env bash
set -e -u

if [[ $# -gt 0 ]];
then
    echo "cleaning up old data..."
    rm -rf vtk_data/ *.pvd
    echo "running 0D-CalculiX miniLIMO example..."
    julia --project=../../../ $1 ../precice-config.xml
fi