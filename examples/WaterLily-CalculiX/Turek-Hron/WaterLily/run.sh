#!/usr/bin/env bash
set -e -u

if [[ $# -gt 0 ]];
then
    # run the case
    julia --project=../../../ $1 ../precice-config.xml
fi