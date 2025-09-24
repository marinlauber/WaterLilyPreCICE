#!/usr/bin/env bash
set -e -u

if [[ $# -gt 0 ]];
then
    julia --project=../../../ $1 ../precice-config.xml
fi