#!/bin/sh
PWD=`pwd`
export PYTHONPATH="$PYTHONPATH:$PWD/src/:$PWD/test"

python -m unittest discover -v
