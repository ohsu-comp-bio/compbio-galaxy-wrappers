#!/bin/sh
PWD=`pwd`
export PYTHONPATH="$PYTHONPATH:$PWD/src/"
cd test
python -m unittest discover -v
