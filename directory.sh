#!/bin/bash
# script to check if data directory already exists
if [ ! -d outputs ]; then mkdir outputs; fi
cd outputs
[ -d $1 ] && echo "Data directory exists" || mkdir $1
cd $1
pwd>../path.txt
