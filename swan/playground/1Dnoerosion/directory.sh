#!/bin/bash
# script to check if data directory already exists
[ -d $1 ] && echo "Directory Exists" || mkdir $1
cd $1
pwd>../path.txt