#!/usr/bin/bash
# This script will take a input file, and break it into 4 parts

no_of_parts=4
infile=$1

no_of_lines=$(wc -l ${infile} | awk '{print $1}')
lines_per_file=$()