#!/usr/bin/bash
 
# set tic levels from file at the gnuplot prompt
# examples:
# gnuplot> `./tff.sh data.dat 1 xtics`
# gnuplot> `./tff.sh data.dat 2 ytics`
 
# dataFile from which tic levels will be read
dataFile=$1
 
# which column to use
column=$2
 
# xtics OR ytics 
whichTics=$3
 
sed -e '/^#/d' $dataFile | \
    sed -e '/^$/d' | \
    cut -f $column | \
    tr '\n' ',' | \
    sed -e "s/^/set $whichTics (/" | \
    sed -e 's/,$/)\n/'
