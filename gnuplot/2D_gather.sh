#!/bin/bash
varfile=$1
nproc=$2
for i in `seq -f %02g $1 $1`;
do
  varfile=$i
done
fout=2D_slice$varfile.log
touch $fout
for i in `seq -f %02g 0 $(($nproc-1))`;
do 
  filename=../data/proc$i/2D_slice$varfile.log
  echo $filename
  cat $filename >> $fout
done
