#!/bin/bash
INPUTS=$PWD/inputs/*
EXE=$1
OUTDIR=`basename $EXE`_output

mkdir -p $OUTDIR
for i in $INPUTS; do
    BASEI=`basename $i`
    OUTFILE=$OUTDIR/${BASEI%.in}.out
    $EXE < $i > $OUTFILE
done
