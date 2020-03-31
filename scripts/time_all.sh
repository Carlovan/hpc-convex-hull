#!/bin/bash

inputs=$PWD/inputs/*
NPROC=`nproc`

echo "executable,input,time,threads"

for e in $PWD/build/*; do
    baseEx=`basename $e`
    for i in $inputs; do
        baseIn=`basename $i`
        maxThreads=1
        if [[ $baseEx = omp-* ]]; then
            maxThreads=$NPROC
        fi
        for threads in `seq 1 $maxThreads`; do
            echo -n "$baseEx,$baseIn,"
            OMP_NUM_THREADS=$threads $e < $i 2>&1 > /dev/null | grep -E "Elapsed|Threads" |
                sed -E 's/^.*time: (.*)$/\1/' |
                sed -E 's/^.*Threads: (.*)$/\1/' |
                paste -sd ','
        done
    done
done
