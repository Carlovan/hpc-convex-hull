#!/bin/bash

inputs=$PWD/inputs/*
NPROC=`nproc`
EVALS=5

echo "executable,input,time,threads"

for ex in $PWD/build/*; do
    baseEx=`basename $ex`
    for inFile in $inputs; do
        baseIn=`basename $inFile`
        maxThreads=1
        if [[ $baseEx = omp-* || $baseEx = mpi-* ]]; then
            maxThreads=$NPROC
        fi
        for threads in `seq 1 $maxThreads`; do
            if [[ $baseEx = omp-* ]]; then
                export OMP_NUM_THREADS=$threads
            fi
            if [[ $baseEx = mpi-* ]]; then
                ex="mpirun -n $threads $ex"
            fi

            echo -n "$baseEx,$baseIn,"
            for _ in `seq 1 $((EVALS + 2))`; do
                $ex < $inFile 2>&1 > /dev/null | grep -E "Elapsed" | sed -E 's/^.*time: (.*)$/\1/'
            done | sort -n | tail -n +2 | head -n -1 | awk '{x+=$1; next} END{printf "%.7f", x/NR}' | tr -d '\n'
            echo ",$threads"
        done
    done
done
