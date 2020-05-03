#!/bin/bash

inputs=$PWD/inputs/*
executables=$PWD/build/*
NPROC=`nproc`
EVALS=10

export OMP_PROC_BIND=true
echo "executable,input,time,threads"

for _ in `seq 1 $EVALS`; do
for ex in $executables; do
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
            cmd=$ex
            if [[ $baseEx = mpi-* ]]; then
                cmd="mpirun -n $threads $cmd"
            fi

            echo -n "$baseEx,$baseIn,"
                $cmd < $inFile 2>&1 > /dev/null | grep -E "Elapsed" | sed -E 's/^.*time: +(.*)$/\1/' | tr -d '\n'
            echo ",$threads"
        done
    done
done
done
