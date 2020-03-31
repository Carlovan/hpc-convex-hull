#!/bin/bash

inputs=$PWD/inputs/*

echo "executable,input,time,threads"

for e in $PWD/build/*; do
    baseEx=`basename $e`
    for i in $inputs; do
        baseIn=`basename $i`
        echo -n "$baseEx,$baseIn,"
        $e < $i 2>&1 > /dev/null | grep -E "Elapsed|Threads" |
            sed -E 's/^.*time: (.*)$/\1/' |
            sed -E 's/^.*Threads: (.*)$/\1/' |
            paste -sd ','
    done
done
