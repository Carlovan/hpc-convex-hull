#!/bin/bash

inputs=$PWD/inputs/*

echo "executable,input,time"

for e in $PWD/build/*; do
    baseEx=`basename $e`
    for i in $inputs; do
        baseIn=`basename $i`
        echo -n "$baseEx,$baseIn,"
        $e < $i 2>&1 > /dev/null | grep -E "Elapsed" |
            sed -E 's/^.*time: (.*)$/\1/' |
            paste -sd ','
    done
done
