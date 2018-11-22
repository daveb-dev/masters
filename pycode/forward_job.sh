#!/bin/bash

case=0
for D in .01 .1 1 10 
do
    for g in .01 .1 1 10
    do
        for b in .01 .1 1 10
        do
            for k in .01 .1 1 10
            do
                printf "Case number %i: \n" $case >> ./output/rat05hesense/info.log
                /bin/time -f '%e' -a -o ./output/rat05hesense/info.log python forward_sensitivity.py $D $g $b $k $case 
                case=$(($case + 1))
            done
        done
    done
done
