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
                printf "Case number %i: \n" $case >> ./output/rat05hesense/log.txt
		printf "D0 = %f, gammaD = %f, beta = %f, k0 = %f \n" $D $g $b $k >> ./output/rat05hesense/log.txt 
                /bin/time -f '%e' -a -o ./output/rat05hesense/log.txt python2.7 forward_sensitivity.py $D $g $b $k $case &>> ./output/rat05hesense/log.txt
                case=$(($case + 1))
            done
        done
    done
done
