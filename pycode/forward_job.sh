#!/bin/bash

logfile=./output/rat05lesense/log.txt
case=1
for D in .01 .1 1 10 
do
    for g in .01 .1 1 10
    do
        for b in .01 .1 1 10
        do
            for k in .01 .1 1 10
            do
                printf "Case number %i: \n" $case >> $logfile
		printf "D0 = %f, gammaD = %f, beta = %f, k0 = %f \n" $D $g $b $k >> $logfile
                /bin/time -f '%e' -a -o $logfile python2.7 forward_sensitivity.py $D $g $b $k $case &>> $logfile
                case=$(($case + 1))
            done
        done
    done
done
