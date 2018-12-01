#!/bin/bash

D0=.5
gammaD=.5
k0=2.
beta=.1

for lin_hyp in 0 1
do
    for day in 4 5 6 9
    do
	printf "lin_hyp = %i \n day = %i \n D0 = %f, gammaD = %f, k0 = %f, beta = %f \n" $lin_hyp $day $D0 $gammaD $k0 $beta 
	python2.7 optimization_Rat05.py $lin_hyp $D0 $gammaD $k0 $beta $day
    done
done


