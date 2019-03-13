#!/bin/bash

for lin_hyp in 0 1
do
    if [ $lin_hyp -eq 0 ]
    then
	gammaD=0.212689163148
	D0=0.831665484255
	beta=0.01
	k0=1.
    else
	gammaD=1.9918771142
	D0=1.01849141246
	beta=0.034372054592
	k0=1.
    fi
    
    for day in 4 5 6 9
    do
	printf "lin_hyp = %i \n day = %i \n D0 = %f, gammaD = %f, k0 = %f, beta = %f \n" $lin_hyp $day $D0 $gammaD $k0 $beta 
	python2.7 optimization_Rat05.py $lin_hyp $D0 $gammaD $k0 $beta $day
    done
done


