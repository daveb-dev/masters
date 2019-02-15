#!/bin/bash

rat_num=(01 02 05 06 09 12)
for i in ${!rat_num[@]}; 
do
    for day_index in 1 2 3
    do
	#echo "$i, ${rat_num[$i]}, $day_index"
	python2.7 all_rats_optimization.py ${rat_num[$i]} $i $day_index
    done
done
