#!/bin/bash

for rat_num in 01 02 05 06 09 12
do
    for day_index in 1 2 3
    do
	python all_rats_optimization.py $rat_num $day_index
    done
done
