#!/bin/bash

lehe=0
#python2.7 forward_sensitivity.py $lehe .1 .1 .1 1 1

# Run full number of cases and time
#: <<'EOF'
tmp=./tmp.txt
logfile=./output/sensitivity/le_log.txt
> $logfile
case=1
for D in .01 .1 1 10 
do
    for g in .01 .1 1 10
    do
        for b in .01 .1 1 10
        do
            for k in .01 .1 1 10
            do
                > $tmp
		printf "Case number %i: \n" $case >> $logfile
		printf "D0 = %f, gammaD = %f, beta = %f, k = %f \n" $D $g $b $k >> $logfile
                /bin/time -f '%e' -a -o $logfile python2.7 forward_sensitivity.py $lehe $D $g $k $b $case &>> $tmp
		num_lines=$(cat $tmp | wc -l)
		if [ $num_lines -ne 0 ]
		then
		    printf "FAILED \n" >> $logfile
		fi
                case=$(($case + 1))
            done
        done
    done
done
#EOF

# Debugging reduced number of cases
: <<'END'
logfile=./output/sensitivity/le_log.txt
> $logfile
case=1
D=.01
g=.01
for b in .01
do
    for k in .01
    do
	printf "Case number %i: \n" $case >> $logfile
	printf "D0 = %f, gammaD = %f, k0 = %f, beta = %f \n" $D $g $k $b >> $logfile
	/bin/time -f '%e' -a -o $logfile python2.7 forward_sensitivity.py $lehe $D $g $k $b $case &>> $tmp
	num_lines=$(cat $tmp | wc -l)
	echo $num_lines
    done
done
END
