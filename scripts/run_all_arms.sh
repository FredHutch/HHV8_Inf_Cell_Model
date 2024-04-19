#!/bin/bash
#rm arm_runs.log
model=$1
shift
user=$1
shift

echo "For model $model ..." >> arm_runs.log


for arm in "$@"
do
    echo "For $arm ..." >> arm_runs.log
    cp arms.in $arm\_$model.in

    ../scripts/arm_scoring.sh $arm 0 $model init_$model.txt
    high=`cut -f82 -d',' top100_$arm\_$model.csv|head -n 2|tail -n 1`
    low=`cut -f82 -d',' top100_$arm\_$model.csv|tail -n 1`
    echo "Top 100 scores for $arm model=$model initially range from $high to $low" >> arm_runs.log
    cp top100_$arm\_$model.csv init100_$arm\_$model.csv
done

for arm in "$@"
do
    echo "For $arm ..." >> arm_runs.log

    attempts=2
    while [ $attempts -gt 0 ]
    do
        echo -n "Launching next set of runs at " >> arm_runs.log
	date >> arm_runs.log
	../scripts/run_next_params.sh $arm $model 1000
	sleep 30
	running=`squeue -u $user | wc -l`
	while [ $running -gt 5 ]
	do
	    running=`expr $running - 1`
	    echo -n "Waiting for $running runs to complete at " >> arm_runs.log
	    date >> arm_runs.log
	    sleep 30
	    running=`squeue -u $user | wc -l`
	done
	if [ $running -gt 1 ]
	then
	    running=`expr $running - 1`
	    echo -n "Killed $running runs at " >> arm_runs.log
	    date >> arm_runs.log
	    scancel -u $user
	fi
	cat `grep -l "Total Score" set_*/slu*` > $arm\_$model.txt
	rm -r set_*/
	../scripts/arm_scoring.sh $arm 0 $model $arm\_$model.txt
	high=`cut -f82 -d',' top100_$arm\_$model.csv|head -n 2|tail -n 1`
	low=`cut -f82 -d',' top100_$arm\_$model.csv|tail -n 1`
	echo "Top 100 scores for $arm model=$model now range from $high to $low" >> arm_runs.log
	attempts=`expr $attempts - 1`
    done
done

