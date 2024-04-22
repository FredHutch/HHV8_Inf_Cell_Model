#!/bin/bash
#rm ptid_runs.log
model=$1
shift

echo "For model $model ..." >> ptid_runs.log
cp ptid.in ptid_$model.in

for i in "$@"
do
    ptid="ptid_""$i"
    echo "For $ptid ..." >> ptid_runs.log
    cp ptid.in $ptid\_$model.in

    ../scripts/shed_scoring.sh $ptid $i $model init_$model.txt
    high=`cut -f82 -d',' top100_$ptid\_$model.csv|head -n 2|tail -n 1`
    low=`cut -f82 -d',' top100_$ptid\_$model.csv|tail -n 1`
    echo "Top 100 scores for $ptid model=$model initially range from $high to $low" >> ptid_runs.log
    cp top100_$ptid\_$model.csv init100_$ptid\_$model.csv
done
rm -r set_*/

for i in "$@"
do
    ptid="ptid_""$i"
    echo "For $ptid ..." >> ptid_runs.log
    attempts=2
    while [ $attempts -gt 0 ]
    do
        echo -n "Launching next set of runs at " >> ptid_runs.log
	date >> ptid_runs.log
	../scripts/run_next_params.sh $ptid $model 1000
	sleep 30
	running=`squeue -u dswan | wc -l`
	while [ $running -gt 5 ]
	do
	    running=`expr $running - 1`
	    echo -n "Waiting for $running runs to complete at " >> ptid_runs.log
	    date >> ptid_runs.log
	    sleep 30
	    running=`squeue -u dswan | wc -l`
	done
	if [ $running -gt 1 ]
	then
	    running=`expr $running - 1`
	    echo -n "Killed $running runs at " >> ptid_runs.log
	    date >> ptid_runs.log
	    scancel -u dswan
	fi
	cat `grep -l "Total Score" set_*/slu*` > $ptid\_$model.txt
	rm -r set_*/
	../scripts/shed_scoring.sh $ptid $i $model $ptid\_$model.txt
	high=`cut -f82 -d',' top100_$ptid\_$model.csv|head -n 2|tail -n 1`
	low=`cut -f82 -d',' top100_$ptid\_$model.csv|tail -n 1`
	echo "Top 100 scores for $ptid model=$model now range from $high to $low" >> ptid_runs.log
	attempts=`expr $attempts - 1`
    done
done

