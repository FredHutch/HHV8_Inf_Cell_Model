#!/bin/bash
rm init_arms_runs.log
model=$1
shift
user=$1
shift

echo "For model $model ..." >> init_arms_runs.log
cp arms.in arms_$model.in

echo -n "Launching 1st set of runs at " >> init_arms_runs.log
date >> init_arms_runs.log
../scripts/run_first_arms.sh $model 10000
sleep 30
running=`squeue -u $user | wc -l`
while [ $running -gt 5 ]
do
    running=`expr $running - 1`
    echo -n "Waiting for $running initial runs to complete at " >> init_arms_runs.log
    date >> init_arms_runs.log
    sleep 30
    running=`squeue -u $user | wc -l`
done
if [ $running -gt 1 ]
then
    running=`expr $running - 1`
    echo -n "Killed $running initial runs at " >> init_arms_runs.log
    date >> init_arms_runs.log
    scancel -u $user
fi
cat `grep -l "Total Score" set_*/slu*` > init_arms_$model.txt
rm -r set_*/
