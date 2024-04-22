#!/bin/bash
rm -f arm_reruns.log
model=$1
shift
for arm in "$@"
do
    echo "For $arm ..." >> arm_reruns.log
    if [ ! -d $arm\_$model ]
    then
	mkdir $arm\_$model
    fi
    cd $arm\_$model
    cp ../$arm\_$model.in .
    echo "N_runs 76" >> $arm\_$model.in

    ../../scripts/run_pdf_set.pl 1 $arm\_$model.in ../../data/$arm.crit
    ../../scripts/distrib_runs.pl hhv8_sim.dat1 > ../$arm\_$model.dat1.csv
    cd ..
done

