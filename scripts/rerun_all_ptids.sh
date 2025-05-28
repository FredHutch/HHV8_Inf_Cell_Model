#!/bin/bash
rm -f single_runs.log
model=$1
shift
critfile="../../data/high_shedders.crit"
for i in "$@"
do
    ptid="ptid_""$i"
    echo "For $ptid ..." >> single_runs.log
    if [ ! -d $ptid\_$model ]
    then
	mkdir $ptid\_$model
    fi
    cd $ptid\_$model
    cp ../$ptid\_$model.in .
    echo "N_runs 10" >> $ptid\_$model.in
    echo "beta_std 0" >> $ptid\_$model.in
    echo "fpos_std 0" >> $ptid\_$model.in
    echo "an_std 0" >> $ptid\_$model.in
    echo "betae_std 0" >> $ptid\_$model.in
    echo "log_p_std 0" >> $ptid\_$model.in
    echo "latent_inf_std 0" >> $ptid\_$model.in
    echo "r_std 0" >> $ptid\_$model.in
    echo "exp_days_std 0" >> $ptid\_$model.in
    echo "alpha_std 0" >> $ptid\_$model.in
    echo "kappa_std 0" >> $ptid\_$model.in
# loosen limit on log_p incase log_p_mean > log_p_high!
    echo "log_p_high 7.5" >> $ptid\_$model.in
    if [ "$model" -eq "4" ]
    then
	echo "Regions 1" >> $ptid\_$model.in
    fi

    ../../scripts/run_pdf_set.pl 1 $ptid\_$model.in $critfile
    ../../scripts/distrib_runs.pl hhv8_sim.dat1 > ../$ptid\_$model.dat1.csv
    cd ..
done

