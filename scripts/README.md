Scripts directory files include:

A. Files for performing fitting by PTID
    run_init_sims.sh - launch with model# and user name (for slurm); runs 10000 initial sims
    run_all_ptids.sh - launch with model# and user name (for slurm); fits each PTID
    run_first_params.sh - launches sims in parallel using run_first_sets.pl 
    run_first_sets.pl - draws parameters from wide distributions and launches 1 simulation
    run_next_params.sh - launches subsequent fitting rounds using run_single_set.pl
    run_single_set.pl - draws parameters from distributions based on prev fit and launches sim
    indiv_rollup.pl - rolls up log files from all parallel simulations
    shed_scoring.sh - creates a .csv score file from the composite logs
    top100_dist.pl - creates parameter distributions from top 100 fits

B. Files for performing fitting by study ARM (Both_pos, Both_neg, HIV_pos, KS_pos)
    run_arms_init.sh - launch with model# and user name (for slurm); runs 10000 initial sims
    run_all_arms.sh - launch with model# and user name (for slurm); fits stats for each ARM
    run_first_arms.sh - launches sims in parallel using run_first_sets.pl
    arm_scoring.sh - analogous to shed_scoring.sh

C. Files for running multiple sims with the "best" parameter set
    rerun_all_ptids.sh - runs each listed PTID 10 times with its best parameters distribution
    rerun_arms.sh - same only by ARM rather than PTID
    run_pdf_set.pl - pulls from the best distrib and runs 1 sim.
    distrib_runs.pl - consolidates output runs into single csv file (for plotting)

D. Files for plotting outputs
    act_boxplots.R 
    model_outputs.R
    params_by_arm.R
    ptid_params.R
    results_by_param.R
    results.R
    sim_act_boxplots_alt.R
    sim_act_boxplots.R
    sim_act_correl_alt.R
    sim_act_correl_arms.R
    sim_act_correl.R
    sim_ptid_boxplots.R
    sim_sim_correl.R
