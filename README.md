# HHV8_Inf_Cell_Model
Mathematical model for HHV8 viral spread

Sub-directories are:

    include - header files
    src - C++ source files
    Linux - build/run directory (Linux-based)
    OSX - build/run directory (Mac-based)
    scripts - shell and perl scripts for launching runs, scoring fits, etc.
    data - PTID-specific criteria files (used for fitting)
    inputs - input files and PTID lists (all and by ARM)

The simulation executable, hhv8_sim, can be built and run in either the Linux or OSX (Mac) sub-directory.
It is designed to be run with an input file that specifies parameter values or a distribution (mean & stddev) for them.  Originally, the simulation outputs were scored against values from a criteria file containing summary statistics from multiple participants using binned levels of shedding, episode peaks, etc.  The scores were written to a log file for post-processing into score files (csv) by the appropriate perl script from the scripts sub-directory.

To score runs against data from an individual, the simulation must be run with parameters drawn for that individual and the output of the simuation saved for special post-processing by a separate script (indiv_rollup.pl) launched in conjunction with that participant's criteria file.  The criteria file in this case contains shedding percentage, mean and peak log10 viral levels for that individual from the 28-day observation period.  All criteria files are contained in the data sub-directory.

Several "models" were evaluated when fitting to participant data.  A model in this case is defined as a set of parameters to be fitted and their associated mechanistic effects.  The best model from our testing was one that varied infectivity of a cell (beta), the death rate of productively infected cells (an), the killing rate of each T cell (fpos) and the reactivation rate from latency (latent_inf). It is identified as model 8 in the scripts.
