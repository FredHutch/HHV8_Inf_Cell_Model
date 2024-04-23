# HHV8_Inf_Cell_Model
Mathematical model for HHV8 viral spread

Sub-directories are:

    include - header files
    src - C++ source files
    Linux - build/run directory (Linux-based)
    OSX - build/run directory (Mac-based)
    scripts - shell and perl scripts for launching runs, scoring fits, etc.
    data - PTID-specific criteria files (used for fitting)

The simulation executable, hhv8_sim, can be built and run in either the Linux or OSX (Mac) sub-directory.
It is designed to be run with an input file that specifies parameter values or their distribution (mean & stddev).  Originally, the simulation outputs were scored against values from a criteria file containing summary statistics from multiple participants using binned levels of shedding, episode peaks, etc.  The scores are written to a log file for post-processing into score files (csv) by the appropriate perl script from the scripts subdirectory.

To score runs against data from an individual, the simulation must be run with parameters drawn for that individual and the output of the simuation saved for special post-processing by a separate script (indiv_rollup.pl) launched in conjunction with that participan's criteria file.  The criteria file in this case contains sheding percentage, mean and peak log10 viral levels for that individual from the 28-day onservation period.  All criteria files are contained in the data sub-directory.
