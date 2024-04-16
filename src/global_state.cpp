#include <cstdlib>
using namespace std;

#include <stdio.h>
#include <string.h>

#include "global_state.h"

globalState::globalState(void) 
{
	Max_x = MAX_X;
	Max_y = MAX_Y;
	Max_z = MAX_Z;

	plotVe = 1;
	plotVi = 0;
	plotInf = 0;
	plotCd8 = 0;

	plotStyle1 = 1;
	plotStyle2 = 0;

	plotLogs = 1;
	plotColor = 1;
	plotRegions = 1;
	writeOn = 0;
	scrollAxes = 1;

	plotOpt1 = 1;
	plotOpt2 = 0;
	plotOpt3 = 0;
	plotOpt4 = 0;
	plotOpt5 = 0;
	plotOpt6 = 0;

	max_vl = MAX_VIRAL_LOAD;
	max_inf = MAX_INF_CELLS;
	max_cd8s = MAX_CD8CELLS;
	max_time = MAX_TIME;
	max_ACV = 2.5;

	use_rho = 1;
	PDF_on = 0;
	plot_bias=0;
	hex_time_bias=0;
	plot_span=10;
	sampling=0.01;
	statInterval=365.;
	refresh=1.0;
	gamma_hrs = 0;
	absorb_hrs = 0;

	exp_days_init = 10;
	exp_days_mean = 10;
	exp_days_std = 0;

	density_killing = 0;
	hill = 0.01;
	hill_mean = 0.01;
	hill_std = 0;

	alpha_init = 0.01;
	alpha_high = 0.01;
	alpha_low = 0.01;
	alpha_mean = 0.01;
	alpha_std = 0;

	Max_steps = 1;
	Stop_walk = 1;
	Bvstop_walk = 1;
	Printmax = 1;
	Tolerance = 0.1;
	Threading = 25;
	Fit_model = 1;
	Rand_start = 1;
	Calc_T0 = 1;
	Model = 5;
	Regions = 300;
	Crit_mask = 255;
	Match_strategy = 2;
	Tdelay_on = 0;
	Search_order = 0;
	Pulse_regions = 0;
	Cluster_pulses = 0;
	infThreshold = 10;
	Sig_test = 0;
	writeOn = 0;
	yy = 1.0;
	refresh = 0.5;
	Param_mask = 1;
	Episode_limit = 0;
	Verbose = 0;
	Episode_limit=0;
	Size_limit=0;
	Sig_test=0;
	Pulse_regions=0;
	Cluster_pulses=0;
	Total_doses=0;
	infThreshold=0;

	time_st = 0;

	dataF1 = NULL;
	dataF2 = NULL;
	dataF3 = NULL;
	dataF4 = NULL;
	dataF5 = NULL;
	dataF6 = NULL;
	dataF7 = NULL;
	dataF8 = NULL;
	dataF9 = NULL;
	dataF10 = NULL;
	dataF11 = NULL;
	dataF12 = NULL;
	dataF13 = NULL;

	points = NULL;
	stopFlag = 0;
	pauseFlag = 0;

	Param_mask=0;
	Input_refresh=0;
	inp_file=NULL;
	alt_inp_file=NULL;

	Model_0 = 0;
	Model_2 = 0;
	Model_3 = 0;
	CritOn = 1;
	crit_start = 0;
	AutoSnapshot=0;
	SnapshotInterval=5;

	Total_epis = 0;
	T0_files = 10;
	T2_Sim = 0.;
	T3_Sim = 0.;
	Cmax_0 = 0.0;

	time=0;

	sample_index=0;

	for (int i=0; i < MAX_CRIT_CATEGORIES;i++)
	    critWeight[i] = 1;

	for (int i=0; i < MAX_HEXCELLS;i++)
	{
	    vet[i] = (unsigned long int *)malloc(MAX_SAMPLES * sizeof (unsigned long int));
	    vet[i][sample_index]=0;

	    vit[i] = (unsigned long int *)malloc(MAX_SAMPLES * sizeof (unsigned long int));
	    vit[i][sample_index]=0;

	    inf[i] = (unsigned long int *)malloc(MAX_SAMPLES * sizeof (unsigned long int));
	    inf[i][sample_index]=0;

	    cd8[i] = (unsigned long int *)malloc(MAX_SAMPLES * sizeof (unsigned long int));
	    cd8[i][sample_index]=0;

	    color[i] = (int *)malloc(MAX_SAMPLES * sizeof (int));
	    color[i][sample_index]=0;

	    diam[i] = (double *)malloc(MAX_SAMPLES * sizeof (double));
	    diam[i][sample_index]=0.0;

	    repro[i] = (double *)malloc(MAX_SAMPLES * sizeof (double));
	    repro[i][sample_index]=0.0;
	}
}

globalState::~globalState(void) 
{
	if (points != NULL)
	    delete points;
}
