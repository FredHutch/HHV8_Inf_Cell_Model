/* The state class holds parameters that have been set either on the command line or 
 * in the input file.  
 */
#ifndef GLOBALSTATE_H
#define GLOBALSTATE_H

#define NUM_PARAMS 18
#define FILE_PARAMS 52
#define MAX_CRIT_CATEGORIES 3

#define VL_BINS 6
#define EPI_BINS 5
#define DUR_BINS 6

#define NUM_CRITERIA (VL_BINS+EPI_BINS+DUR_BINS)

#define MAX_TIME          10
#define MAX_VIRAL_LOAD    10000000000
#define MAX_INF_CELLS     100000
#define MAX_CD8CELLS      1000000

#define MAX_SAMPLES      20000 /* keep 10 days worth of hex measures */

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "hexcell.h"

typedef struct critType {
	    double low;
	    double high;
	    double mean;
	} critType;

class globalState {
    public:
	gsl_rng * ur;

/* set from input file */
	int infect_by_virus;
	double TSim;
	double T0_Sim;
	double T2_Sim;
	double T3_Sim;
	double tstep;
	double Tmindel2;
	double Io;
	double To;
	double Yo;
	double Zo;
	double N0;
	double exp_days_init;
	double exp_days_mean;
	double exp_days_std;
	double beta_init;
	double betae_init;
	double log_p_init;
	double log_p_low;
	double log_p_high;
	double latent_inf_init;
	double theta_init;
	double alpha_init;
	double alpha_mean;
	double alpha_std;
	double kappa_init;
	double kappa_mean;
	double kappa_std;
	double r_init;
	double cd8_ic50_mean;
	double cd8_ic50_init;
	double cd8_ic50_std;
	double delta_init;
	double c_init;
	double inf_init;
	double rinf_init;
	double rho_init;
	double eclipse_init;
	double db;
	double an;
	double fpos;
	int density_killing;
	double hill;
	double hill_mean;
	double hill_std;
	double swabInterval;
	double statInterval;
	double To_mean;
	double To_std;
	double beta_mean;
	double beta_std;
	double betae_mean;
	double betae_std;
	double latent_inf_mean;
	double latent_inf_std;
	double c_mean;
	double c_std;
	double delta_mean;
	double delta_std;
	double fpos_mean;
	double fpos_std;
	double an_mean;
	double an_std;
	double rho_mean;
	double rho_std;
	double theta_mean;
	double theta_std;
	double inf_mean;
	double inf_std;
	double rinf_mean;
	double rinf_std;
	double r_mean;
	double r_std;
	double eclipse_mean;
	double eclipse_std;
	double log_p_mean;
	double log_p_std;

	int PDF_on;
	int gamma_hrs;
	int absorb_hrs;

	/* model 6 & 7 vars */
	double gamma_init;
	double gamma_low;
	double gamma_high;
	double gamma_mean;
	double gamma_std;
	double bolus;
	double absorb_init;
	double absorb_low;
	double absorb_high;
	double absorb_mean;
	double absorb_std;
	double Cmax_init;
	double Cmax_mean;
	double Cmax_std;
	double IC50_init;
	double IC50_mean;
	double IC50_std;
	double m_init;
	double m_low;
	double m_high;
	double m_mean;
	double m_std;

	double Cmax_0;

	/* model 8 vars */
	double kD;
	int cT;
	int nT;

	int N_runs;
	int Total_epis;
	int Tdelay_on;
	int T0_files;
	int use_rho;
	int Tmin;
	int del1;
	int Pulse_regions;
	int Cluster_pulses;
	int yy;
	int Total_doses;
	int infThreshold;

	int Episode_limit;
	int Sig_test;

	double Size_limit;
	double Input_refresh;

	char *inp_file;
	char *alt_inp_file;

	int Calc_T0;

	int Model;
	int Model_0;
	int Model_2;
	int Model_3;

	int Regions;
	int Crit_mask;
	int Match_strategy;
	int Search_order; /* 0 - low to high, 1 - start w/ middle then low then high, 2 - high to low */

/* set from criteria file */
	critType crit[NUM_CRITERIA];
	double critWeight[MAX_CRIT_CATEGORIES];

/* used for graphing */
   	int Max_x;
   	int Max_y;
   	int Max_z;

   	int writeOn;
 	long int max_vl;
 	long int max_inf;
	long int max_cd8s;
	double max_ACV;

	double time;
	double max_time;
	double sampling;
	double hex_time_bias;
	double crit_start;

	int Verbose;
	int CritOn;

	FILE *dataF1;
	FILE *dataF2;
	FILE *dataF3;
	FILE *dataF4;
	FILE *dataF5;
	FILE *dataF6;
	FILE *dataF7;
	FILE *dataF8;
	FILE *dataF9;
	FILE *dataF10;
	FILE *dataF11;
	FILE *dataF12;
	FILE *dataF13;

	long time_st;

 	hexcell **cells;

	int sample_index;
	unsigned long int *vet[MAX_HEXCELLS];
	unsigned long int *vit[MAX_HEXCELLS];
	unsigned long int *inf[MAX_HEXCELLS];
	unsigned long int *cd8[MAX_HEXCELLS];

	double *repro[MAX_HEXCELLS];
	double *diam[MAX_HEXCELLS];
	int *color[MAX_HEXCELLS];

	globalState(void);
	~globalState(void);
};
#endif
