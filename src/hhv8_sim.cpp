/*** This program is a C++ code for Joshua Schiffer model ***/
/*** basically it simulates the stochastic pattern of Hsv shedding ***/ 	 
/*** for an individual over a period of time this is called a run ***/
/*** multiple runs are generated for a given parameter set and the ***/
/*** results are summarized for that set in five summary measures ***/
/*** a grid of parameter values is set up at the beginning of the ***/
/*** code and thus summary measures are obtained for the space of ***/
/*** parameter values ***********************************************/

/************ Ramzi Alsallaq in collaboration with ****************** 
 ************ Joshua Schiffer Amalia Magaret  January 2009 *********/
/************ Extended to include model5 and matching crtiteria ***** 
 ************ Dave Swan and Joshua Schiffer May/June 2010 *********/
/************ Extended to include model6 and matching crtiteria ***** 
 ************ Dave Swan and Joshua Schiffer January 2011 *********/


/** Program expects criteria in file hhv8_sim.crit and parameters from stdin (i.e. use < operator)

**/

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include <map>
#include <sys/types.h>
#include <unistd.h>
#include <stdarg.h>

#ifdef __sun
#include <strings.h>
#else
#include <string.h>
#endif

// Some STL Headers
#include <vector>
#include <stdlib.h>

// Using The STL Exception Library Increases The
// Chances That Someone Else Using Our Code Will Correctly
// Catch Any Exceptions That We Throw.
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include "state_var.h"  
#include "global_state.h"  

#include "hexcell.h"

using namespace std;

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#ifdef __sun
#define DEFAULT_TITLE  "HHV8 simulation (for SunOS)"
#else
#define DEFAULT_TITLE  "HHV8 simulation"
#endif
#define PROGRAM_NAME  "hhv8_sim"
#define PACKAGE_VERSION  "1.0"
#define PACKAGE_BUGREPORT  "dswan@fredhutch.org"

#define MAX_SWABS 10000

static globalState *theState=NULL;
static int stoch;

static string outDir;

void abort_(const char * s, ...)
{
        va_list args;
        va_start(args, s);
        vfprintf(stderr, s, args);
        fprintf(stderr, "\n");
        va_end(args);
        abort();
}
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define MAX_LINE 80

#define SET_PARAMETER(name,param) \
    if (inputs.find(name) != inputs.end()) { vars->param = inputs[name]; } else if (!refresh) { \
	cerr << "Missing setting for double precision parameter "<<name<<".  Exiting!\n"; \
	exit(1); \
    }

#define SET_OPT_PARAMETER(name,param) \
    if (inputs.find(name) != inputs.end()) { vars->param = inputs[name]; } else if (!refresh) { \
	cerr << "No setting for optional parameter "<<name<<".\n"; \
    }

#define SET_INT_PARAMETER(name,param) \
    if (inputs.find(name) != inputs.end()) { vars->param = static_cast<int>(inputs[name]); } else if (!refresh) { \
	cerr << "Missing setting for integer parameter "<<name<<".  Exiting!\n"; \
	exit(1); \
    }

#define SET_OPT_INT_PARAMETER(name,param) \
    if (inputs.find(name) != inputs.end()) { vars->param = static_cast<int>(inputs[name]); } else if (!refresh) { \
	cerr << "No setting for optional parameter "<<name<<".\n"; \
    }

//
//comparator necessary for qsort
int compare_doubles (const void * a, const void * b)
{
  if (*(double *)a > *(double *)b)
     return 1;
  if (*(double *)a < *(double *)b)
     return -1;
  return 0;
}


//max fuction on sequence
double getMax(double *arr,int size)
{
    double the_max=0.;

    for (int i=0; i < size; i++)
	if (*(arr+i) > the_max)
	    the_max = *(arr+i);

    return the_max;
}

//mean fuction on sequence
double getMean(double *arr,int size)
{
    double mean;
    double total=0.;

    for (int i=0; i < size; i++)
	total += *(arr+i);

    mean = total/size;

    return mean;
}

//stddev fuction on sequence
double getStddev(double *arr,int size)
{
    double mean;
    double stddev;
    double total=0.;
    double err=0.;

    for (int i=0; i < size; i++)
	total += *(arr+i);

    mean = total/size;

    for (int i=0; i < size; i++)
	err+=pow((mean-*(arr+i)),2.0);

    stddev=sqrt(err/(size));
    return stddev;
}
//stddev fuction on sequence
double getStateStddev(state_var *svar[],int size)
{
    double mean;
    double stddev;
    double total=0.;
    double err=0.;

    for (int i=0; i < size; i++)
	total += svar[i]->val();

    mean = total/size;

    for (int i=0; i < size; i++)
	err+=pow((mean-svar[i]->val()),2.0);

    stddev=sqrt(err/(size));
    return stddev;
}

//median fuction after sorting sequence
double getMedian(double *sorted_arr,int size)
{
    int middle = size/2;
    double median;
    if (size%2==0) 
	median = (*(sorted_arr+middle-1)+*(sorted_arr+middle))/2.;
    else 
	median = *(sorted_arr+middle);

    return median;
    //cout << "Median of sorted_array is: " << average << endl;
}

double calcDecline (double time, double tstep, double delta, int region, double T, gsl_matrix *pastTps, globalState *vars) 
{
	double delTp;

	int slot, next_slot;

	if (vars->Tdelay_on == 0)
	    return delta;

	if (time < vars->del1)
	    slot = 0;
	else
	    slot = (int)(time - ((int)(time/vars->del1)*vars->del1));

	/* read out past value */
        delTp = gsl_matrix_get(pastTps,region,slot);

	/* put in current for del1 days from now 
           but only once per day (just before slot changes)!*/
	next_slot = (int)(time+tstep - ((int)((time+tstep)/vars->del1)*vars->del1));
	if (next_slot != slot && time > vars->del1)
	{
	    /*fprintf(stderr,"oldT=%lf,newT=%lf,time=%lf\n",delTp,T,time);*/
	    gsl_matrix_set(pastTps,region,slot,T);
	}

	if (delTp <= T)
	    return delta;
	
	if (T < vars->Tmindel2)
	    return delta/10.;

	return delta/2.0;
}

double pulse(double volume, double start_time, double freq, double time, double deltat, double *next_pulse)
{
    if (*next_pulse == 0.)
	*next_pulse = start_time;	

    if (time > *next_pulse)
    {
	*next_pulse += freq;
	return volume;
    }

    return 0.0;
}

int pick_neighbor(int cellIndex, globalState *vars)
{
    int picked = (int)gsl_rng_uniform_int(vars->ur,MAX_NEIGHBORS);
    while (vars->cells[cellIndex]->neighbors[picked] == NULL)
	picked = (int)gsl_rng_uniform_int(vars->ur,MAX_NEIGHBORS);
    return vars->cells[cellIndex]->neighbors[picked]->cell;
}

/* 
 * model5 is actually called for all spatial models.  These include
 * 	model 4 - no neighborhood effects
 * 	model 5 - localized spread
 * 	model 6 - ACV 
 * 	model 7 - ACV & 
 */
bool model5 (double time, state_var* Sp, state_var* Ssub[], 
	state_var* Ec[], state_var* Ip[], 
	state_var* Tp[], gsl_matrix *pastTps, 
	state_var* Vi[], state_var*  Ve[], 
	state_var*  Veadj[], unsigned int *Iplaquebirths, 
	int *totalPlaques, int plaqColors[], 
	globalState *vars, 
	double Repro[], double logRepro[], double p, double latent_inf, 
	state_var* Vethis[], state_var* Vithis[], state_var* Ithis[],state_var* Id_this[], double start_inf_time[],
	double *repro_max,
	double *repro_mean,
	double *log_repro_mean,
	double *repro_std,
	double *log_repro_std,
	double *newTinf,
	double *Tbump,
	double *Tdump,
	double *Tkills,
	double *Ideaths,
	unsigned long int *newTcells, 
	unsigned long int *newInfcells, 
	unsigned long int *newInfevents, 
	unsigned long int *newVirons, 
	unsigned long int *IthisTot, unsigned long int *IdthisTot, unsigned long int *VithisTot, 
	unsigned long int *VethisTot, int selectedRegions[]) 
{
        double tstep;   
	double db, beta, beta_e, hill, fpos,  theta, FIp, FInf, delta, c, 
		rho, r, an, inf, rinf, eclipse, cd8_ic50, alpha, kappa, expand_days; 
	double N0, decline;

	unsigned long int Inf_e[MAX_HEXCELLS], Inf_neu;

	long int Vedelta;
	unsigned long int Vesaved;

	long int Videlta;
	unsigned long int Visaved;

	unsigned long int Iptot, Vetot, Vitot;

	long int Birth, Death, Inf, TDeath_Death, TDeath, Tinf,Exhaust;

	/* set locals from global stat vars */

	tstep = vars->tstep; 
	db = vars->db; 
	an = vars->an; 
	fpos = vars->fpos; 
	hill = vars->hill; 
	cd8_ic50 = vars->cd8_ic50_init; 
	kappa = vars->kappa_init; 
	alpha = vars->alpha_init; 
	expand_days = vars->exp_days_init; 
	N0 = vars->N0; 

	/* set parameters using passed param vector */

	/* NOTE: latent_inf and p passed in to allow for drug effects (in models 6 & 7) */

	beta=vars->beta_init;
	p = pow(10,vars->log_p_init);
	c=vars->c_init;
	theta=vars->theta_init;
	inf=vars->inf_init;
	r=vars->r_init;
	rinf=vars->rinf_init;
	delta=vars->delta_init;
	beta_e=vars->betae_init;
	rho=vars->rho_init;
	eclipse=vars->eclipse_init;

////////////////////////////////////////////////////////////////////////
// update time-dependent parameters

	Vetot=0;
	Vitot=0;
	unsigned long int Velast=0;
	unsigned long int Vilast=0;

	int highest=0;
	unsigned long int highCnt=0;
	int highest2=0;
	unsigned long int highCnt2=0;


	for (int i =0; i < vars->Regions; i++)
	{
	    Velast=Vetot;
	    Vetot += Ve[i]->val();
	    /* Watch for roll-over! */
	    if (Vetot < Velast)
	    {
		fprintf(stderr,"Error: Vetot rolled over!\n");
		fprintf(stderr,"time=%g,i=%d, Vetot=%lu, Ve[i]=%lu\n, Velast=%lu",
			time,i,Vetot,Ve[i]->val(),Velast);
		fprintf(stderr,
		    "beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
		    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse);
		exit(1);
	    }
	    if (theState->Verbose > 1 && Vetot > 4e9 && Velast < 4e9)
	    {
		fprintf(stderr,"Vetot Over 4 billion (%lu) at t=%lf! (was %lu)\n",Vetot, time, Velast);
		fprintf(stderr,
		    "beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
		    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse);
	    }

	    Vilast=Vitot;
	    Vitot += Vi[i]->val();

	    /* Watch for Vi roll-over! */
	    if (Vitot < Vilast)
	    {
		fprintf(stderr,"Error: Vitot rolled over!\n");
		fprintf(stderr,"time=%g,i=%d, Vitot=%lu, Vi[i]=%lu\n, Vilast=%lu",
			time,i,Vitot,Vi[i]->val(),Vilast);
		fprintf(stderr,
		    "beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
		    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse);
		exit(1);
	    }
	    if (Vi[i]->val() > highCnt)
	    {
		highest2=highest;
		highCnt2=highCnt;
		highest=i;
		highCnt=Vi[i]->val();
	    }
	    else if (Vi[i]->val() > highCnt2)
	    {
		highest2=i;
		highCnt2=Vi[i]->val();
	    }
	    if (theState->Verbose > 1 && Vitot > 4e9 && Vilast < 4e9)
	    {
		fprintf(stderr,"Vitot Over 4 billion (%lu) at t=%lf! (was %lu)\n",Vitot, time, Vilast);
		fprintf(stderr,"\tHighest cell is %d (%lu)\n",highest,highCnt);
		fprintf(stderr,"\t2nd Highest cell is %d (%lu)\n",highest2,highCnt2);
		fprintf(stderr,
		    "beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
		    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse);
	    }
	}

	// update system space

	// HHV8+ partner

	// Susceptible cells
	long int Sinf[MAX_HEXCELLS];
	long int Ecl_tot = 0, Sbirthstot=0, Sdeathtot=0, Sinftot=0;
	long int Sripe[MAX_HEXCELLS];

	*newTinf=0;
	*Tbump=0;
	*Tdump=0;
	*Tkills=0;
	*Ideaths=0;
	*newTcells=0;
	*newInfcells=0; 
	*newInfevents=0; 
	*newVirons=0; 

	for (int i =0; i < vars->Regions; i++)
	{
	    Birth = (stoch == 1) ? gsl_ran_poisson (vars->ur, db*N0*tstep/vars->Regions): 
		(db*N0*tstep/vars->Regions); 

	    Sbirthstot += Birth;


	    Death = (stoch == 1) ? gsl_ran_poisson (vars->ur, db*tstep*Ssub[i]->val()):
		(db*tstep*Ssub[i]->val());

	    Sdeathtot += Death;

	    if (vars->infect_by_virus==1) {
		//Sinf[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, beta*Vi[i]->val()*tstep, Sp->val() / vars->Regions): (beta*Vi[i]->val()*tstep* Sp->val() / vars->Regions);
		Sinf[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, beta*tstep, Vi[i]->val()): (beta*Vi[i]->val()*tstep);
		//Inf_e[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, beta_e*Ve[i]->val()*tstep, Sp->val() / vars->Regions): (beta_e*Ve[i]->val()*tstep* Sp->val() / vars->Regions);
		Inf_e[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, beta_e*tstep, Ve[i]->val()): (beta_e*Ve[i]->val()*tstep);
	    } else {
		Sinf[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, beta*tstep, Ip[i]->val()): (beta*tstep* Ip[i]->val());
		Inf_e[i] = 0;
	    }
	    Sinftot += Sinf[i] + Inf_e[i];

	    Ssub[i]->update_with (Birth - Death- Sinf[i] - Inf_e[i]);
        }

	Inf_neu = (stoch == 1) ? gsl_ran_poisson (vars->ur, latent_inf*tstep):latent_inf*tstep;

	Sp->update_with (Sbirthstot - Sdeathtot - Inf_neu -Sinftot);

	*newInfevents = Sinftot;

	/* limit Sp to N0 */
	if (Sp->val() > N0)
	    Sp->initialz_with(N0);

	// New plaques from neuronal virus
	//*Iplaquebirths = (*Iplaquebirths > 0) ? 0 : Inf_neu;
	*Iplaquebirths = Inf_neu;

	// ******************************NEW CODE (11/2) *****************************************
	//
	// Plaque determinations
	unsigned int Inewplaques[MAX_HEXCELLS];

	// set all plaques to zero
	for (int i =0; i < vars->Regions; i++)
	    Inewplaques[i]=Sinf[i];

	// for any infected cell (Inf_e > 0), pick a neighbor and start a plaque there
	for (int i =0; i < vars->Regions; i++)
	    if (Inf_e[i] > 0)
	    {
		int index;

		/* when using "neighbor effects" plaques only spread to neighbors */
		if (vars->Regions > 1 && vars->Model >= 5)
		    index = pick_neighbor(i,vars);
		else if (vars->Regions > 1)
		    index = (int)gsl_rng_uniform_int(vars->ur,vars->Regions);
		else
		    index = 0;

		Inewplaques[index] += Inf_e[i];

		if (theState->Verbose > 1)
		    fprintf(stderr,"Spread plaque (color=%d,to=%d,from=%d) at t=%g\n",
			plaqColors[index],index,i,time);

		if (Ec[index]->val()==0 && Ip[index]->val()==0)
		    (*totalPlaques)++;
	    }

	int plaqueRegion;
	if (vars->Regions > 1 && vars->Pulse_regions == 0)
	{
	    plaqueRegion = (int)gsl_rng_uniform_int(vars->ur,vars->Regions);
	    /* when using "neighbor effects" plaques shouldn't start at the edges */
	    if (vars->Model >= 5)
	    {
		while (vars->cells[plaqueRegion]->num_neighbors != 6)
		    plaqueRegion = (int)gsl_rng_uniform_int(vars->ur,vars->Regions);
	    }
	}
	else if (vars->Regions > 1)
	{
	    plaqueRegion = (int)gsl_rng_uniform_int(vars->ur,vars->Pulse_regions);
	    if (!vars->Cluster_pulses)
		plaqueRegion=selectedRegions[plaqueRegion];
	}
	else
	    plaqueRegion = 0;

	for (int i =0; i < vars->Regions; i++)
	{
	    if (i == plaqueRegion && *Iplaquebirths > 0)
	    {
	    	Inewplaques[i]+=*Iplaquebirths;
		if (theState->Verbose > 1)
		    fprintf(stderr,
			"Created plaque (color=%d,reg=%d) at t=%g\n",plaqColors[i],i,time);

		if (Ec[i]->val()==0 && Ip[i]->val()==0)
		    (*totalPlaques)++;
	    }	

	    if (eclipse > 0.001)
	    {
		Sripe[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, (1.0/eclipse)*tstep,Ec[i]->val()): ((1.0/eclipse)*tstep*Ec[i]->val());
	    }
	    else
		Sripe[i] = Inewplaques[i];

	    Ec[i]->update_with (Inewplaques[i] - Sripe[i]);	
	}

	Ecl_tot=0;
	for (int i =0; i < vars->Regions; i++)
	    Ecl_tot += Ec[i]->val();

	for (int i =0; i < vars->Regions; i++)
	{
	    Birth = Sripe[i];
	    *IthisTot += Birth;
	    Ithis[i]->update_with (Birth);
	    if (Birth>0)
		*newInfcells += Birth;

	    if (vars->density_killing) 
	    {
		TDeath_Death = (stoch == 1) ? gsl_ran_binomial (vars->ur, (an+pow(Ip[i]->val(),(1+hill))+fpos*Tp[i]->val())*tstep, Ip[i]->val()): ((an*(1+pow(Ip[i]->val(),hill))+fpos*Tp[i]->val())*tstep*Ip[i]->val());

		TDeath = (stoch == 1) ? gsl_ran_binomial (vars->ur, fpos*Tp[i]->val()/(an+pow(Ip[i]->val(),(1+hill))+fpos*Tp[i]->val()), TDeath_Death): (fpos*Tp[i]->val()/(an+pow(Ip[i]->val(),(1+hill))+fpos*Tp[i]->val())* TDeath_Death);

	    } else {
		TDeath_Death = (stoch == 1) ? gsl_ran_binomial (vars->ur, (an+fpos*Tp[i]->val())*tstep, Ip[i]->val()): ((an+fpos*Tp[i]->val())*tstep* Ip[i]->val());

		TDeath = (stoch == 1) ? gsl_ran_binomial (vars->ur, fpos*Tp[i]->val()/(an+fpos*Tp[i]->val()), TDeath_Death): (fpos*Tp[i]->val()/(an+fpos*Tp[i]->val())* TDeath_Death);
	    }
	    Death = TDeath_Death - TDeath;

	    Id_this[i]->update_with (TDeath_Death);
	    *IdthisTot += TDeath_Death;
	    *Tkills += TDeath;
	    *Ideaths += Death;

	    //Ip[i]->update_with (Inewplaques[i] + Birth - TDeath_Death);	
	    Ip[i]->update_with (Birth - TDeath_Death);	
	}

	Iptot=0;
	for (int i =0; i < vars->Regions; i++)
	    Iptot += Ip[i]->val();

	// Virus, intercellular (Vi)
	//
	long int Videath[MAX_HEXCELLS];
	long int Vitcelldeath[MAX_HEXCELLS];

	for (int i =0; i < vars->Regions; i++)
	{
	    // protect against poisson hangups!
	    if (Ip[i]->val()*p*tstep < 4e9)
	    {
		Birth = (stoch == 1) ? gsl_ran_poisson (vars->ur, p*Ip[i]->val()*tstep)
		    : (p*Ip[i]->val()*tstep);
	    }
	    else
		Birth = (p*Ip[i]->val()*tstep);

	    *VithisTot += Birth;
	    Vithis[i]->update_with (Birth);

	    // protect against binomial blowups!
	    if (Vi[i]->val() < 4e11)
	    {
		Inf = (stoch == 1) ? 
		    gsl_ran_binomial (vars->ur, beta*tstep, Vi[i]->val()): 
			(beta*tstep* Vi[i]->val());

		if (vars->density_killing) 
		{
		    Videath[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, an*pow(Ip[i]->val(),(1+hill))*tstep,Vi[i]->val()-Inf): (an*pow(Ip[i]->val(),(1+hill))*tstep*Vi[i]->val()-Inf);
		} else {
		    Videath[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, an*tstep,Vi[i]->val()-Inf): (an*tstep*Vi[i]->val()-Inf);
		}

		Vitcelldeath[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, tstep*fpos*Tp[i]->val(), Vi[i]->val()-Videath[i]): (tstep*fpos*Tp[i]->val()* (Vi[i]->val()-Videath[i]));

	    }
	    else
	    {
		Inf = (beta*tstep* Vi[i]->val());

		if (vars->density_killing) 
		{
		    Videath[i] = (an*pow(Ip[i]->val(),(1+hill))*tstep*Vi[i]->val()-Inf);
		} else {
		    Videath[i] = (an*tstep*Vi[i]->val()-Inf);
		}

		Vitcelldeath[i] = (tstep*fpos*Tp[i]->val()* (Vi[i]->val()-Videath[i]));
	    }


	    Videlta = (Birth - Videath[i] - Vitcelldeath[i] - Inf);

	    Visaved=Vi[i]->val();

	    if (Videlta < 0 && Vi[i]->val() < -Videlta)
		Videlta=-Vi[i]->val();

	    Vi[i]->update_with (Videlta);

	    if (Videlta < 0 && Vi[i]->val() > Visaved)
	    {
		fprintf(stderr,"Error: Vi declined by more than exists!\n");
		fprintf(stderr,"time=%g,i=%d, Vi[i]=%lu (was %lu) , Videlta=%ld\n",time,i,Vi[i]->val(),Visaved,Videlta);
		exit(1);
	    }

	    /* watch for roll-over! */
	    if (Videlta > 0 && Vi[i]->val() < Visaved)
	    {
		fprintf(stderr,"Error: Vi got too big!\n");
		fprintf(stderr,"time=%g,i=%d, Vi[i]=%lu, Videlta=%ld\n, Visaved=%lu",
			time,i,Vi[i]->val(),Videlta,Visaved);
		exit(1);
	    }
	    // invalidate score if model produced more than 4 billion virons!
	    if (Vi[i]->val() > 4e11)
	    {
		/*if (vars->Verbose)*/
		    fprintf(stdout,"Error: Vi[%d] got big at t=%lf! (%lu).  bigger than 4e11.\n",
			i,time,Vi[i]->val());
		/*if (vars->Verbose)*/
		    fprintf(stdout,
			"beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf\n",
			beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho);
		    exit(1);
#ifdef INVALIDATE
		Vi[i]->initialz_with(4e11);
#endif
	    }
	}
	//
	// Virus, extracellular (Ve and Veadj)
	//
	unsigned int AdjDeath;

	for (int i =0; i < vars->Regions; i++)
	{
	    Birth = Videath[i] + Vitcelldeath[i]; // Vi from inf cell death or killing
	    *VethisTot += Birth;
	    Vethis[i]->update_with (Birth);

	    if (Birth>0)
		*newVirons += Birth;

	    if (Ve[i]->val() < 4e11)
	    {
		Death = (stoch == 1) ? 
		    gsl_ran_poisson (vars->ur, c*tstep*Ve[i]->val()): 
			(c*tstep*Ve[i]->val());

	    }
	    else
	    {
		Death = (c*tstep*Ve[i]->val());

	    }

	    Vedelta = (Birth - Death);

	    Vesaved=Ve[i]->val();

	    if (Vedelta < 0 && Ve[i]->val() < -Vedelta)
		Vedelta=-Ve[i]->val();

	    Ve[i]->update_with (Vedelta);

	    //limit Ve[i] to 4 billion!
	    if (Vedelta < 0 && Ve[i]->val() > Vesaved)
	    {
		fprintf(stderr,"Error: Ve declined by more than exists!\n");
		fprintf(stderr,"time=%g,i=%d, Ve[i]=%lu (was %lu) , Vedelta=%ld (B=%ld,D=%ld,I=%ld)\n",
		    time,i,Ve[i]->val(),Vesaved,
		    Vedelta,Birth,Death,Inf);
		exit(1);
	    }

	    /* watch for roll-over! */
	    if (Vedelta > 0 && Ve[i]->val() < Vesaved)
	    {
		fprintf(stderr,"Error: Ve got too big!\n");
		fprintf(stderr,"time=%g,i=%d, Ve[i]=%lu (was %lu) , Vedelta=%ld (B=%ld,D=%ld,I=%ld)\n",
		    time,i,Ve[i]->val(),Vesaved,
		    Vedelta,Birth,Death,Inf);
		exit(1);
	    }
	    // invalidate score if model produced more than 4 billion virons!
	    if (Ve[i]->val() > 4e11)
	    {
		/*if (vars->Verbose)*/
		fprintf(stdout,"Warning: Ve[%d] bigger than 4e11.\n", i);
		fprintf(stdout,"time=%g,i=%d, Ve[i]=%lu (was %lu) , Vedelta=%ld (B1=%ld,B2=%ld,D=%ld,I=%ld)\n",
		    time,i,Ve[i]->val(),Vesaved,
		    Vedelta,Videath[i],Vitcelldeath[i],Death,Inf);
	    
		/*if (vars->Verbose)*/
		fprintf(stdout,
		    "beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf\n",
		    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho);
		exit(1);
#ifdef INVALIDATE
		Ve[i]->initialz_with(4e11);
#endif
	    }
	    
	    AdjDeath = (stoch == 1) ? gsl_ran_binomial (vars->ur, c*tstep,Veadj[i]->val()): (c*tstep*Veadj[i]->val());

	    Veadj[i]->update_with (Birth - AdjDeath);
	}

	// T-cells
	//


	for (int i =0; i < vars->Regions; i++)
	{
	    // *****************NEW CODE (11/2) **********************
	    //
	    double Tmean=Tp[i]->val();
	    int cnt=1;
	    for (int j =0; j < MAX_NEIGHBORS; j++)
		if (vars->cells[i]->neighbors[j] != NULL)
		{
		    Tmean += Tp[vars->cells[i]->neighbors[j]->cell]->val();
		    cnt++;
		}

	    Tmean = Tmean / cnt;
	    // Only allow expansion for expand_days days!
	    if (expand_days == 0 || start_inf_time[i] == 0 || time - start_inf_time[i] < expand_days)
		FIp = (double)Ip[i]->val()/(double)(Ip[i]->val()+r);
	    else
		FIp = 0;

	    if (Ip[i]->val() > 0)
		//FInf = Iptot/(Iptot+rinf); //old way
		FInf = Ip[i]->val()/(Ip[i]->val()+rinf); // revised 10/25
	    else
		FInf = 0.;

	    decline = calcDecline (time, tstep, delta, i, Tp[i]->val(), pastTps, vars);

	    Tinf = (stoch == 1) ? gsl_ran_poisson (vars->ur,inf*FInf*tstep)
	     : (inf*FInf*tstep);
	    *newTinf += Tinf;

	    Birth = (stoch == 1) ? gsl_ran_poisson (vars->ur,(alpha+theta*FIp*Tp[i]->val())*tstep)
	    : ((alpha+theta*FIp*Tp[i]->val())*tstep);

	    if (Birth>0)
		*newTcells += Birth;

	    Death = (stoch == 1) ? gsl_ran_poisson (vars->ur,decline*Tp[i]->val()*tstep)
	    : (decline*Tp[i]->val()*tstep);

	    double FExh=(double)Tp[i]->val()/((double)Tp[i]->val()+cd8_ic50);
	    Exhaust = (stoch == 1) ? gsl_ran_poisson (vars->ur,kappa*FExh*tstep*Tp[i]->val())
	    : (kappa*FExh*tstep*Tp[i]->val());

	    if (vars->use_rho && Inewplaques[i] > 0)
	    {
		long unsigned int prevTcells = Tp[i]->val();
		Tp[i]->initialz_with(Tp[i]->val()*rho + 
		((stoch == 1)?gsl_ran_poisson(vars->ur,Tmean): Tmean)*(1.0-rho));
		if (Tp[i]->val() > prevTcells)
		   *Tbump+=(double)(Tp[i]->val())-(double)(prevTcells);
		else
		   *Tdump+=(double)(prevTcells)-(double)(Tp[i]->val());
	    }
	    else
	    {
		Tp[i]->update_with (Tinf + Birth - Death - Exhaust);
		//Tp[i]->update_with (Birth - Death);
	    }
	    if (alpha == 0 && Tp[i]->val() < vars->Tmin)
		Tp[i]->initialz_with(vars->Tmin);

	    if (Tp[i]->val() > 1e7)
	    {
		fprintf(stderr,"Error: Tp[%d] bigger than 1e7. Exiting\n", i);
		fprintf(stderr,"time=%g,i=%d, Tp[i]=%lu Ip=%lu (B=%ld,D=%ld,I=%ld,E=%ld)\n",
		    time,i,Tp[i]->val(),Ip[i]->val(),Birth,Death,Tinf,Exhaust);
	    
		/*if (vars->Verbose)*/
		    fprintf(stderr,
			"beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf kappa=%lf alpha=%lf cd8_ic50=%lf fpos=%lf\n",
			beta,latent_inf,p,c,theta,delta,r,inf,rinf,kappa,alpha,cd8_ic50,fpos);
		exit(1);
	    }
		
	}


	// reproduction number
	for (int i =0; i < vars->Regions; i++)
	{
	    Repro[i] = (beta / vars->Regions) * p / pow((an + (Tp[i]->val() * fpos)),2.0);

	    if (plaqColors[i] > 0 && Ve[i]->val() < 1 && Vi[i]->val() < 1 && Ec[i]->val() < 1)
	    {
		if (theState->Verbose > 1)
		    fprintf(stderr,"Cleared plaque in %d (old color=%d) at t=%g\n",
			i, plaqColors[i],time);
		plaqColors[i] = 0;
	    }
	}
	*repro_max = getMax(Repro,vars->Regions);
	*repro_mean = getMean(Repro,vars->Regions);
	*repro_std = getStddev(Repro,vars->Regions);
	*log_repro_mean = getMean(logRepro,vars->Regions);
	*log_repro_std = getStddev(logRepro,vars->Regions);
	return true;
}

inline double shed_cpy (int n)
{return pow(10.,n);} 
inline unsigned int shedbin (unsigned int VL)
{
	int k;
	if(VL<=100) 
		return 1;//[0,100] interval 
	else
	{	
		for (k=2; k<=10; ++k)
		{
			if(double(VL) > shed_cpy(k) && double(VL) <= shed_cpy(k+1)) return k;

		}
	}		

}

/* Function to compute Siliciano's cruitical subset */
double siliciano(int n, int c, double k, double D)
{
    double numer=0, denom=0;

    /* sum (nt i)(D/k)**i */
    for (int i=(n-c+1); i<=n; i++)
    {
	numer += gsl_sf_choose(n,i) * pow((D/k),i);
    }

    /* sum (nt i)(D/k)**i */
    for (int i=0; i<=n; i++)
    {
	denom += gsl_sf_choose(n,i) * pow((D/k),i);
    }
    return numer/denom;
}
/* Function to compute adjustment to drug concentration based on bolus schedule and decay param */
double deltaACV(double time, double tstep, double Cmax, double absorb, 
	double ACV, double gamma, double bolus, double *infuse, 
	int *firstBolus, double *lastBolus, int *doses, int total_doses)
{
    double decay;
    
    decay = gamma*ACV*tstep;	

    /* reset lastBolus if it is time for a boost */
    if (time - *lastBolus >= bolus)
    {
	*lastBolus = time;
	(*doses)++;
	*infuse = (Cmax/absorb) * tstep;
    }

    /* if not yet absorbed, calculate infusion */
    if (time - *lastBolus < absorb && 
	(total_doses==0 || *doses <= total_doses))
	*infuse = (Cmax/absorb) * tstep;
    else
    {
	*infuse = 0.;
	*firstBolus=0;
    }

    return *infuse - decay;
}

void read_input_file(int refresh, char *inp_file, globalState *vars)
{
    ////////////////////////////////////////////////////////////////////////
    ///// read input parameters through the input file 
    char tmpline[MAX_LINE];
    char *valuep;
    int i=0;
    FILE *inf;

    map<string,double> inputs;
    string par;
    double parv;

    if ((inf = fopen (inp_file,"r")) == NULL) {
	cerr << "Could not open input file "<<inp_file<<"\n";
	exit(1);
    }
    while (fgets(tmpline, MAX_LINE-1,inf) != NULL) {
	i++;
	tmpline[MAX_LINE-1] = '\0';
	valuep = rindex(tmpline,' ');
	if (valuep == NULL) {
	    cerr << "Error while reading parameter name from "<<inp_file<<" at parameter #"<<i<<"\n";
	    exit(1);
	}
	*valuep = '\0';
	par = tmpline;
	
	if (sscanf(valuep+1,"%lf",&parv) != 1) {
	    cerr << "Error while reading value for "<<par<<" in "<<inp_file<<" (parameter #"<<i<<")\n";
	    exit(1);
	}
	inputs[par] = parv;
	cout <<"#"<< i <<" "<<par <<" = "<< inputs[par];
	cout<<endl;
    }
    cout<<endl;
    cout <<"Finished reading input file "<<inp_file<<".\n";
    fclose (inf);

#ifdef STREAMS
    ifstream inpF;
    inpF.exceptions(fstream::eofbit |ifstream::failbit | ifstream::badbit );

    try {
	inpF.open(inp_file);
    } catch (ifstream::failure e)
    {
	cerr << "Could not open input file "<<inp_file<<"\n";
	exit(1);
    }

    try {
	int i=0;

	while(inpF.peek() != EOF){
	    inpF >> par;
	    inpF >> parv;
	    inputs[par] = parv;
	    i++;
	    cout <<"#"<< i <<" "<<par <<" = "<< inputs[par];
	    cout<<endl;
	}	
	cout<<endl;
	cout <<"Finished reading input file "<<inp_file<<".\n";
	inpF.close();
    } catch (ifstream::failure e)
    {
	cerr << "Error while reading input file "<<inp_file<<"\n";
	exit(1);
    }
#endif

    ////////////////////////////////////////////////////////////////////////
    ///// read parameters through piping a file that has them //////////////

    //assign parameters:

    SET_OPT_INT_PARAMETER("infect_by_virus",infect_by_virus);

    SET_PARAMETER("TSim",TSim);
    SET_PARAMETER("T0_Sim",T0_Sim);
    SET_OPT_PARAMETER("T2_Sim",T2_Sim);
    SET_OPT_PARAMETER("T3_Sim",T3_Sim);
    SET_PARAMETER("tstep",tstep);
    SET_PARAMETER("Tmindel2",Tmindel2);

    SET_PARAMETER("Io",Io);
    SET_PARAMETER("To",To);
    SET_OPT_PARAMETER("N0",N0);

    SET_OPT_PARAMETER("alpha_init",alpha_init);
    SET_OPT_PARAMETER("alpha_mean",alpha_mean);
    SET_OPT_PARAMETER("alpha_std",alpha_std);

    SET_PARAMETER("beta_init",beta_init);
    SET_PARAMETER("log_p_init",log_p_init);
    SET_PARAMETER("log_p_low",log_p_low);
    SET_PARAMETER("log_p_high",log_p_high);
    SET_PARAMETER("latent_inf_init",latent_inf_init);
    SET_PARAMETER("theta_init",theta_init);
    SET_PARAMETER("r_init",r_init);
    SET_PARAMETER("delta_init",delta_init);
    SET_PARAMETER("c_init",c_init);
    SET_PARAMETER("inf_init",inf_init);
    SET_PARAMETER("inf_mean",inf_mean);
    SET_PARAMETER("inf_std",inf_std);
    SET_PARAMETER("rinf_init",rinf_init);
    SET_PARAMETER("rinf_mean",rinf_mean);
    SET_PARAMETER("rinf_std",rinf_std);
    SET_PARAMETER("rho_init",rho_init);
    SET_PARAMETER("db",db);
    SET_PARAMETER("To_mean",To_mean);
    SET_PARAMETER("To_std",To_std);
    SET_PARAMETER("an",an);
    SET_PARAMETER("fpos",fpos);
    SET_OPT_INT_PARAMETER("density_killing",density_killing);
    SET_OPT_PARAMETER("hill",hill);
    SET_PARAMETER("swabInterval",swabInterval);
    SET_OPT_PARAMETER("statInterval",statInterval);

    SET_OPT_INT_PARAMETER("del1",del1);
    SET_OPT_INT_PARAMETER("Tdelay_on",Tdelay_on);
    SET_INT_PARAMETER("Tmin",Tmin);
    SET_INT_PARAMETER("N_runs",N_runs);

    SET_INT_PARAMETER("Calc_T0",Calc_T0);
    SET_INT_PARAMETER("writeOn",writeOn);
    SET_INT_PARAMETER("Regions",Regions);
    SET_INT_PARAMETER("Crit_mask",Crit_mask);
    SET_INT_PARAMETER("Match_strategy",Match_strategy);
    SET_OPT_INT_PARAMETER("Episode_limit",Episode_limit);
    SET_OPT_PARAMETER("Size_limit",Size_limit);
    SET_OPT_PARAMETER("crit_start",crit_start);

    SET_INT_PARAMETER("Model",Model);
    SET_OPT_INT_PARAMETER("Model_2",Model_2);
    SET_OPT_INT_PARAMETER("Model_3",Model_3);
    SET_OPT_INT_PARAMETER("Verbose",Verbose);
    SET_OPT_INT_PARAMETER("Total_epis",Total_epis);
    SET_OPT_INT_PARAMETER("T0_files",T0_files);
    SET_PARAMETER("Sampling",sampling);
    SET_OPT_PARAMETER("Input_refresh",Input_refresh);
    SET_OPT_INT_PARAMETER("Total_doses",Total_doses);
    SET_OPT_INT_PARAMETER("CritOn",CritOn);
    SET_OPT_INT_PARAMETER("PDF_on",PDF_on);
    SET_OPT_INT_PARAMETER("use_rho",use_rho);
    if (vars->Model >= 4) 
    {
	SET_OPT_INT_PARAMETER("Pulse_regions",Pulse_regions);
	SET_OPT_INT_PARAMETER("Cluster_pulses",Cluster_pulses);
	SET_OPT_INT_PARAMETER("Sig_test",Sig_test);
	SET_OPT_INT_PARAMETER("yy",yy);
	SET_OPT_PARAMETER("betae_init",betae_init);
	SET_PARAMETER("eclipse_init",eclipse_init);
    }
    if (vars->Model >= 5) 
    {
	SET_OPT_INT_PARAMETER("Model_0",Model_0);
	SET_PARAMETER("gamma_init",gamma_init);
	SET_PARAMETER("gamma_high",gamma_high);
	SET_PARAMETER("gamma_low",gamma_low);
	SET_PARAMETER("gamma_mean",gamma_mean);
	SET_PARAMETER("Cmax_init",Cmax_init);
	SET_PARAMETER("Cmax_mean",Cmax_mean);
	SET_OPT_PARAMETER("Cmax_0",Cmax_0);
	SET_PARAMETER("IC50_init",IC50_init);
	SET_PARAMETER("IC50_mean",IC50_mean);
	SET_PARAMETER("m_init",m_init);
	SET_PARAMETER("m_low",m_low);
	SET_PARAMETER("m_high",m_high);
	SET_PARAMETER("m_mean",m_mean);
	SET_PARAMETER("bolus",bolus);
	SET_PARAMETER("absorb_init",absorb_init);
	SET_PARAMETER("absorb_low",absorb_low);
	SET_PARAMETER("absorb_high",absorb_high);
	SET_PARAMETER("absorb_mean",absorb_mean);

	SET_OPT_INT_PARAMETER("absorb_hrs",absorb_hrs);
	SET_OPT_INT_PARAMETER("gamma_hrs",gamma_hrs);

	SET_PARAMETER("beta_mean",beta_mean);
	SET_OPT_PARAMETER("betae_mean",betae_mean);
	SET_PARAMETER("c_mean",c_mean);
	SET_PARAMETER("latent_inf_mean",latent_inf_mean);
	SET_PARAMETER("delta_mean",delta_mean);
	SET_PARAMETER("fpos_mean",fpos_mean);
	SET_OPT_PARAMETER("hill_mean",hill_mean);
	SET_PARAMETER("an_mean",an_mean);
	SET_PARAMETER("rho_mean",rho_mean);
	SET_PARAMETER("theta_mean",theta_mean);
	SET_PARAMETER("r_mean",r_mean);
	SET_PARAMETER("eclipse_mean",eclipse_mean);
	SET_PARAMETER("log_p_mean",log_p_mean);

	SET_PARAMETER("cd8_ic50_mean",cd8_ic50_mean);
	SET_PARAMETER("cd8_ic50_init",cd8_ic50_init);
	SET_PARAMETER("cd8_ic50_std",cd8_ic50_std);

	SET_OPT_PARAMETER("exp_days_init",exp_days_init);
	SET_OPT_PARAMETER("exp_days_mean",exp_days_mean);
	SET_OPT_PARAMETER("exp_days_std",exp_days_std);

	SET_PARAMETER("kappa_init",kappa_init);
	SET_PARAMETER("kappa_mean",kappa_mean);
	SET_PARAMETER("kappa_std",kappa_std);

	SET_PARAMETER("beta_std",beta_std);
	SET_OPT_PARAMETER("betae_std",betae_std);
	SET_PARAMETER("c_std",c_std);
	SET_PARAMETER("latent_inf_std",latent_inf_std);
	SET_PARAMETER("delta_std",delta_std);
	SET_PARAMETER("fpos_std",fpos_std);
	SET_OPT_PARAMETER("hill_std",hill_std);
	SET_PARAMETER("an_std",an_std);
	SET_PARAMETER("rho_std",rho_std);
	SET_PARAMETER("theta_std",theta_std);
	SET_PARAMETER("r_std",r_std);
	SET_PARAMETER("eclipse_std",eclipse_std);
	SET_PARAMETER("log_p_std",log_p_std);
    }
    if (vars->Model == 8) 
    {
	SET_PARAMETER("kD",kD);
	SET_INT_PARAMETER("cT",cT);
	SET_INT_PARAMETER("nT",nT);
    }
}
void scaleTcells(state_var* Tp[], double scaleFactor, bool globally, globalState *vars) 
{
    if (globally)
    {
	unsigned long long Ttot=0;
	for (int i=0; i < vars->Regions; i++)
	    Ttot+=Tp[i]->val();

	double avgChange= Ttot*(scaleFactor-1.0)/vars->Regions;
	for (int i=0; i < vars->Regions; i++)
	    Tp[i]->update_with(avgChange);
	
    }
    else
	for (int i=0; i < vars->Regions; i++)
	    Tp[i]->initialz_with(Tp[i]->val()*scaleFactor);
}

/* This function returns a"score" for the simulations runs based on a given
 * set of parameter values (for now this is the % of 50 criteria in 95% CI ranges */
double ScoreFunction(int *valid, 
	void *data, FILE *results)
{ 
    double score =0.;

    double score1=0.,score2=0.,score3=0.;

    *valid=0;

    globalState *vars = (globalState *)data;

    int Nepis = 10000;

    double shed_thresh=100.; /* level of detection for virus shedding episode */

    int model;

    double Episode_limit, T1Sim, T0_Sim;
    double TSim, tstep, timep, period_p;

    double swabT=0., max_T=0., First_T=0., Last_T=0.,cont_First_T=0.;
    bool inContEpisode=false;

    double start_inf_time[MAX_HEXCELLS];
    double start_R0[MAX_HEXCELLS];

    double totPlaques[Nepis]; 

    double episodDur[Nepis], 
	    maxFirstReg[Nepis];

    double period_episodDur[Nepis]; 

    double Repro[MAX_HEXCELLS];
    double logRepro[MAX_HEXCELLS];
    double ReproMax[MAX_HEXCELLS];
    double ReproMin[MAX_HEXCELLS];
	    
    double ViLifeMax[MAX_HEXCELLS];
    double ViLifeMin[MAX_HEXCELLS];
    double ViLifeSpan;

    int swabs[VL_BINS];
    int criteria1[VL_BINS];
    double crit1perc[VL_BINS];

    int peaks[EPI_BINS];
    int criteria2[EPI_BINS];
    double crit2perc[EPI_BINS];

    int criteria3[DUR_BINS];
    double crit3perc[DUR_BINS];

    double timeTo100[MAX_HEXCELLS];
    double timeTo1000[MAX_HEXCELLS];
    double timeTo10k[MAX_HEXCELLS];
    double timeTo100k[MAX_HEXCELLS];
    double timeToPeak[MAX_HEXCELLS];

    double timeTo100Deaths[MAX_HEXCELLS];
    double timeTo1000Deaths[MAX_HEXCELLS];
    double timeTo10kDeaths[MAX_HEXCELLS];
    double timeTo100kDeaths[MAX_HEXCELLS];

    double timeAtPeak[MAX_HEXCELLS];
    double timeFromPeak[MAX_HEXCELLS];
    unsigned long int iThisAtPeak[MAX_HEXCELLS];

    double timeAtGlobalPeak=-1;


    int VL_bin, maxVL_bin;

    double beta, beta_e, theta, delta, c, r, latent_inf, p, avgVL; 
    double rho, inf, rinf, eclipse;

    unsigned long int pastVL;
    unsigned long int currentVL=0;
    unsigned long int measuredVL;

    unsigned long int  T0=0;
    unsigned long int  I0=0;
    unsigned long int  Ve0=0;

    int del1;

    int  counter, swabsThisEpi; 
    int  cont_counter;
    int  period_counter;
    double  pos_this_period;

    int totalEpis;
    int cont_totalEpis;
    int totalSwabs;
    int totalPosSwabs;
    int posSwabs;
    int withPosSwabs=0;
    int thisSubjSwabs=0;
    int shedderSwabs=0;

    double patSwabMeds[Nepis];	// Track median pos swab VL for each patient (report avg)
    double patSwabVars[Nepis];	// Track variance of swabs for each patient (report avg)

    double totPosSwabs[MAX_SWABS];
    double totPeaks[MAX_SWABS];
    double totDurations[MAX_SWABS];
    double thisPosSwabs[MAX_SWABS];

    double time_in_days;
    double total_stats_time = 0;
    double stats_time = 0;

    unsigned int Iplaquebirths =0; /* used for model 5 plaques */
    int plaque_births =0; /* used for model 5 plaques */

    /* state variables */
    state_var *Sp;
    state_var *Ssub[MAX_HEXCELLS];
    state_var *Ip[MAX_HEXCELLS], *Tp[MAX_HEXCELLS];
    state_var *pre_Tp[MAX_HEXCELLS];
    state_var *Ve[MAX_HEXCELLS];
    state_var *Ec[MAX_HEXCELLS], *Veadj[MAX_HEXCELLS], *Vi[MAX_HEXCELLS];
    state_var *Ithis[MAX_HEXCELLS], *Vithis[MAX_HEXCELLS], *Vethis[MAX_HEXCELLS];
    state_var *Id_this[MAX_HEXCELLS];

    double *actArray;

    int episodeStop=0;

    int totalPlaques; /* plaques generated by an episode */
    int plaqColors[MAX_HEXCELLS];
    int selectedRegions[MAX_HEXCELLS];

    gsl_matrix *pastTps;

    /* model 6 & 7 drug related parameters */

    double gamma=0;
    double Cmax=0;
    double Cmax_0=0;
    double IC50=0;
    double m=0;

    /* model 6 & 7 drug related constants */
    double absorb=0;	/* time to absorb bolus */
    double bolus=0;	/* bolus interval (ex. 12 hrs or 0.5 days) */

    /* calculated drug vars */
    double ACV=0;
    double p_ACV=0;
    double latent_inf_ACV=0;

    double infuse=0;
    double lastBolus=0;
    int firstBolus=0;
    int doses=0;

    int Nr_count, N_runs, Total_epis, T0_files;
    int pcounter = 0;

    int epi_1mm = 0;
    int epi_2mm = 0;

    double time_over_1 = 0;
    double time_over_2 = 0;

    char t0_default[] = "hhv8_sim.T0";
    char t0_file[100];
    long time_t0;
    long time_now;

    FILE *t0Fp;
    char tmpline[MAX_LINE];
    int reg=0;
    unsigned int tp_val;
    unsigned long int Tptot = 0;

    char I0_file[] = "hhv8_sim.I0";
    FILE *I0Fp;
    unsigned long int Ip_val;
    unsigned long int Iptot = 0;

    int startRegion=-1;
    unsigned long int maxStartRegionVe = 0;	/* max extracellular virons in start region */
    double maxStartTime = 0;	

    int firstRegion; /* tracks region 1st infected for an episode */
    int infRegions=0; /* used when threshold is set */

    int firstPass=1;

    unsigned long int VethisTot = 0;	/* extracellular viral births this episode */
    unsigned long int VithisTot = 0;	/* intracellular viral births this episode */
    unsigned long int IthisTot = 0;	/* infected cells this episode */
    unsigned long int IdthisTot = 0;	/* dead cells this episode */
    unsigned long int Itot_p = 0;	/* infected cells this episode (past value)*/
    double Tinf = 0;	/* Tcell infused */
    double Tbump = 0;	/* Tcell jump due to rho */
    double Tdump = 0;	/* Tcell loss due to rho */
    double Tkills = 0;	/* Tcell kills */
    double Ideaths = 0;	/* natural deaths */
    unsigned long int newTcells = 0;	/* newly created Tcells */
    unsigned long int newInfcells = 0;	/* newly infected cells (productive)*/
    unsigned long int newInfevents = 0;	/* newly infected cells (from spread)*/
    unsigned long int newVirons = 0;	/* new Virons (Ve)*/
    long unsigned int maxVL[Nepis];

    long unsigned int cd8size=0;
    long unsigned int vet=0;
    long unsigned int vet_p=0;
    long unsigned int vit=0;
    long unsigned int infTot=0;
    long unsigned int infGlobalMax=0;
    long unsigned int infLocalMax[MAX_HEXCELLS];
    long unsigned int Ip_p[MAX_HEXCELLS];
    long unsigned int Ithis_p[MAX_HEXCELLS];
    long unsigned int Id_p[MAX_HEXCELLS];
    long unsigned int Ve_to_date = 0;
    long unsigned int Inf_to_date = 0;
    long unsigned int Spread_to_date = 0;
    long unsigned int T_to_date = 0;
    long unsigned int T_inf_tot = 0;
    long unsigned int T_bump_tot = 0;
    long unsigned int T_dump_tot = 0;
    long unsigned int T_kills_tot = 0;
    long unsigned int Ideaths_tot = 0;
    double running_cd8s_max=0;
    double running_cd8s_mean=0;
    double running_cd8s_std=0;
    double running_repro_under1=0;
    double running_repro_max=0;
    double running_repro_mean=0;
    double running_log_repro_mean=0;
    double running_repro_std=0;
    double running_log_repro_std=0;
    double cd8s_std=0;
    double repro_under1=0;
    double repro_max=0;
    double repro_mean=0;
    double log_repro_mean=0;
    double repro_std=0;
    double log_repro_std=0;
    int hhv8_regs=0;
    int inf_regs=0;
    int curr_hhv8_regs=0;
    int curr_inf_regs=0;
    double avg_hhv8_regs=0;
    double avg_inf_regs=0;
    int extra_plaque_starts=0;
    int steps=0;
    int cd8_reexpansions=0;

    double period_AUC=0; /* cumulative area under log Ve curve (this period)*/


    if (results == NULL)
	results=stdout;

    /* set state variables using passed state */

    Episode_limit = vars->Episode_limit; 
    tstep = vars->tstep;
    N_runs = vars->N_runs;
    Total_epis = vars->Total_epis;
    T0_files = vars->T0_files;
    T0 = vars->To; 
    del1 = vars->del1; 
    T1Sim = 0.0;

    strcpy(t0_file, t0_default);

    vars->sample_index = 0;

    actArray = (double *)malloc(N_runs*sizeof(double));

    /* set parameters using passed param vector */
    beta=vars->beta_init;
    latent_inf=vars->latent_inf_init;
    p = pow(10,vars->log_p_init);
    c=vars->c_init;
    theta=vars->theta_init;
    inf=vars->inf_init;
    r=vars->r_init;
    rinf=vars->rinf_init;
    delta=vars->delta_init;
    beta_e=vars->betae_init;
    rho=vars->rho_init;
    eclipse=vars->eclipse_init;

    if (vars->Model >= 5) /* get drug related params */
    {
	Cmax=vars->Cmax_init;
	IC50=vars->IC50_init;
	m=vars->m_init;
	if (vars->gamma_hrs)
	    gamma=log(2) / (vars->gamma_init/ 24.);
	else
	    gamma=log(2)/vars->gamma_init; 

	if (vars->absorb_hrs)
	    absorb=(vars->absorb_init/ 24.);
	else
	    absorb=vars->absorb_init; 

	bolus = vars->bolus;
	if (vars->Cmax_0 == 0.0)
	    Cmax_0 = 0.0;
	else
	    Cmax_0 = vars->Cmax_0;
    }

    pastTps=gsl_matrix_alloc(vars->Regions,del1);

    for (int i=0; i < vars->Regions; i++)
	for (int j=0; j < del1; j++)
	    gsl_matrix_set(pastTps,i,j,T0/vars->Regions);
	
    Sp = new state_var();
    Sp->initialz_with(vars->N0); 
    for (int j=0; j <vars->Regions && j < MAX_HEXCELLS; j++)
    {
	// initialize all state variables
	Ssub[j] = new state_var();
	Ssub[j]->initialz_with(vars->N0/vars->Regions); 
	Ec[j] = new state_var();
	Ip[j] = new state_var();
	Tp[j] = new state_var();
	pre_Tp[j] = new state_var();
	pre_Tp[j]->initialz_with(0);
	Ithis[j] = new state_var();
	Id_this[j] = new state_var();
	Vithis[j] = new state_var();
	Vethis[j] = new state_var();
	Vi[j] = new state_var();
	Ve[j] = new state_var();
	Veadj[j] = new state_var();
	plaqColors[j] = 0;
	if (vars->Regions > 0 && j < vars->Pulse_regions && !vars->Cluster_pulses)
	{
	    selectedRegions[j] = (int)gsl_rng_uniform_int(vars->ur,vars->Regions);
	    /* when using "neighbor effects" plaques shouldn't start at the edges */
	    if (vars->Model >= 5)
	    {
		while (vars->cells[selectedRegions[j]]->num_neighbors != 6)
		    selectedRegions[j] = (int)gsl_rng_uniform_int(vars->ur,vars->Regions);
	    }
	}
    }
    Nr_count = 1;

    for (int i=0; i < N_runs; i++)
	actArray[i]=0;

    totalEpis=0; 
    cont_totalEpis=0; 
    totalSwabs=0; 
    totalPosSwabs=0;
    posSwabs=0; 
    shedderSwabs=0; 
    thisSubjSwabs=0;

    //
    //empty cumulative sums for each summary measure
    for(int j=0;j< VL_BINS;j++)
	criteria1[j] = 0;
    
    for(int j=0;j< EPI_BINS;j++)
	criteria2[j] = 0;
    
    for(int j=0;j< DUR_BINS;j++)
	criteria3[j] = 0;
    
    //go through patients (runs) - stats should be for combined totals!
    while (Nr_count <= N_runs && ((Total_epis == 0) || 
	   (Total_epis > 0 && totalEpis < Total_epis)))
    {
	if (Nr_count > 1 && (vars->PDF_on==2 || vars->alt_inp_file != NULL))
	{
	/* if we read in an alternate inpuit file last run, then re-read the original input file! */
	    if (vars->alt_inp_file != NULL)
	    {
		fprintf(results,"Re-reading original input file!\n");

		read_input_file(1, vars->inp_file,vars);
	    }

	    // (no change allowed if PDF_on==1)
	    if (vars->PDF_on==2)
	    {
		do
		    vars->beta_init = vars->beta_mean+gsl_ran_gaussian (vars->ur, vars->beta_std);
		while (vars->beta_init < 0);
		do
		    vars->betae_init = vars->betae_mean+gsl_ran_gaussian (vars->ur, vars->betae_std);
		while (vars->betae_init < 0);
		do
		    vars->latent_inf_init = vars->latent_inf_mean+gsl_ran_gaussian (vars->ur, vars->latent_inf_std);
		while (vars->latent_inf_init < 0);
		do
		    vars->log_p_init = vars->log_p_mean+gsl_ran_gaussian (vars->ur, vars->log_p_std);
		while (vars->log_p_init < vars->log_p_low || vars->log_p_init > vars->log_p_high);
		do
		    vars->c_init = vars->c_mean+gsl_ran_gaussian (vars->ur, vars->c_std);
		while (vars->c_init < 0);
		do
		    vars->theta_init = vars->theta_mean+gsl_ran_gaussian (vars->ur, vars->theta_std);
		while (vars->theta_init < 0);
		do
		    vars->r_init = vars->r_mean+gsl_ran_gaussian (vars->ur, vars->r_std);
		while (vars->r_init < 0);
		do
		    vars->rinf_init = vars->rinf_mean+gsl_ran_gaussian (vars->ur, vars->rinf_std);
		while (vars->rinf_init < 0);
		do
		    vars->inf_init = vars->inf_mean+gsl_ran_gaussian (vars->ur, vars->inf_std);
		while (vars->inf_init < 0);

		do
		    vars->delta_init = vars->delta_mean+gsl_ran_gaussian (vars->ur, vars->delta_std);
		while (vars->delta_init < 0);

		do
		    vars->eclipse_init = vars->eclipse_mean+gsl_ran_gaussian (vars->ur, vars->eclipse_std);
		while (vars->eclipse_init < 0);

		do
		    vars->fpos = vars->fpos_mean + gsl_ran_gaussian(vars->ur,vars->fpos_std);
		while (vars->fpos < 0);

		do
		    vars->hill = vars->hill_mean + gsl_ran_gaussian(vars->ur,vars->hill_std);
		while (vars->hill < 0);

		do
		    vars->an = vars->an_mean + gsl_ran_gaussian(vars->ur,vars->an_std);
		while (vars->an < 0);

		do
		    vars->alpha_init = vars->alpha_mean + gsl_ran_gaussian(vars->ur,vars->alpha_std);
		while (vars->alpha_init < 0);

		do
		    vars->kappa_init = vars->kappa_mean + gsl_ran_gaussian(vars->ur,vars->kappa_std);
		while (vars->kappa_init < 0);

		do
		    vars->exp_days_init = vars->exp_days_mean + gsl_ran_gaussian(vars->ur,vars->exp_days_std);
		while (vars->exp_days_init < 0);

		do
		    vars->cd8_ic50_init = vars->cd8_ic50_mean + gsl_ran_gaussian(vars->ur,vars->cd8_ic50_std);
		while (vars->cd8_ic50_init < 0);

		do
		    vars->To = vars->To_mean + gsl_ran_gaussian(vars->ur,vars->To_std);
		while (vars->To < 0);

		if (vars->Model > 5 || vars->Model_2 > 5)
		{
		    do
			vars->Cmax_init = vars->Cmax_mean+gsl_ran_gaussian (vars->ur, vars->Cmax_std);
		    while (vars->Cmax_init < 0);

		    do
			vars->IC50_init = vars->IC50_mean+gsl_ran_gaussian (vars->ur, vars->IC50_std);
		    while (vars->IC50_init < 0);

		    do
			if (vars->m_mean != 0)
			    vars->m_init = vars->m_mean+gsl_ran_gaussian (vars->ur, vars->m_std);
			else
			    vars->m_init = vars->m_low +gsl_rng_uniform(vars->ur)* (vars->m_high-vars->m_low);
		    while (vars->m_init < 1);

		    do
			if (vars->gamma_mean != 0)
			    vars->gamma_init = vars->gamma_mean+gsl_ran_gaussian (vars->ur, vars->gamma_std);
			else
			    vars->gamma_init = vars->gamma_low +gsl_rng_uniform(vars->ur)* (vars->gamma_high-vars->gamma_low);
		    while (vars->gamma_init < 0);

		    /* absorb handled with range rather than distribution */
		    do
			if (vars->absorb_mean != 0)
			    vars->absorb_init = vars->absorb_mean+gsl_ran_gaussian (vars->ur, vars->absorb_std);
			else
			    vars->absorb_init = vars->absorb_low +gsl_rng_uniform(vars->ur)* (vars->absorb_high-vars->absorb_low);
		    while (vars->absorb_init < 0);
		}
	    }

	    beta=vars->beta_init;
	    latent_inf=vars->latent_inf_init;
	    p = pow(10,vars->log_p_init);
	    c=vars->c_init;
	    theta=vars->theta_init;
	    inf=vars->inf_init;
	    r=vars->r_init;
	    rinf=vars->rinf_init;
	    delta=vars->delta_init;
	    beta_e=vars->betae_init;
	    rho=vars->rho_init;
	    eclipse=vars->eclipse_init;


	    if (vars->Model > 5) /* adjust drug related params */
	    {
		Cmax=vars->Cmax_init;
		IC50=vars->IC50_init;
		m=vars->m_init;
		if (vars->gamma_hrs)
		    gamma=log(2) / (vars->gamma_init/ 24.);
		else
		    gamma=log(2)/vars->gamma_init; 

		if (vars->absorb_hrs)
		    absorb=(vars->absorb_init/ 24.);
		else
		    absorb=vars->absorb_init; 

		bolus = vars->bolus;
		if (vars->Cmax_0 == 0.0)
		    Cmax_0 = 0.0;
		else
		    Cmax_0 = vars->Cmax_0;
	    }
	}
	Iptot = 0;
	episodeStop=0;
	reg = 0;
	// Read T0 values from saved file (if available)
	//
	Tptot = 0;
	reg = 0;

	// Read random T0 file if Calc_T0 == -1 
	if (vars->Calc_T0 < 0)
	{
	    int randfile=(int)gsl_rng_uniform_int(vars->ur,T0_files);
	    sprintf(t0_file, "%s.%d",t0_default,randfile);
	}

	if ((vars->Calc_T0 <= 0 || vars->Calc_T0==3) && 
		(t0Fp = fopen (t0_file,"r")) != NULL)
	{
	    char *cp = NULL;
	    while ((cp = fgets(tmpline, MAX_LINE-1,t0Fp)) != NULL && 
		    reg < vars->Regions && reg < MAX_HEXCELLS) {
		tmpline[MAX_LINE-1] = '\0';
		
		if (sscanf(tmpline,"%u",&tp_val) != 1) {
		    cerr << "Error while reading T0 value for region #"<<reg<<"t0 file "<<t0_file<<"\n";
		    exit(1);
		}
		Tp[reg]->initialz_with(tp_val) ; 
		Tptot += Tp[reg]->val();
		if (vars->Verbose)
		{
		    cout <<"T0["<< reg <<"]="<<tp_val;
		    cout<<endl;
		}
		reg++;
	    }
	    if (reg != vars->Regions)
		fprintf(stderr,
		    "Warning: number of T0 values (%d) in file %s not the same as number of regions (%d)!\n",
		    reg,t0_file,vars->Regions);
	    else if (cp != NULL)
		fprintf(stderr,
		    "Warning: additional T0 values in file %s should be same as #regions (%d)!\n",
		    t0_file,vars->Regions);
			
		    
	    cout<<endl;
	    cout <<"Finished reading T0 init file "<<t0_file<<".\n";
	    fprintf(stderr,"Total T0 over %d regions read in is %lu\n",vars->Regions,Tptot);
	    fclose (t0Fp);
	}
	else
	{
	    // initialize T0 to random count between 3000 and 10000 (per region)
	    if (vars->Model >= 4)
	    {
		for (int j=0; j <vars->Regions; j++)
		{
		    Tp[j]->initialz_with((int)gsl_rng_uniform_int(vars->ur,(vars->To/vars->Regions)*0.7) + (vars->To/vars->Regions)*0.3);
		    Tptot += Tp[j]->val();
		}
		if (vars->Verbose)
		    fprintf(stderr,"Total T0 over %d regions randomly generated is %lu\n",vars->Regions,Tptot);
	    }
	    else
		Tp[0]->initialz_with(vars->To);

	    if (vars->Calc_T0 != 0)
	    {
		Sp->initialz_with(vars->N0); 
		// initialize all state variables
		for (int j=0; j <vars->Regions; j++)
		{
		    Ssub[j]->initialz_with(vars->N0/vars->Regions); 
		    Ec[j]->initialz_with(0); 
		    Ip[j]->initialz_with(0); 
		    Ve[j]->initialz_with(0); 
		    Vi[j]->initialz_with(0); 
		    Vethis[j]->initialz_with(0); 
		    Vithis[j]->initialz_with(0); 
		    Ithis[j]->initialz_with(0); 
		    Id_this[j]->initialz_with(0); 
		    Veadj[j]->initialz_with(0); 
		    Repro[j] = 0.0;
		    logRepro[j] = 0.0;
		}
		time_in_days=0.; 
		stats_time=0.; 
		time_over_1=0.; 
		time_over_2=0.; 
		swabT=0.;
		VethisTot = 0;
		VithisTot = 0;
		IthisTot = 0;
		IdthisTot = 0;
		newTcells=0;
		Tinf=0;
		Tbump=0;
		Tkills=0;
		Ideaths=0;
		newInfcells=0;
		newInfevents=0;
		newVirons=0;
		Itot_p = 0;

		if (vars->Model >= 4)
		{
		    ACV = 0;
		    lastBolus = -vars->bolus;
		    firstBolus = 1;
		    infuse = 0;
		    doses=0;
		}
		// run multi day (year) simulation to get mean Tp to use as To
		// actually uses T0_sim +- 40% (randomnly determined)

		// if initialization model is 0, use runtime model
		if (vars->Model_0 == 0)
		    vars->Model_0 = vars->Model;

		T0_Sim = 0.6*vars->T0_Sim +gsl_rng_uniform(vars->ur)*vars->T0_Sim*0.8;
		if (theState->Verbose > 1)
		{
		    fprintf(stderr,"[%d]:Calculating T0 by running for %lf days using Model %d and...\n",Nr_count,T0_Sim, vars->Model_0);
		    if (vars->Model_0 > 5)
			fprintf(stderr,
			    "To=%lf, beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf absorb=%lf gamma=%lf Cmax=%lf IC50=%lf m=%lf regions=%d swabInt=%lf\n",
			    vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,absorb,gamma,Cmax,IC50,m,vars->Regions,vars->swabInterval);
		    else
			fprintf(stderr,
			    "To=%lf, beta=%e latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf time=%lf\n",
			    vars->To,beta,latent_inf,vars->log_p_init,c,theta,delta,r,inf,rinf,rho,eclipse,vars->fpos,vars->an,time_in_days);
		}

		if (vars->Calc_T0 < 2)
		{
		    while (time_in_days<= T0_Sim) 
		    {
			if (vars->Model_0 == 4 || vars->Model_0 == 5)
			{
			    if (!model5(time_in_days,Sp, Ssub, Ec, Ip, 
				    Tp, pastTps, Vi, Ve, 
				    Veadj, &Iplaquebirths, &totalPlaques,
				    plaqColors, vars, 
				    Repro, logRepro, p, latent_inf,
				    Vethis, Vithis, Ithis, Id_this, start_inf_time,
				    &repro_max, &repro_mean, &log_repro_mean,
				    &repro_std, &log_repro_std,
				    &Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newInfevents, &newVirons,
				    &IthisTot,&IdthisTot,&VithisTot,&VethisTot,selectedRegions))
				return 0;
			}
			else if (vars->Model_0 == 6 || vars->Model_0 == 7)
			{

			    double tempCmax;
			    if (firstBolus && Cmax_0 != 0)
			    {
				tempCmax = Cmax_0;
			        ACV += deltaACV(time_in_days, tstep, tempCmax, absorb, ACV, gamma, bolus, &infuse, &firstBolus, &lastBolus, &doses,vars->Total_doses);
				if (vars->Verbose)
				    fprintf(stderr, "ACV level is %lf (Cmax=%lf)\n",ACV,Cmax_0);
			    }
			    else
			    {
				tempCmax = Cmax;
			        ACV += deltaACV(time_in_days, tstep, tempCmax, absorb, ACV, gamma, bolus, &infuse, &firstBolus, &lastBolus, &doses,vars->Total_doses);
			    }

			    if (ACV < 0)
				ACV = 0.;

			    double denom = 1 + pow((ACV/IC50),m);
			    p_ACV = p / denom; 

			    if (vars->Model_0 == 7)
				latent_inf_ACV = latent_inf / denom;
			    else
				latent_inf_ACV = latent_inf;

			    if (!model5(time_in_days,Sp, Ssub, Ec, Ip, 
				    Tp, pastTps, Vi, Ve, Veadj, &Iplaquebirths, 
				    &totalPlaques, plaqColors, vars, 
				    Repro,logRepro,p_ACV,latent_inf_ACV,
				    Vethis, Vithis, Ithis, Id_this, start_inf_time,
				    &repro_max, &repro_mean, &log_repro_mean,
				    &repro_std, &log_repro_std,
				    &Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newInfevents, &newVirons,
				    &IthisTot,&IdthisTot,&VithisTot,&VethisTot,selectedRegions))
				return 0;
			}
			else if (vars->Model_0 == 8)
			{

			    double tempCmax;
			    if (firstBolus && Cmax_0 != 0)
				tempCmax = Cmax_0;
			    else
				tempCmax = Cmax;
			    ACV += deltaACV(time_in_days, tstep, tempCmax, absorb, ACV, gamma, bolus, &infuse, &firstBolus, &lastBolus, &doses,vars->Total_doses);

			    if (ACV < 0)
				ACV = 0.;

			    double factor = 1 - siliciano(vars->nT,vars->cT,vars->kD,ACV);

			    if (firstPass)
			    {
				fprintf(stderr,"Siliciano factor = %lf\n",
				    factor);
				firstPass=0;
			    }
			    p_ACV = p * factor; 
			    latent_inf_ACV = latent_inf * factor;

			    if (!model5(time_in_days,Sp, Ssub, Ec, Ip, 
				    Tp, pastTps, Vi, Ve, Veadj, &Iplaquebirths, 
				    &totalPlaques, plaqColors, vars, 
				    Repro,logRepro,p_ACV,latent_inf_ACV,
				    Vethis, Vithis, Ithis, Id_this, start_inf_time,
				    &repro_max, &repro_mean, &log_repro_mean,
				    &repro_std, &log_repro_std,
				    &Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newInfevents, &newVirons,
				    &IthisTot,&IdthisTot,&VithisTot,&VethisTot,selectedRegions))
				return 0;
			}

			time_in_days = time_in_days + tstep;
			if (time_in_days - swabT > vars->swabInterval) 
			{
			    /* keep T0 as running avg */
			    swabT=time_in_days;
			}
			T0 = 0;
			I0 = 0;
			Ve0 = 0;
			for (int j=0; j <vars->Regions; j++)
			{
			    T0 += Tp[j]->val();
			    I0 += Ip[j]->val();
			    Ve0 += Ve[j]->val();
			    if (Ip[j]->val() > 0)
				inf_regs++;
			    if (Ve[j]->val() > 0)
				hhv8_regs++;
			}
		    }
		    if (vars->Verbose > 1)
			fprintf(stderr,
			    "[%d]: Calculated T0=%lu (after %lf days)\n",
			    Nr_count,T0,time_in_days);
		    if (vars->Calc_T0 < 0)
		    {
			if( (t0Fp = fopen(t0_file,"wt")) == NULL){
			    cerr << "Could not open T0 file "<<t0_file<<" for writing.\n";
			}
			else
			{
			    for (int i=0; i <vars->Regions; i++)
				fprintf(t0Fp,"%lu\n",Tp[i]->val());
			    fclose(t0Fp);
			    time(&time_t0);
			    if (vars->Verbose > 1)
				fprintf(stderr,"(wrote T0 file %s (elapsed=%ld)\n",
					t0_file,time_t0-vars->time_st);
			}
		    }
		}
	    }
	}

	curr_inf_regs=0;
	curr_hhv8_regs=0;
	for (int j=0; j <vars->Regions; j++)
	{
	    if (Ip[j]->val() > 0)
		curr_inf_regs++;
	    if (Ve[j]->val() > 0)
		curr_hhv8_regs++;
	}


	if (T0 > 1e8)
	{
	    fprintf(results,
		"Initial T0 is too high %lu (inf reg=%d (%d at end), hhv8 reg=%d (%d at end))! Skipping simulation!\n",
		T0,inf_regs,curr_inf_regs,hhv8_regs,curr_hhv8_regs);
	    return 0.;
	}
	else
	    fprintf(results,"Initial T0 is %lu (inf reg=%d (%d at end), hhv8 reg=%d (%d at end)). Beginning simulation!\n",T0,inf_regs,curr_inf_regs,hhv8_regs,curr_hhv8_regs);
	/* reset T0 values in past value matrix
	for (int i=0; i < vars->Regions; i++)
	    for (int j=0; j < del1; j++)
		gsl_matrix_set(pastTps,i,j,Tp[j]->val());
	*/
	
	vars->sample_index = 0;

	Sp->initialz_with(vars->N0); 
	// initialize all state variables
	// (Leave Ve, Ip and Tp at burn-in levels)
	for (int j=0; j <vars->Regions; j++)
	{
	    Ssub[j]->initialz_with(vars->N0/vars->Regions); 
	    Ec[j]->initialz_with(0); 
	    Ip[j]->initialz_with(0); 
	    Ve[j]->initialz_with(0); 
	    Vi[j]->initialz_with(0); 
	    Vethis[j]->initialz_with(Ve[j]->val()); 
	    Vithis[j]->initialz_with(Vi[j]->val()); 
	    Ithis[j]->initialz_with(Ip[j]->val()); 
	    Id_this[j]->initialz_with(0); 
	    Veadj[j]->initialz_with(0); 
	}
	ACV = 0;
	VethisTot = 0;
	VithisTot = 0;
	IthisTot = 0;
	IdthisTot = 0;

	running_repro_under1=0;
	running_cd8s_max=0;
	running_cd8s_mean=0;
	running_cd8s_std=0;
	running_repro_max=0;
	running_repro_mean=0;
	running_log_repro_mean=0;
	running_repro_std=0;
	running_log_repro_std=0;

	cd8s_std=0;
	repro_mean=0;
	log_repro_mean=0;
	repro_std=0;
	log_repro_std=0;

	hhv8_regs=0;
	inf_regs=0;
	extra_plaque_starts=0;
	avg_hhv8_regs=0;
	avg_inf_regs=0;
	steps=0;

	Tinf=0;
	Tbump=0;
	Tdump=0;
	Tkills=0;
	Ideaths=0;
	newTcells=0;
	newInfcells=0;
	newInfevents=0;
	newVirons=0;
	cd8size=0;
	vet=0;
	vet_p=0;
	vit=0;
	totalPlaques=0;
	infTot=0;
	infGlobalMax=0;
	firstRegion = -1;

	if (vars->Model >= 4)
	{
	    lastBolus = -vars->bolus;
	    firstBolus = 1;
	    infuse = 0;
	}

	counter=0; 
	cont_counter=0; 
	period_counter=0; 
	pos_this_period=0; 

	pcounter=0; 
	time_in_days=0.; 
	stats_time=0.; 

	timep=0.; 
	period_p=0.; 
	for(int i=0;i<Nepis;i++)
	{
	    maxVL[i] = 0;
	    maxFirstReg[i] = 0.;
	    episodDur[i] = 0.;
	    period_episodDur[i] = 0.;
	    totPlaques[i] = 0.;
	}
	for(int i=0;i<MAX_HEXCELLS;i++)
	{
	    timeTo100[i] = -1;
	    timeTo1000[i] = -1;
	    timeTo10k[i] = -1;
	    timeTo100k[i] = -1;
	    timeTo100Deaths[i] = -1;
	    timeTo1000Deaths[i] = -1;
	    timeTo10kDeaths[i] = -1;
	    timeTo100kDeaths[i] = -1;
	    timeToPeak[i] = -1;
	    timeFromPeak[i] = -1;
	    start_inf_time[i] = 0;
	    start_R0[i] = 0;
	    Ip_p[i] = 0;
	    Ithis_p[i] = 0;
	    iThisAtPeak[i] = 0;
	    Id_p[i] = 0;
	}

	for(int i=0;i<VL_BINS;i++)
	    swabs[i] = 0;

	for(int i=0;i<EPI_BINS;i++)
	{
	    peaks[i] = 0;
	}
	
	//propagating in time loop
	swabT=0;
	pastVL = 0.;
	swabsThisEpi=0;

	bool inEpisode=false;
	bool any_neg=false;
	double max_dur = 0;
	Iplaquebirths=0;
	plaque_births=0;
	avgVL=0.;

	for (int j=0; j <vars->Regions; j++) {
	    ReproMax[j] = 0;
	    ReproMin[j] = 100.;
	    ViLifeMax[j] = 0;
	    ViLifeMin[j] = 100.;
	}
	if (vars->writeOn && vars->Model >= 4)
	{
	    if (vars->dataF1 != NULL) 
	    {
		fprintf(vars->dataF1,"time");
		fprintf(vars->dataF1,",vet");
		for (int j=0; j <vars->Regions; j++)
		    fprintf(vars->dataF1,",ve[%d]",j+1);
		fprintf(vars->dataF1,"\n");
	    }
	    if (vars->dataF2 != NULL) 
	    {
		fprintf(vars->dataF2,"time");
		fprintf(vars->dataF2,",VethisTot");
		for (int j=0; j <vars->Regions; j++)
		    fprintf(vars->dataF2,",Vethis[%d]",j+1);
		fprintf(vars->dataF2,"\n");
	    }
	    if (vars->dataF3 != NULL) 
	    {
		fprintf(vars->dataF3,"time");
		fprintf(vars->dataF3,",vit");
		for (int j=0; j <vars->Regions; j++)
		    fprintf(vars->dataF3,",vi[%d]",j+1);
		fprintf(vars->dataF3,"\n");
	    }
	    if (vars->dataF4 != NULL) 
	    {
		fprintf(vars->dataF4,"time");
		fprintf(vars->dataF4,",VithisTot");
		for (int j=0; j <vars->Regions; j++)
		    fprintf(vars->dataF4,",Vithis[%d]",j+1);
		fprintf(vars->dataF4,"\n");
	    }
	    if (vars->dataF5 != NULL) 
	    {
		fprintf(vars->dataF5,"time");
		fprintf(vars->dataF5,",infTot");
		for (int j=0; j <vars->Regions; j++)
		    fprintf(vars->dataF5,",inf[%d]",j+1);
		fprintf(vars->dataF5,"\n");
	    }
	    if (vars->dataF6 != NULL) 
	    {
		fprintf(vars->dataF6,"time");
		fprintf(vars->dataF6,",IthisTot");
		fprintf(vars->dataF6,",IdthisTot");
		fprintf(vars->dataF6,"\n");
	    }
	    if (vars->dataF8 != NULL) 
	    {
		fprintf(vars->dataF8,"time");
		for (int j=0; j <vars->Regions; j++)
		    fprintf(vars->dataF8,",Repro[%d]",j+1);
		fprintf(vars->dataF8,"\n");
	    }
	    if (vars->dataF9 != NULL) 
	    {
		fprintf(vars->dataF9,"SLE,");
		fprintf(vars->dataF9,"start time,");
		fprintf(vars->dataF9,"end time,");
		fprintf(vars->dataF9,"region,");
		fprintf(vars->dataF9,"start R0,");
		fprintf(vars->dataF9,"birth100,");
		fprintf(vars->dataF9,"birth1000,");
		fprintf(vars->dataF9,"birth10k,");
		fprintf(vars->dataF9,"birth100k,");
		fprintf(vars->dataF9,"death100,");
		fprintf(vars->dataF9,"death1000,");
		fprintf(vars->dataF9,"death10k,");
		fprintf(vars->dataF9,"death100k,");
		fprintf(vars->dataF9,"toPeak,");
		fprintf(vars->dataF9,"fromPeak,");
		fprintf(vars->dataF9,"duration,");
		fprintf(vars->dataF9,"Local peak,");
		fprintf(vars->dataF9,"Global peak,");
		fprintf(vars->dataF9,"Global peak time,");
		fprintf(vars->dataF9,"Ibirth at local peak\n");
	    }
	    if (vars->dataF10 != NULL) 
	    {
		fprintf(vars->dataF10,"episode,");
		fprintf(vars->dataF10,"duration,");
		fprintf(vars->dataF10,"overall peak,");
		fprintf(vars->dataF10,"plaques,");
		fprintf(vars->dataF10,"first region peak\n");
	    }
	    if (vars->dataF11 != NULL) 
	    {
		fprintf(vars->dataF11,"time,");
		fprintf(vars->dataF11,"T_max");
		fprintf(vars->dataF11,",T_tot");
		fprintf(vars->dataF11,",Ve_max");
		fprintf(vars->dataF11,",Ve_tot");
		fprintf(vars->dataF11,",R_max");
		fprintf(vars->dataF11,",R_mean");
		fprintf(vars->dataF11,",R_stddev");
		fprintf(vars->dataF11,"\n");
	    }
	    if (vars->dataF12 != NULL) 
	    {
		fprintf(vars->dataF12,"time");
		fprintf(vars->dataF12,",model");
		fprintf(vars->dataF12,",bolus");
		fprintf(vars->dataF12,",doses");
		fprintf(vars->dataF12,",ACV");
		fprintf(vars->dataF12,",Infus");
		fprintf(vars->dataF12,",lastBolus");
		fprintf(vars->dataF12,",p");
		fprintf(vars->dataF12,",p_ACV");
		fprintf(vars->dataF12,",latent_inf");
		fprintf(vars->dataF12,",latent_inf_ACV");
		fprintf(vars->dataF12,",vet");
		fprintf(vars->dataF12,"\n");
	    }
	}

	int sig_episodes=0;
	infTot=0;
	infGlobalMax=0;

	firstPass=1;
	// Read I0 values from saved file (if available)
	//
	if (vars->Calc_T0 != 2 && (I0Fp = fopen (I0_file,"r")) != NULL)
	{
	    char *cp = NULL;
	    while ((cp = fgets(tmpline, MAX_LINE-1,I0Fp)) != NULL && 
		    reg < vars->Regions && reg < MAX_HEXCELLS) {
		tmpline[MAX_LINE-1] = '\0';
		
		if (sscanf(tmpline,"%lu",&Ip_val) != 1) {
		    cerr << "Error while reading I0 value for region #"<<reg<<"I0 file "<<I0_file<<"\n";
		    exit(1);
		}
		Ip[reg]->initialz_with(Ip_val) ; 
		Iptot += Ip[reg]->val();
		if (vars->Verbose)
		{
		    cout <<"I0["<< reg <<"]="<<Ip_val;
		    cout<<endl;
		}
		reg++;
	    }
	    if (reg != vars->Regions)
		fprintf(stderr,
		    "Warning: number of I0 values (%d) in file %s not the same as number of regions (%d)!\n",
		    reg,I0_file,vars->Regions);
	    else if (cp != NULL)
		fprintf(stderr,
		    "Warning: additional I0 values in file %s should be same as #regions (%d)!\n",
		    I0_file,vars->Regions);
			
		    
	    if (vars->Verbose)
	    {
		cout<<endl;
		cout <<"Finished reading I0 init file "<<I0_file<<".\n";
		fprintf(stderr,"Total of I0 over %d regions is %lu\n",vars->Regions,Iptot);
	    }
	    fclose (I0Fp);
	}
	else 
	{
	    for (int j=0; j <vars->Regions; j++)
		Ip[j]->initialz_with(vars->Io/vars->Regions) ; 
	}


	if (vars->Calc_T0 >= 2)
	{
	    // include all time intervals just in case
	    if (vars->Calc_T0 == 3)
	    {
		T0_Sim = vars->T0_Sim;
		TSim = vars->T3_Sim + vars->T2_Sim + vars->TSim + T1Sim;
	    }
	    else
		TSim = vars->T3_Sim + vars->T2_Sim + vars->TSim + T1Sim + T0_Sim;

	    if (vars->Model_0 != 0)
	        model = vars->Model_0;
	    else
		model = vars->Model;
	}
	else
	{
	    model = vars->Model;
	    T0_Sim = 0;
	    TSim = vars->T3_Sim + vars->T2_Sim + vars->TSim + T1Sim + T0_Sim;

	}
	if (vars->Verbose > 1)
	    fprintf(stderr, "Run %d, Model %d, Total run will be %lf days (from %lf)\n", Nr_count, model, TSim, time_in_days);

	double AUC=0; /* cumulative area under log Ve curve per episode */
	double total_AUC=0; /* cumulative area under log Ve curve */

	inContEpisode=false;
	inEpisode=false;
	ACV = 0;
	lastBolus = -vars->bolus;
	firstBolus = 1;
	infuse = 0;
	doses = 0;
	Ve_to_date = 0;
	Inf_to_date = 0;
	Spread_to_date = 0;
	T_to_date = 0;
	T_inf_tot = 0;
	T_bump_tot = 0;
	T_dump_tot = 0;
	T_kills_tot = 0;
	Ideaths_tot = 0;

	while (((Total_epis > 0 && totalEpis + counter < Total_epis) || 
			Total_epis==0) &&
		(((Episode_limit > 0 && counter < Episode_limit) || 
		  (Episode_limit == 0 && time_in_days<=TSim )) && counter<Nepis && !episodeStop))
	{
	    if (vars->alt_inp_file != NULL &&
		vars->Input_refresh > 0 && time_in_days>=vars->Input_refresh &&
		time_in_days<=vars->Input_refresh+tstep)
	    {
		fprintf(results,"Reading in %s at t=%lf\n",
			vars->alt_inp_file,time_in_days);
		if (vars->Model > 5)
		    fprintf(results,
			"beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf,absorb=%lf, gamma=%lf,Cmax=%lf,IC50=%lf,m=%lf,time=%lf\n",
			beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,absorb,gamma,Cmax,IC50,m,stats_time);
		else
		    fprintf(results,
			"To=%lf,beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf ,time=%lf\n",
			vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,stats_time);

		read_input_file(1, vars->alt_inp_file,vars);

		// if we are to use mean/stddev then redo parameters after file load in case they changed
		// (no change allowed if PDF_on==1)
		if (vars->PDF_on==2)
		{
		    do
			vars->beta_init = vars->beta_mean+gsl_ran_gaussian (vars->ur, vars->beta_std);
		    while (vars->beta_init < 0);
		    do
			vars->betae_init = vars->betae_mean+gsl_ran_gaussian (vars->ur, vars->betae_std);
		    while (vars->betae_init < 0);
		    do
			vars->latent_inf_init = vars->latent_inf_mean+gsl_ran_gaussian (vars->ur, vars->latent_inf_std);
		    while (vars->latent_inf_init < 0);
		    do
			vars->log_p_init = vars->log_p_mean+gsl_ran_gaussian (vars->ur, vars->log_p_std);
		    while (vars->log_p_init < vars->log_p_low || vars->log_p_init > vars->log_p_high);
		    do
			vars->c_init = vars->c_mean+gsl_ran_gaussian (vars->ur, vars->c_std);
		    while (vars->c_init < 0);
		    do
			vars->theta_init = vars->theta_mean+gsl_ran_gaussian (vars->ur, vars->theta_std);
		    while (vars->theta_init < 0);
		    do
			vars->r_init = vars->r_mean+gsl_ran_gaussian (vars->ur, vars->r_std);
		    while (vars->r_init < 0);
		    do
			vars->rinf_init = vars->rinf_mean+gsl_ran_gaussian (vars->ur, vars->rinf_std);
		    while (vars->rinf_init < 0);
		    do
			vars->inf_init = vars->inf_mean+gsl_ran_gaussian (vars->ur, vars->inf_std);
		    while (vars->inf_init < 0);

		    do
			vars->delta_init = vars->delta_mean+gsl_ran_gaussian (vars->ur, vars->delta_std);
		    while (vars->delta_init < 0);

		    do
			vars->eclipse_init = vars->eclipse_mean+gsl_ran_gaussian (vars->ur, vars->eclipse_std);
		    while (vars->eclipse_init < 0);

		    do
			vars->fpos = vars->fpos_mean + gsl_ran_gaussian(vars->ur,vars->fpos_std);
		    while (vars->fpos < 0);

		    do
			vars->hill = vars->hill_mean + gsl_ran_gaussian(vars->ur,vars->hill_std);
		    while (vars->hill < 0);

		    do
			vars->an = vars->an_mean + gsl_ran_gaussian(vars->ur,vars->an_std);
		    while (vars->an < 0);

		    do
			vars->alpha_init = vars->alpha_mean + gsl_ran_gaussian(vars->ur,vars->alpha_std);
		    while (vars->alpha_init < 0);

		    do
			vars->kappa_init = vars->kappa_mean + gsl_ran_gaussian(vars->ur,vars->kappa_std);
		    while (vars->kappa_init < 0);

		    do
			vars->exp_days_init = vars->exp_days_mean + gsl_ran_gaussian(vars->ur,vars->exp_days_std);
		    while (vars->exp_days_init < 0);

		    do
			vars->cd8_ic50_init = vars->cd8_ic50_mean + gsl_ran_gaussian(vars->ur,vars->cd8_ic50_std);
		    while (vars->cd8_ic50_init < 0);

		    do
			vars->To = vars->To_mean + gsl_ran_gaussian(vars->ur,vars->To_std);
		    while (vars->To < 0);

		    if (vars->Model > 5 || vars->Model_2 > 5)
		    {
			do
			    vars->Cmax_init = vars->Cmax_mean+gsl_ran_gaussian (vars->ur, vars->Cmax_std);
			while (vars->Cmax_init < 0);

			do
			    vars->IC50_init = vars->IC50_mean+gsl_ran_gaussian (vars->ur, vars->IC50_std);
			while (vars->IC50_init < 0);

			do
			    if (vars->m_mean != 0)
				vars->m_init = vars->m_mean+gsl_ran_gaussian (vars->ur, vars->m_std);
			    else
				vars->m_init = vars->m_low +gsl_rng_uniform(vars->ur)* (vars->m_high-vars->m_low);
			while (vars->m_init < 1);

			do
			    if (vars->gamma_mean != 0)
				vars->gamma_init = vars->gamma_mean+gsl_ran_gaussian (vars->ur, vars->gamma_std);
			    else
				vars->gamma_init = vars->gamma_low +gsl_rng_uniform(vars->ur)* (vars->gamma_high-vars->gamma_low);
			while (vars->gamma_init < 0);

			/* absorb handled with range rather than distribution */
			do
			    if (vars->absorb_mean != 0)
				vars->absorb_init = vars->absorb_mean+gsl_ran_gaussian (vars->ur, vars->absorb_std);
			    else
				vars->absorb_init = vars->absorb_low +gsl_rng_uniform(vars->ur)* (vars->absorb_high-vars->absorb_low);
			while (vars->absorb_init < 0);
		    }
		}
		fprintf(results,"After Reading ...\n");

		/* reset parameter vectors incase any changed */
		beta=vars->beta_init;
		latent_inf=vars->latent_inf_init;
		p = pow(10,vars->log_p_init);
		c=vars->c_init;
		theta=vars->theta_init;
		inf=vars->inf_init;
		r=vars->r_init;
		rinf=vars->rinf_init;
		delta=vars->delta_init;
		beta_e=vars->betae_init;
		rho=vars->rho_init;
		eclipse=vars->eclipse_init;

		if (vars->Model > 5) /* get drug related params */
		{
		    Cmax=vars->Cmax_init;
		    IC50=vars->IC50_init;
		    m=vars->m_init;
		    if (vars->gamma_hrs)
			gamma=log(2) / (vars->gamma_init/ 24.);
		    else
			gamma=log(2)/vars->gamma_init; 

		    if (vars->absorb_hrs)
			absorb=(vars->absorb_init/ 24.);
		    else
			absorb=vars->absorb_init; 

		    bolus = vars->bolus;
		    if (vars->Cmax_0 == 0.0)
			Cmax_0 = 0.0;
		    else
			Cmax_0 = vars->Cmax_0;
		}
		if (vars->Model > 5)
		    fprintf(results,
			"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf absorb=%lf, gamma=%lf,Cmax=%lf,IC50=%lf,m=%lf\n",
			beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,absorb,gamma,Cmax,IC50,m);
		else
		    fprintf(results,
			"To=%lf,beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf \n",
			vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse);
	    }
	    /* should be outside of an episode (pastVL < shed_thresh)*/

	    if ((vars->Calc_T0 >= 2 && time_in_days>=T0_Sim) || vars->Calc_T0 < 2)
	    {
		if (vars->Model_3 != 0 && 
			time_in_days>=(T0_Sim + vars->TSim + T1Sim+vars->T2_Sim))
		{
		    if (model != vars->Model_3)
		    {
			fprintf(stdout, "Model set to %d at %lf days\n", 
			    vars->Model_3, time_in_days);
			if (vars->Model_3 < 0 && Nr_count==1)
			    scaleTcells(Tp,1.1,false,vars);
			else if (vars->Model_3 < 0 && Nr_count==2)
			    scaleTcells(Tp,1.25,false,vars);
			else if (vars->Model_3 < 0 && Nr_count==3)
			    scaleTcells(Tp,1.5,false,vars);
			else if (vars->Model_3 < 0 && Nr_count==4)
			    scaleTcells(Tp,1.1,true,vars);
			else if (vars->Model_3 < 0 && Nr_count==5)
			    scaleTcells(Tp,1.25,true,vars);
			else if (vars->Model_3 < 0 && Nr_count==6)
			    scaleTcells(Tp,1.5,true,vars);
		    }
		    model = vars->Model_3;
		}
		else if (vars->Model_2 != 0 && 
			time_in_days>=T0_Sim + vars->TSim)
		{
		    if (vars->Calc_T0 == -2 && T1Sim == 0.)
		    {
			T1Sim = 0.005 * gsl_ran_poisson (vars->ur, 100.0);
			fprintf(stdout, "The model will be set to %d in %lf days\n", 
			    vars->Model_2, T1Sim);
		    }

		    if (model != vars->Model_2 && 
			time_in_days >= T1Sim + vars->TSim)
		    {
			fprintf(stdout, "Model set to %d at %lf days\n", 
			    vars->Model_2, time_in_days);
			model = vars->Model_2;
			if (vars->Model_2 < 0 && Nr_count==1)
			    scaleTcells(Tp,1.1,false,vars);
			else if (vars->Model_2 < 0 && Nr_count==2)
			    scaleTcells(Tp,1.25,false,vars);
			else if (vars->Model_2 < 0 && Nr_count==3)
			    scaleTcells(Tp,1.5,false,vars);
			else if (vars->Model_2 < 0 && Nr_count==4)
			    scaleTcells(Tp,1.1,true,vars);
			else if (vars->Model_2 < 0 && Nr_count==5)
			    scaleTcells(Tp,1.25,true,vars);
			else if (vars->Model_2 < 0 && Nr_count==6)
			    scaleTcells(Tp,1.5,true,vars);
		    }
		}
		else
		{
		    if (model != vars->Model)
			fprintf(stdout, "Model set to %d at %lf days\n", 
			    vars->Model, time_in_days);
		    model = vars->Model;
		}
	    }
	    if (abs(model) <= 5)
	    {
		ACV = 0;
		lastBolus = -vars->bolus;
		firstBolus = 1;
		infuse = 0;
		doses = 0;
	    }

	    if (model == 4 || abs(model) == 5)
	    {
		pastVL = 0;
		for (int j=0; j <vars->Regions; j++)
		    pastVL += Ve[j]->val();

		if (!model5(time_in_days,Sp, Ssub, Ec,Ip,Tp,pastTps,  
			Vi, Ve, Veadj, &Iplaquebirths, &totalPlaques, plaqColors, 
			vars, Repro,logRepro,p,
			latent_inf,
			Vethis, Vithis, Ithis, Id_this, start_inf_time,
			&repro_max, &repro_mean, &log_repro_mean,
			&repro_std, &log_repro_std,
			&Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newInfevents, &newVirons,
			&IthisTot,&IdthisTot,&VithisTot,&VethisTot,selectedRegions))
		    return 0;

		currentVL = 0;
		for (int j=0; j <vars->Regions; j++)
		    currentVL += Ve[j]->val();
	    }
	    else if (abs(model) == 6 || abs(model) == 7)
	    {
		pastVL = 0;
		for (int j=0; j <vars->Regions; j++)
		    pastVL += Ve[j]->val();

		double tempCmax;
		if (firstBolus && Cmax_0 != 0)
		{
		    tempCmax = Cmax_0;
		    ACV += deltaACV(time_in_days, tstep, tempCmax, absorb, ACV, gamma, bolus, &infuse, &firstBolus, &lastBolus, &doses,vars->Total_doses);
		}
		else
		{
		    tempCmax = Cmax;
		    ACV += deltaACV(time_in_days, tstep, tempCmax, absorb, ACV, gamma, bolus, &infuse, &firstBolus, &lastBolus, &doses,vars->Total_doses);
		}

		/*if (ACV > Cmax)
		    ACV = Cmax;*/

		if (ACV < 0)
		    ACV = 0.;


		double denom = 1 + pow((ACV/IC50),m);
		p_ACV = p / denom; 

		if (abs(model) == 7)
		    latent_inf_ACV = latent_inf / denom;
		else
		    latent_inf_ACV = latent_inf;

		if (!model5(time_in_days,Sp, Ssub, Ec,Ip,Tp,pastTps,  
			Vi, Ve, Veadj, &Iplaquebirths, &totalPlaques, plaqColors, 
			vars, Repro,logRepro,p_ACV,
			latent_inf_ACV,
			Vethis, Vithis, Ithis, Id_this, start_inf_time,
			&repro_max, &repro_mean, &log_repro_mean,
			&repro_std, &log_repro_std,
			&Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newInfevents, &newVirons,
			&IthisTot,&IdthisTot,&VithisTot,&VethisTot,selectedRegions))
		    return 0;

		currentVL = 0;
		for (int j=0; j <vars->Regions; j++)
		    currentVL += Ve[j]->val();
	    }
	    else if (abs(model) == 8)
	    {
		pastVL = 0;
		for (int j=0; j <vars->Regions; j++)
		    pastVL += Ve[j]->val();

		double tempCmax;
		if (firstBolus && Cmax_0 != 0)
		    tempCmax = Cmax_0;
		else
		    tempCmax = Cmax;
		ACV += deltaACV(time_in_days, tstep, tempCmax, absorb, ACV, gamma, bolus, &infuse, &firstBolus, &lastBolus,&doses,vars->Total_doses);

		if (ACV < 0)
		    ACV = 0.;

		double factor = 1 - siliciano(vars->nT,vars->cT,vars->kD,ACV);
		if (firstPass)
		{
		    fprintf(stderr,"Siliciano factor = %lf\n",
			factor);
		    firstPass=0;
		}
		p_ACV = p * factor; 
		latent_inf_ACV = latent_inf * factor;

		if (!model5(time_in_days,Sp, Ssub, Ec,Ip,Tp,pastTps,  
			Vi, Ve, Veadj, &Iplaquebirths, &totalPlaques, plaqColors, 
			vars, Repro,logRepro,p_ACV,
			latent_inf_ACV,
			Vethis, Vithis, Ithis, Id_this, start_inf_time,
			&repro_max, &repro_mean, &log_repro_mean,
			&repro_std, &log_repro_std,
			&Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newInfevents, &newVirons,
			&IthisTot,&IdthisTot,&VithisTot,&VethisTot,selectedRegions))
		    return 0;

		currentVL = 0;
		for (int j=0; j <vars->Regions; j++)
		    currentVL += Ve[j]->val();
	    }
	    if (pastVL > 0 || currentVL > 0)
		AUC+=((((pastVL > 0)?log10(pastVL):0)+((currentVL > 0)?log10(currentVL):0))/2.0) * vars->tstep;

	    infTot=0;
	    infRegions=0;
	    vet = 0;
	    vit = 0;
	    cd8size = 0;
	    hhv8_regs = 0;
	    inf_regs = 0;
	
	    cd8s_std = getStateStddev(Tp,vars->Regions);

	    for (int j=0; j <vars->Regions; j++)
	    {
	        if (Ip[j]->val() > 0)
		    inf_regs++;
	        if (Ve[j]->val() > 0)
		    hhv8_regs++;
		if (Ithis[j]->val() > 0 && Ithis_p[j] <= 0)
		{
		    start_inf_time[j] = time_in_days;
		    start_R0[j] = Repro[j];
		}

		infTot += Ip[j]->val();
		cd8size += Tp[j]->val();
		if (abs(model) >= 4)
		{
		    vet += Ve[j]->val();
		    vit += Vi[j]->val();

		    if (Repro[j] > ReproMax[j])
			ReproMax[j] = Repro[j];
		    if (Repro[j] < ReproMin[j])
			ReproMin[j] = Repro[j];
		    if (Repro[j]<=1.0)
			repro_under1=repro_under1+1.0;

		    if (vars->density_killing) 
		    {
			ViLifeSpan = 24. / (vars->an*pow(Ip[j]->val(),vars->hill) + (vars->fpos*Tp[j]->val()));
		    } else {
			ViLifeSpan = 24. / (vars->an + (vars->fpos*Tp[j]->val()));
		    }
		    if (ViLifeMax[j] < ViLifeSpan)
			ViLifeMax[j] = ViLifeSpan;
		    if (ViLifeMin[j] > ViLifeSpan)
			ViLifeMin[j] = ViLifeSpan;

		    if ( Vi[j]->val() == 0)
		    {
			Vithis[j]->initialz_with(0); 
		    }

		    if ( Ve[j]->val() == 0)
		    {
			Vethis[j]->initialz_with(0); 
		    }

		    if (j==startRegion && Ve[j]->val() > maxStartRegionVe)
		    {
			maxStartRegionVe=Ve[j]->val();
			maxStartTime=time_in_days;
		    }

		    if (Ip[j]->val() > infLocalMax[j])
		    {
			infLocalMax[j] = Ip[j]->val();
			timeToPeak[j] = time_in_days-start_inf_time[j];
			timeAtPeak[j] = time_in_days;
			iThisAtPeak[j] = Ithis[j]->val();
		    }
		    if ( Ip[j]->val() == 0)
		    {
			Ithis[j]->initialz_with(0); 
			Id_this[j]->initialz_with(0); 
			start_inf_time[j] = 0;
		    }
		    else 
		    {
			if (Ithis_p[j] <= 0)
			{
			    start_inf_time[j] = time_in_days;
			    start_R0[j] = Repro[j];
			}
			else
			    start_inf_time[j] = 0;
			infRegions++;
			if(vars->Sig_test )
			{
			    if (Ithis[j]->val() > 0 && Ithis_p[j] <= 0)
			    {
				start_inf_time[j] = time_in_days;
				start_R0[j] = Repro[j];
			    }

			    if(Ithis[j]->val()>100000 && Ithis_p[j] <= 100000)
			    {
				if (timeTo100k[j] < 0)
				    timeTo100k[j] = time_in_days-start_inf_time[j];
				if (timeTo10k[j] < 0)
				{
				    timeTo10k[j] = time_in_days-start_inf_time[j];
				    timeTo1000[j] = time_in_days-start_inf_time[j];
				    timeTo100[j] = time_in_days-start_inf_time[j];
				}
			    }
			    else if(Ithis[j]->val()>10000 && Ithis_p[j] <= 10000)
			    {
				if (timeTo10k[j] < 0)
				    timeTo10k[j] = time_in_days-start_inf_time[j];
				if (timeTo1000[j] < 0)
				{
				    timeTo1000[j] = time_in_days-start_inf_time[j];
				    timeTo100[j] = time_in_days-start_inf_time[j];
				}
			    }
			    else if(Ithis[j]->val()>1000 && Ithis_p[j] <= 1000)
			    {
				if (timeTo1000[j] < 0)
				    timeTo1000[j] = time_in_days-start_inf_time[j];
				if (timeTo100[j] < 0)
				{
				    timeTo100[j] = time_in_days-start_inf_time[j];
				}
			    }
			    else if(Ithis[j]->val()>100 && Ithis_p[j] <= 100)
			    {
				if (timeTo100[j] < 0)
				    timeTo100[j] = time_in_days-start_inf_time[j];
			    }
			    if (Ithis[j]->val() > 0 && Ithis_p[j] <= 0)
			    {
				start_inf_time[j] = time_in_days;
			    }

			    /* track death counts for this region by episode*/
			    if(Id_this[j]->val()>100000 && Id_p[j] <= 100000)
			    {
				if (timeTo100kDeaths[j] < 0)
				    timeTo100kDeaths[j] = time_in_days-start_inf_time[j];
				if (timeTo10kDeaths[j] < 0)
				{
				    timeTo10kDeaths[j] = time_in_days-start_inf_time[j];
				    timeTo1000Deaths[j] = time_in_days-start_inf_time[j];
				    timeTo100Deaths[j] = time_in_days-start_inf_time[j];
				}
			    }
			    else if(Id_this[j]->val()>10000 && Id_p[j] <= 10000)
			    {
				if (timeTo10kDeaths[j] < 0)
				    timeTo10kDeaths[j] = time_in_days-start_inf_time[j];
				if (timeTo1000Deaths[j] < 0)
				{
				    timeTo1000Deaths[j] = time_in_days-start_inf_time[j];
				    timeTo100Deaths[j] = time_in_days-start_inf_time[j];
				}
			    }
			    else if(Id_this[j]->val()>1000 && Id_p[j] <= 1000)
			    {
				if (timeTo1000Deaths[j] < 0)
				    timeTo1000Deaths[j] = time_in_days-start_inf_time[j];
				if (timeTo100Deaths[j] < 0)
				{
				    timeTo100Deaths[j] = time_in_days-start_inf_time[j];
				}
			    }
			    else if(Id_this[j]->val()>100 && Id_p[j] <= 100)
			    {
				if (timeTo100Deaths[j] < 0)
				    timeTo100Deaths[j] = time_in_days-start_inf_time[j];
			    }
			}
		    }
		    if (Ip[j]->val() <= 0 && Ip_p[j] > 0)
		    {
			if (vars->Sig_test &&
			    vars->dataF9 != NULL && 
			    timeTo10k[j] > 0 && 
			    sig_episodes < 100)
			{
			    if (timeAtPeak[j] >= 0)
				timeFromPeak[j] = time_in_days - timeAtPeak[j] ;

			    fprintf(vars->dataF9,"%d,",sig_episodes+1);
			    fprintf(vars->dataF9,"%lf,",start_inf_time[j]);
			    fprintf(vars->dataF9,"%lf,",time_in_days);
			    fprintf(vars->dataF9,"%d,",j);
			    fprintf(vars->dataF9,"%lf,",start_R0[j]);
			    fprintf(vars->dataF9,"%lf,",timeTo100[j]);
			    fprintf(vars->dataF9,"%lf,",timeTo1000[j]);

			    fprintf(vars->dataF9,"%lf,",timeTo10k[j]);
			    if (timeTo100k[j] > 0)
				fprintf(vars->dataF9,"%lf,",timeTo100k[j]);
			    else
				fprintf(vars->dataF9,"*,");
			    fprintf(vars->dataF9,"%lf,",timeTo100Deaths[j]);
			    fprintf(vars->dataF9,"%lf,",timeTo1000Deaths[j]);

			    fprintf(vars->dataF9,"%lf,",timeTo10kDeaths[j]);

			    if (timeTo100kDeaths[j] > 0)
				fprintf(vars->dataF9,"%lf,",timeTo100kDeaths[j]);
			    else
				fprintf(vars->dataF9,"*,");

			    fprintf(vars->dataF9,"%lf,",timeToPeak[j]);
			    fprintf(vars->dataF9,"%lf,",timeFromPeak[j]);
			    fprintf(vars->dataF9,"%lf,",time_in_days-start_inf_time[j]);
			    fprintf(vars->dataF9,"%lu,",infLocalMax[j]);
			    fprintf(vars->dataF9,"%lu,", infGlobalMax);
			    fprintf(vars->dataF9,"%lf,", timeAtGlobalPeak);
			    fprintf(vars->dataF9,"%lu\n", iThisAtPeak[j]);

			    sig_episodes++;
			    if (sig_episodes == 100)
			    {
				    episodeStop = 1;
			    }
			}
			timeTo100[j]=-1;
			timeTo1000[j]=-1;
			timeTo10k[j]=-1;
			timeTo100k[j]=-1;
			timeTo100Deaths[j]=-1;
			timeTo1000Deaths[j]=-1;
			timeTo10kDeaths[j]=-1;
			timeTo100kDeaths[j]=-1;
			timeToPeak[j]=-1;
			timeAtPeak[j]=-1;
			iThisAtPeak[j]=0;
			infLocalMax[j]=0;
			start_inf_time[j]=0;
			start_R0[j]=0;
		    }
		    Ip_p[j] = Ip[j]->val();
		    Ithis_p[j] = Ithis[j]->val();
		    Id_p[j] = Id_this[j]->val();
		}
	    }
	    repro_under1=repro_under1/vars->Regions;
	    running_repro_under1 = ((running_repro_under1 * steps)+repro_under1)/(steps+1);
	    running_cd8s_max = max(running_cd8s_max,(double)cd8size);
	    running_cd8s_mean = ((running_cd8s_mean * steps)+cd8size)/(steps+1);
	    running_cd8s_std = ((running_cd8s_std * steps)+cd8s_std)/(steps+1);
	    running_repro_max = max(running_repro_max, repro_max);
	    running_repro_mean = ((running_repro_mean * steps)+repro_mean)/(steps+1);
	    running_repro_std = ((running_repro_std * steps)+repro_std)/(steps+1);

	    running_log_repro_mean = ((running_log_repro_mean * steps)+log_repro_mean)/(steps+1);
	    running_log_repro_std = ((running_log_repro_std * steps)+log_repro_std)/(steps+1);

	    avg_hhv8_regs = ((avg_hhv8_regs * steps)+hhv8_regs)/(steps+1);
	    avg_inf_regs = ((avg_inf_regs * steps)+inf_regs)/(steps+1);

	    if (Iplaquebirths > 0 && Itot_p > 0)
		extra_plaque_starts++;
	    if (Iplaquebirths > 0)
		plaque_births++;

	    steps++;
	
	    /* end of infected cell episode -> zero counters! */
	    if (infTot > infGlobalMax)
	    {
		infGlobalMax = infTot;
		timeAtGlobalPeak = time_in_days;
	    }

	    if (vet > 0 && vet_p <= 0)
	    {
		double err=0;
		double global_mean=(double)cd8size/vars->Regions;
		if (vars->Verbose > 1)
		    fprintf(results,"At infection start (t=%lf), vet=(%lu), Tcell total=%lu\n",time_in_days,vet,cd8size);
		for (int j=0; j <vars->Regions; j++)
		{
		    err += pow((global_mean-Tp[j]->val()),2.0);
		    if(Ve[j]->val() > 0)
		    {
			if (vars->Verbose > 1)
			    fprintf(results,"At infection site (cell %d), Tcell count=%lu\n",j,Tp[j]->val());
			double Tmean=Tp[j]->val();
			int cnt=1;
			for (int k =0; k < MAX_NEIGHBORS; k++)
			    if (vars->cells[j]->neighbors[k] != NULL)
			    {
				Tmean += Tp[vars->cells[j]->neighbors[k]->cell]->val();
				cnt++;
			    }

			Tmean = Tmean / cnt;
			if (vars->Verbose > 1)
			    fprintf(results,"In infection neighborhood, Tcell mean=%lf\n",Tmean);
			startRegion=j;
			maxStartRegionVe=0;
			break;
		    }
		}
		//fprintf(results,"At infection time, Tcell stddev=%lf\n",sqrt(err/(vars->Regions-1)));
	    }
	    vet_p=vet;

	    if (vet <= 0 )
	    {
		VethisTot = 0;
	    }

	    if (vit <= 0 )
		VithisTot = 0;

	    if (infTot <= 0 && Itot_p > 0)
	    {
		if (infGlobalMax > 10000)
		{
		    fprintf(results,"global peak of %lu at %lf\n",infGlobalMax,timeAtGlobalPeak);
		}
		infGlobalMax=0;
	    }
	    if (infTot <= 0)
	    {
		IthisTot = 0;
		IdthisTot = 0;
		Itot_p = 0;
	    }
	    else
		Itot_p=infTot;

	    Ve_to_date += newVirons;
	    Inf_to_date += newInfcells;
	    Spread_to_date += newInfevents;
	    T_to_date += newTcells;
	    T_inf_tot += Tinf;
	    T_bump_tot += Tbump;
	    T_dump_tot += Tdump;
	    T_kills_tot += Tkills;
	    Ideaths_tot += Ideaths;

	    /* swabs done at specified intervals (6hr or 1 day) */
	    // if Calc_T0 ==2 the only swab AFTER the "warmup" period
	    // if Calc_T0 ==3 swab during the "warmup" period (can be in episode)
	    if ((time_in_days == 0 || vars->swabInterval == 0 || time_in_days - swabT >= vars->swabInterval) &&
		   (vars->Calc_T0 != 2 || time_in_days > T0_Sim))
	    {
		measuredVL = currentVL;
		if (vars->CritOn && time_in_days >= vars->crit_start)
		    thisSubjSwabs++;

		if (time_in_days >= vars->crit_start && measuredVL<= shed_thresh )
		    any_neg=true;

		//starting an episode?
		if (vars->CritOn && time_in_days >= vars->crit_start &&
		    !inEpisode && measuredVL>shed_thresh /*&& infRegions >= vars->infThreshold*/)
		{
		    fprintf(results, 
			"Episode start at %lf days\n", time_in_days);
		    if (vars->Verbose || vars->Size_limit > 0)
			fprintf(results, 
			    "Episode start at %lf days, ACV=%lf elapsed=%lf cmax=%lf absorb=%lf gamma=%lf use_rho=%d\n", 
			    time_in_days,ACV,time_in_days-lastBolus,
			    Cmax,absorb,gamma,vars->use_rho);

		    for (int j=0; j <vars->Regions; j++)
		    {
			pre_Tp[j]->initialz_with(Tp[j]->val());
		    }
		    //first of episode
		    First_T = time_in_days;

		    inEpisode=true;

		    for (int j=0; j <vars->Regions; j++)
		    {
			if (abs(model) >= 4 && Ve[j]->val() > 0)
			{
			    firstRegion = j;
			    break;
			}
		    }
		}

		//ending an episode? 
		if (inEpisode && measuredVL<= shed_thresh )
		{
		    maxVL_bin = MIN(EPI_BINS-1,MAX(0,(int)(log10(maxVL[counter])-log10(shed_thresh)))); 
		    peaks[maxVL_bin]++;
		    totPeaks[totalEpis+counter] = log10(maxVL[counter]);

		    fprintf(results,"Start region %d has max of %lu virons(at t=%lf)\n",startRegion,maxStartRegionVe,maxStartTime);
		    startRegion=-1;
		    maxStartRegionVe=0;

		    // for durations, track days using first and last
		    // positive sample times plus random portions of the
		    // sampling interval before and after those sample times
		    if (First_T == Last_T)
			episodDur[counter]=1;
		    else
			episodDur[counter]=Last_T - First_T; // + vars->swabInterval*(gsl_rng_uniform(vars->ur)+gsl_rng_uniform(vars->ur));

		    episodDur[counter]=swabsThisEpi;
		    totDurations[totalEpis+counter] = episodDur[counter];

		    //if (vars->Verbose || vars->Size_limit > 0)
			fprintf(results,
			    "Episode end at %lf days: Ve peak=%lu (bin=%d), Episode duration=%lf days (peak at t=%lf)\n",
			    time_in_days,maxVL[counter], maxVL_bin, episodDur[counter], max_T);



		    if (episodDur[counter] > max_dur) max_dur = episodDur[counter];

		    totPlaques[counter] = totalPlaques;
		    totalPlaques = 0;
		    firstRegion = -1;

		    for (int j=0; j <vars->Regions; j++)
		    {
			if((double)Tp[j]->val()-(double)pre_Tp[j]->val() > 50)
			    cd8_reexpansions++;
		    }
		    if (vars->writeOn && abs(model) >= 4 &&
		       (vars->Calc_T0 != 2 || time_in_days > T0_Sim))
		    {
			if (vars->dataF10 != NULL) {
			    fprintf(vars->dataF10,"%d,",counter);
			    fprintf(vars->dataF10,"%lf,",episodDur[counter]);
			    fprintf(vars->dataF10,"%lu,",maxVL[counter]);
			    fprintf(vars->dataF10,"%lf,",totPlaques[counter]);
			    fprintf(vars->dataF10,"%lf,",maxFirstReg[counter]);
			    fprintf(vars->dataF10,"\n");
			}
		    }

		    /* bump episode count at end of episode */
		    counter++;
		    inEpisode=false;
		    swabsThisEpi=0.;
		    avgVL=0.;

		    pcounter++;

		    if (pcounter == 100) 
		    {
			long time_now;
			time(&time_now);
			if (vars->Verbose)
			    fprintf(stderr,"(Patient %d hit %d episodes at t=%lf (elapsed=%ld)\n",
				Nr_count,counter,time_in_days,time_now-vars->time_st);
			pcounter = 0;
		    }
		}

		// in the midst of an episode?
		if (inEpisode)
		{
		    if(measuredVL>maxVL[counter])
		    {
			maxVL[counter] = measuredVL;
			max_T = time_in_days;
		    }
		    if (abs(model) >= 4 && firstRegion >= 0 && Ve[firstRegion]->val() > maxFirstReg[counter])
		    {
			maxFirstReg[counter] = Ve[firstRegion]->val();
		    }
		
		    /* calc a running average VL per episode for criteria10 */
		    avgVL = ((avgVL * swabsThisEpi) + measuredVL) / (swabsThisEpi + 1);
		    swabsThisEpi++;

		    //VL_bin = MIN(VL_BINS-1,MAX(0,(int)((log10(measuredVL)-log10(shed_thresh))+0.5))); 
		    VL_bin = MIN(VL_BINS-1,MAX(0,(int)(1.0 + log10(measuredVL)-log10(shed_thresh)))); 

		    swabs[VL_bin]++;

		    thisPosSwabs[posSwabs]=log10(measuredVL);
		    if (posSwabs < MAX_SWABS)
			posSwabs++;

		    totPosSwabs[totalPosSwabs]=log10(measuredVL);
		    if (totalPosSwabs < MAX_SWABS)
			totalPosSwabs++;
		    Last_T = time_in_days;

		    //if (First_T == Last_T)
		//	episodDur[counter]=1;
		 //   else
		//	episodDur[counter]=Last_T - First_T; // + vars->swabInterval*(gsl_rng_uniform(vars->ur)+gsl_rng_uniform(vars->ur));
		    episodDur[counter]=swabsThisEpi;

		    if (episodDur[counter] > max_dur) max_dur = episodDur[counter];

		}
		else if (vars->CritOn && time_in_days >= vars->crit_start)
		{
		    swabs[0]++;
		}
		if (vars->CritOn && time_in_days >= vars->crit_start)
		{
		    if (totalSwabs < MAX_SWABS)
			totalSwabs++;

		    swabT=time_in_days;
		}
	    }
	    /* catch non-swabbed episodes for statistical purposes */
	    if(currentVL>shed_thresh /*&& infRegions >= vars->infThreshold */&& 
		!inContEpisode)
	    {
		cont_First_T=time_in_days;
		inContEpisode=true;
		AUC=0;
	    }
	    if (inContEpisode && currentVL<=0)
	    {
		period_episodDur[period_counter]=time_in_days - cont_First_T;

		cont_counter++;
		period_counter++;

		total_AUC += AUC;
		period_AUC += AUC;
		AUC=0;
		inContEpisode=false;
	    }
	    else if (inContEpisode)
		pos_this_period += tstep;

	    pastVL = currentVL;

	    if (time_in_days - period_p > vars->statInterval) 
	    {
		double Mean_duration = getMean(period_episodDur, period_counter);  
		fprintf(results,
		    "Period Shedding rate = %lf%%\n",
		    pos_this_period/vars->statInterval);
		fprintf(results,
		    "Episode rate this period = %lf\n",
		    period_counter*(365./vars->statInterval));
		fprintf(results,
		    "Period Mean episode duration = %lf days\n",
		    Mean_duration);

		period_counter=0;
		period_AUC=0;
		period_p=time_in_days;
		pos_this_period = 0;
	    }

	    // use sampling interval for graphical and file output frequency
	    if (time_in_days-timep > vars->sampling) { 

		if (vars->writeOn && abs(model) >= 4 &&
		   (vars->Calc_T0 != 2 || time_in_days > T0_Sim))
		{
		    if (vars->dataF1 != NULL) {
			fprintf(vars->dataF1,"%lf,",time_in_days);
			fprintf(vars->dataF1,"%lu,",vet);
			for (int j=0; j <vars->Regions; j++)
			    fprintf(vars->dataF1,"%lu,",Ve[j]->val());
			fprintf(vars->dataF1,"\n");
		    }
		    if (vars->dataF2 != NULL) {
			fprintf(vars->dataF2,"%lf,",time_in_days);
			fprintf(vars->dataF2,"%lu,",VethisTot);
			for (int j=0; j <vars->Regions; j++)
			    fprintf(vars->dataF2,"%lu,",Vethis[j]->val());
			fprintf(vars->dataF2,"\n");
		    }
		    if (vars->dataF3 != NULL) {
			fprintf(vars->dataF3,"%lf,",time_in_days);
			fprintf(vars->dataF3,"%lu,",vit);
			for (int j=0; j <vars->Regions; j++)
			    fprintf(vars->dataF3,"%lu,",Vi[j]->val());
			fprintf(vars->dataF3,"\n");
		    }
		    if (vars->dataF4 != NULL) {
			fprintf(vars->dataF4,"%lf,",time_in_days);
			fprintf(vars->dataF4,"%lu,",VithisTot);
			for (int j=0; j <vars->Regions; j++)
			    fprintf(vars->dataF4,"%lu,",Vithis[j]->val());
			fprintf(vars->dataF4,"\n");
		    }
		    if (vars->dataF5 != NULL) {
			fprintf(vars->dataF5,"%lf,",time_in_days);
			fprintf(vars->dataF5,"%lu,",infTot);
			for (int j=0; j <vars->Regions; j++)
			    fprintf(vars->dataF5,"%lu,",Ip[j]->val());
			fprintf(vars->dataF5,"\n");
		    }
		    if (vars->dataF6 != NULL) {
			fprintf(vars->dataF6,"%lf,",time_in_days);
			fprintf(vars->dataF6,"%lu,",IthisTot);
			fprintf(vars->dataF6,"%lu,",IdthisTot);
			fprintf(vars->dataF6,"\n");
		    }
		    if (vars->dataF8 != NULL) {
			fprintf(vars->dataF8,"%lf,",time_in_days);
			for (int j=0; j <vars->Regions; j++)
			    fprintf(vars->dataF8,"%lf,",Repro[j]);
			fprintf(vars->dataF8,"\n");
		    }
		    if (vars->dataF11 != NULL) {
			unsigned long int tmax = 0;
			unsigned long int ttot = 0;
			unsigned long int vemax = 0;
			unsigned long int vetot = 0;
			fprintf(vars->dataF11,"%lf,",time_in_days);
			for (int j=0; j <vars->Regions; j++)
			{
			    if (Tp[j]->val() > tmax)
				tmax = Tp[j]->val();
			    ttot +=Tp[j]->val();
			    if (Ve[j]->val() > vemax)
				vemax = Ve[j]->val();
			    vetot +=Ve[j]->val();
			}
			fprintf(vars->dataF11,"%lu",tmax);
			fprintf(vars->dataF11,",%lu",ttot);
			fprintf(vars->dataF11,",%lu",vemax);
			fprintf(vars->dataF11,",%lu",vetot);
			fprintf(vars->dataF11,",%lf",repro_max);
			fprintf(vars->dataF11,",%lf",repro_mean);
			fprintf(vars->dataF11,",%lf",repro_std);
			fprintf(vars->dataF11,"\n");
		    }
		    if (vars->dataF12 != NULL) 
		    {
			fprintf(vars->dataF12,"%lf,",time_in_days);
			fprintf(vars->dataF12,"%d,",model);
			fprintf(vars->dataF12,"%lf,",bolus);
			fprintf(vars->dataF12,"%d,",doses);
			fprintf(vars->dataF12,"%lf,",ACV);
			fprintf(vars->dataF12,"%lf,",infuse);
			fprintf(vars->dataF12,"%lf,",lastBolus);
			fprintf(vars->dataF12,"%lf,",p);
			fprintf(vars->dataF12,"%lf,",p_ACV);
			fprintf(vars->dataF12,"%lf,",latent_inf);
			fprintf(vars->dataF12,"%lf,",latent_inf_ACV);
			fprintf(vars->dataF12,"%lu",vet);
			fprintf(vars->dataF12,"\n");
		    }
		}

		vars->time = time_in_days;

	      timep = time_in_days;
	  }

	  time_in_days = time_in_days + tstep;

	  if (vars->CritOn && time_in_days >= vars->crit_start &&
	      (vars->Calc_T0 != 2 || time_in_days > T0_Sim))
	  {
	      total_stats_time += tstep;
	      stats_time += tstep;
	  }

	  if (time_in_days>100.0  &&
	      (time_in_days  - (100.0 * (int)(time_in_days/100.0)) <= tstep))
	  { 
		T0 = 0;
		for (int j=0; j <vars->Regions; j++)
		{
		    T0 += Tp[j]->val();
		}
		if (vars->Verbose > 1)
		    fprintf(stderr, 
			"Run %d at %lf days (%g %% positive swabs, %lu tcells)...\n", 
			Nr_count, time_in_days,
			100*(float)totalPosSwabs/(float)totalSwabs, T0);

	  }


	}//end of time simulation
	if (inContEpisode )
	{
	    period_episodDur[period_counter]=time_in_days - cont_First_T;

	    cont_counter++;
	    period_counter++;

	    total_AUC += AUC;
	    period_AUC += AUC;
	    AUC=0;
	    inContEpisode=false;
	}
	if (inEpisode)
	{
	    Last_T = time_in_days;
	    episodDur[counter]=swabsThisEpi;
	    if (episodDur[counter] > max_dur) max_dur = episodDur[counter];

	    totPlaques[counter] = totalPlaques;

	    if (episodDur[counter] > max_dur) max_dur = episodDur[counter];

	    maxVL_bin = MIN(EPI_BINS-1,MAX(0,(int)(log10(maxVL[counter])-log10(shed_thresh)))); 
	    peaks[maxVL_bin]++;
	    totPeaks[totalEpis+counter] = log10(maxVL[counter]);
	    totDurations[totalEpis+counter] = episodDur[counter];
	    fprintf(results,
		"Episode end at %lf days: Ve peak=%lu (bin=%d), Episode duration=%lf days (peak at t=%lf)\n",
		time_in_days,maxVL[counter], maxVL_bin, episodDur[counter], max_T);

	    counter++;
	}
	fprintf(results, "Run %d completed in %lf days with %g positive swabs (max episode = %lf days, total log AUC=%lf)\n", 
		Nr_count, time_in_days, totalPosSwabs/(float)totalSwabs,max_dur,total_AUC);
	fprintf(results,
		"beta=%e latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf cd8_ic50=%lf kappa=%lf alpha=%lf exp_days=%lf hill=%lf betae=%e\n",
		vars->beta_init,vars->latent_inf_init,vars->log_p_init,vars->c_init,vars->theta_init,vars->delta_init,vars->r_init,vars->inf_init,vars->rinf_init,vars->rho_init,vars->eclipse_init,vars->fpos,vars->an,vars->cd8_ic50_init,vars->kappa_init,vars->alpha_init,vars->exp_days_init,vars->hill,vars->betae_init);

	fprintf(results,"Total infect cell ripenings %lu\n",Inf_to_date);
	fprintf(results,"Total infect cell spread events %lu\n",Spread_to_date);
	fprintf(results,"Total infect cell deaths (nat) %lu\n",Ideaths_tot);
	fprintf(results,"Total virons produced %lu\n",Ve_to_date);
	fprintf(results,"Total T cells infused %lu\n",T_inf_tot);
	fprintf(results,"Total T cells produced %lu\n",T_to_date);
	fprintf(results,"Total T cells kills %lu\n",T_kills_tot);
	fprintf(results,"Total T cells added due to rho %lu\n",T_bump_tot);
	fprintf(results,"Total T cells removed due to rho %lu\n",T_dump_tot);

	// if no episode end, then set peak here (missed otherwise!)
	//if (!any_neg)
	//{
	 //   maxVL_bin = MIN(EPI_BINS-1,MAX(0,(int)(log10(maxVL[counter])-log10(shed_thresh)))); 
	  //  peaks[maxVL_bin]++;
	//}

	

	//VL summary
	for(int j=0;j< VL_BINS;j++)
		criteria1[j] += swabs[j];
	
	//Peak summary
	for(int j=0;j< EPI_BINS;j++)
		criteria2[j] += peaks[j];
	
	//Longest episode duration summary
	// no episodes
	if (!any_neg)
	    criteria3[5]++;
	else if (counter == 0)
	    criteria3[0]++;
	else if (max_dur < 2)
	    criteria3[1]++;
	else if (max_dur >= 2 && max_dur < 6)
	    criteria3[2]++;
	else if (max_dur >= 6 && max_dur <= 10)
	    criteria3[3]++;
	else 
	    criteria3[4]++;
	
	totalEpis += counter;
	cont_totalEpis += cont_counter;

	fprintf(results,"%d episodes (tot = %d)for patient %d (cont epis=%d)\n",
	    counter,totalEpis,Nr_count,cont_totalEpis);

	time(&time_now);
	
	if (posSwabs > 0)
	{
	    double Med_swabVL = getMedian(thisPosSwabs, posSwabs);  
	    double Max_swabVL = getMax(thisPosSwabs, posSwabs);  
	    double S_swabVL = getStddev(thisPosSwabs, posSwabs);  
	    patSwabMeds[withPosSwabs] = Med_swabVL;
	    patSwabVars[withPosSwabs] = S_swabVL*S_swabVL;
	    shedderSwabs+=thisSubjSwabs;
	    withPosSwabs++;
	    fprintf(results,"Subject Percent Pos swabs = %lf (%d pos swabs of %d)\n",
		100.0*posSwabs/thisSubjSwabs,posSwabs,thisSubjSwabs);
	    fprintf(results,"Subject Median log swab VL = %lf\n",Med_swabVL);
	    fprintf(results,"Subject Peak log swab VL = %lf\n",Max_swabVL);
	    fprintf(results,"Subject Var log swab VL = %lf (%d pos swabs of %d)\n",S_swabVL*S_swabVL,posSwabs,thisSubjSwabs);
	} else {
	    fprintf(results,"Subject %d had no positive swabs\n",Nr_count);
	}
	posSwabs=0;
	thisSubjSwabs=0;

	fprintf(results,"Max CD8s = %lf\n",running_cd8s_max);
	fprintf(results,"Mean CD8s = %lf\n",running_cd8s_mean);
	fprintf(results,"Mean Stddev CD8s = %lf\n",running_cd8s_std);
	fprintf(results,"Mean percent repro < 1.0 = %lf\n",running_repro_under1);
	fprintf(results,"Max repro = %lf\n",running_repro_max);
	fprintf(results,"Mean repro = %lf\n",running_repro_mean);
	fprintf(results,"Mean log repro = %lf\n",running_log_repro_mean);
	fprintf(results,"Mean repro stddev= %lf\n",running_repro_std);
	fprintf(results,"Mean log repro stddev= %lf\n",running_log_repro_std);
	fprintf(results,"Mean regions with HHV8 (Ve>0) = %lf\n", avg_hhv8_regs);
	fprintf(results,"Mean infected regions = %lf\n", avg_inf_regs);
	fprintf(results,"Starts per plaque = %lf\n", (plaque_births-extra_plaque_starts<=0)?0:(double)(plaque_births)/(plaque_births-extra_plaque_starts));
	fprintf(results,"CD8 reexpansions per year = %lf\n", cd8_reexpansions*(365./total_stats_time));

	Nr_count++;
    }//end of runs loop

    if (totalEpis == 0)
    {
	fprintf(results,"No episodes!\n");
	if (vars->Model > 5)
	    fprintf(results,
		"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf ,absorb=%lf, gamma=%lf,Cmax=%lf,IC50=%lf,m=%lf,time=%lf\n",
		beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,absorb,gamma,Cmax,IC50,m,stats_time);
	else
	    fprintf(results,
		"To=%lf,beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf ,time=%lf\n",
		vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,stats_time);
	score = 0.;
    }
    else
    {
	fprintf(results,"%d episodes for %d patients (%d with pos swabs)\n", totalEpis,N_runs,withPosSwabs);

	double Mean_swabVL = getMean(totPosSwabs, totalPosSwabs);  
	double Mean_med_swabVL = getMean(patSwabMeds, withPosSwabs);  
	double Mean_var_swabVL = getMean(patSwabVars, withPosSwabs);  

	double Mean_peakVL = getMean(totPeaks, totalEpis);  
	double S_peakVL = getStddev(totPeaks, totalEpis);  

	double Mean_duration = getMean(totDurations, totalEpis);  
	qsort (totDurations, totalEpis, sizeof(double), compare_doubles);
	double Med_duration = getMedian(totDurations, totalEpis);  
	double S_duration = getStddev(totDurations, totalEpis);  

	fprintf(results,"Mean log swab VL = %lf\n",Mean_swabVL);
	fprintf(results,"Median log swab VL = %lf\n",Mean_med_swabVL);
	fprintf(results,"Var log swab VL = %lfs\n",Mean_var_swabVL);
	fprintf(results,"Mean log peak VL = %lf\n",Mean_peakVL);
	fprintf(results,"Var peak VL = %lf\n",S_peakVL*S_peakVL);
	fprintf(results,"Continuous episodes per year = %lf\n",cont_totalEpis*(365./total_stats_time));
	fprintf(results,"Mean episode duration = %lf days\n",Mean_duration);
	fprintf(results,"Median episode duration = %lf days\n",Med_duration);
	fprintf(results,"Var episode duration = %lf days\n",S_duration*S_duration);

    }

    // Criteria 1: Quantitative shedding frequency
    // proportion of + swabs at each strata.
    // < 10, < 10^3, < 10^4, < 10^5, < 10^6, < 10^7, < 10^8, < 10^9, >= 10^9.

    int critNum=0;

    double mean_mean_1; /* avg of means in category 1 */
    double mean_mean_2; /* avg of means in category 2 */
    double mean_mean_3; /* avg of means in category 3 */
    
    double totCrit1Perc=0.;
    double subCrit1Perc=0.; // only for runs with at least 1 positve swab!
    
    //VL summary
    for(int j=0;j< VL_BINS;j++)
	crit1perc[j] = 100. * criteria1[j]/totalSwabs;

    for(int j=1;j< VL_BINS;j++)
	subCrit1Perc += 100. * criteria1[j]/shedderSwabs;

    // make percentages cumulative to match CIs
    if (vars->Match_strategy == 0)
    {
	for(int j=1;j< VL_BINS;j++)
	    crit1perc[j] += crit1perc[j-1];
    }
    else 
    {
	mean_mean_1 = 0;
	/* skip 1st bin for category 1 */
	for(int j=1;j< VL_BINS;j++)
	    mean_mean_1 += vars->crit[critNum+j].mean;

	mean_mean_1 = mean_mean_1/(VL_BINS-1);
    }
    for(int j=1;j< VL_BINS;j++)
	totCrit1Perc += crit1perc[j];

    fprintf(results,"Total percentage of swabs that were positive= %lf (%d of %d)\n",totCrit1Perc,totalPosSwabs,totalSwabs);
    fprintf(results,"Total percentage of shedder swabs that were positive= %lf (%d of %d)\n",subCrit1Perc,totalPosSwabs,shedderSwabs);
    fprintf(results,"Rate of episodes with >1mm plaques= %lf\n",epi_1mm* (365./total_stats_time));
    fprintf(results,"Percent Time in episodes with >1mm plaques= %lf\n",100.*time_over_1/total_stats_time);
    fprintf(results,"Rate of episodes with >2mm plaques= %lf\n",epi_2mm* (365./total_stats_time));
    fprintf(results,"Percent Time in episodes with >2mm plaques= %lf\n",100.*time_over_2/total_stats_time);

    /* skip 1st measurement in both targets and actuals*/
    critNum=1;
    for(int j=1;j< VL_BINS;j++)
    {
	if (vars->Match_strategy > 0)
	{
	    double err;
	    // Match_strategy=1 -> match mean and 
	    // combine 1st two bins for criteria 1 & 2
	    if (vars->Match_strategy == 1 && j < 2)
	    {
		err = abs((vars->crit[1].mean - crit1perc[1])+
			    (vars->crit[2].mean - crit1perc[2]));
		err = err/2.0;
	    }
	    else
		err = abs((vars->crit[critNum].mean - crit1perc[j]));

	    double weight;

	    // 1st bin gets 1.0 weighting!
	    weight = vars->critWeight[0]/((VL_BINS-1)*mean_mean_1);
	    fprintf(results,"criteria 1[%d]: %lf vs. mean of %lf - err = %lf * %lf * %lf\n",
		j+1,crit1perc[j], vars->crit[critNum].mean,err,weight,vars->critWeight[0]);
	    score1+= err*weight;
	}
	else 
	{
	    fprintf(results,"criteria 1[%d]: %lf < %lf < %lf - ",
		j+1,vars->crit[critNum].low, crit1perc[j], vars->crit[critNum].high);
	    
	    if (crit1perc[j] < vars->crit[critNum].low )
	    {
		fprintf(results,"no\n");
		score1-= (vars->crit[critNum].low - crit1perc[j])/100.; 
	    }
	    else if (crit1perc[j] > vars->crit[critNum].high)
	    {
		fprintf(results,"no\n");
		score1-= (crit1perc[j] - vars->crit[critNum].high)/100.; 
	    }
	    else
	    {
		fprintf(results,"yes\n");
		score1+=1.0; 
	    }
	}
	critNum++;
    }
    fprintf(results,"Score1 = %lf%s\n",score1,(vars->Crit_mask == 0 || (vars->Crit_mask & 1))?"*":"");

    if (vars->Crit_mask == 0 || (vars->Crit_mask & 1)) score += score1;

    // Criteria 2: Peak copy # (8)
    // percentage with peak in each of given strata
    // < 10^3, < 10^4, < 10^5, < 10^6, < 10^7, < 10^8, < 10^9, >= 10^9.
    if (totalEpis == 0)
    {
	crit2perc[0] = 100.;
	for(int j=1;j< EPI_BINS;j++)
	    crit2perc[j] = 0.;
    }
    else
    {
	for(int j=0;j< EPI_BINS;j++)
	    crit2perc[j] = 100. * criteria2[j]/totalEpis;
    }
    //
    // make percentages cumulative to match CIs
    if (vars->Match_strategy == 0)
    {
	for(int j=1;j< EPI_BINS;j++)
	    crit2perc[j] += crit2perc[j-1];
    }
    else 
    {
	mean_mean_2 = 0;
	for(int j=0;j< EPI_BINS;j++)
	    mean_mean_2 += vars->crit[critNum+j].mean;

	mean_mean_2 = mean_mean_2/EPI_BINS;
    }

    int crit2start = critNum;
    for(int j=0;j< EPI_BINS;j++)
    {
	if (vars->Match_strategy > 0)
	{
	    double err;
	    double weight;
	    // combine 1st two bins for criteria 1 & 2
	    if (vars->Match_strategy == 1 && j < 2)
	    {
		err = abs((vars->crit[crit2start].mean + vars->crit[crit2start+1].mean) - (crit2perc[0]+crit2perc[1]));
		err = err/2.0;
	    }
	    // combine every two bins for criteria 2
	    else if (vars->Match_strategy == 2 && j < EPI_BINS -1 && j % 2 == 0)
	    {
		err = abs((vars->crit[critNum].mean + vars->crit[critNum+1].mean) - (crit2perc[j]+crit2perc[j+1]));
		err = err/2.0;
		fprintf(results,"criteria 2[%d (%d-%d)]: %lf vs. mean of %lf - err = %lf * %lf\n",
		    j+1,j+1,j+2,crit2perc[j], vars->crit[critNum].mean,err,vars->critWeight[1]);
	    }
	    else if (vars->Match_strategy == 2 && j % 2 == 1)
	    {
		err = abs((vars->crit[critNum].mean + vars->crit[critNum-1].mean) - (crit2perc[j]+crit2perc[j-1]));
		err = err/2.0;
		fprintf(results,"criteria 2[%d (%d-%d)]: %lf vs. mean of %lf - err = %lf * %lf\n",
		    j+1,j,j+1,crit2perc[j], vars->crit[critNum].mean,err,vars->critWeight[1]);
	    }
	    else
	    {
		err = abs((vars->crit[critNum].mean - crit2perc[j]));
		fprintf(results,"criteria 2[%d]: %lf vs. mean of %lf - err = %lf * %lf\n",
		    j+1,crit2perc[j], vars->crit[critNum].mean,err,vars->critWeight[1]);
	    }

	    weight = (vars->critWeight[1]/mean_mean_2) / EPI_BINS;
	    score2+= weight *  err;
	}
	else
	{
	    fprintf(results,"criteria 2[%d]: %lf < %lf < %lf - ",
		j+1,vars->crit[critNum].low, crit2perc[j], vars->crit[critNum].high);
	    
	    if (crit2perc[j] < vars->crit[critNum].low )
	    {
		fprintf(results,"no\n");
		score2-= (vars->crit[critNum].low - crit2perc[j])/100.; 
	    }
	    else if (crit2perc[j] > vars->crit[critNum].high)
	    {
		fprintf(results,"no\n");
		score2-= (crit2perc[j] - vars->crit[critNum].high)/100.; 
	    }
	    else
	    {
		fprintf(results,"yes\n");
		score2+=1.0; 
	    }
	}
	critNum++;
    }
    fprintf(results,"Score2 = %lf%s\n",score2,(vars->Crit_mask == 0 || (vars->Crit_mask & 0x2))?"*":"");

    // Criteria 3: Duration category # (6)
    // percentage with longest duration in each of given strata
    // 0 days, 1 day, 2-5 days, 6-10 days, >10 days w/ negs, all days
    if (totalEpis == 0)
    {
	crit3perc[0] = 100.;
	for(int j=1;j< DUR_BINS;j++)
	    crit3perc[j] = 0.;
    }
    else
    {
	for(int j=0;j< DUR_BINS;j++)
	    crit3perc[j] = 100. * criteria3[j]/N_runs;
    }
    //
    // make percentages cumulative to match CIs
    if (vars->Match_strategy == 0)
    {
	for(int j=1;j< DUR_BINS;j++)
	    crit3perc[j] += crit3perc[j-1];
    }
    else 
    {
	mean_mean_3 = 0;
	for(int j=0;j< DUR_BINS;j++)
	    mean_mean_3 += vars->crit[critNum+j].mean;

	mean_mean_3 = mean_mean_3/DUR_BINS;
    }

    for(int j=0;j< DUR_BINS;j++)
    {
	if (vars->Match_strategy > 0)
	{
	    double err;
	    double weight;
	    err = abs((vars->crit[critNum].mean - crit3perc[j]));
	    fprintf(results,"criteria 3[%d]: %lf vs. mean of %lf - err = %lf * %lf\n",
		j+1,crit3perc[j], vars->crit[critNum].mean,err,vars->critWeight[1]);

	    weight = (vars->critWeight[1]/mean_mean_3) / DUR_BINS;
	    score3+= weight *  err;
	}
	else
	{
	    fprintf(results,"criteria 3[%d]: %lf < %lf < %lf - ",
		j+1,vars->crit[critNum].low, crit3perc[j], vars->crit[critNum].high);
	    
	    if (crit3perc[j] < vars->crit[critNum].low )
	    {
		fprintf(results,"no\n");
		score3-= (vars->crit[critNum].low - crit3perc[j])/100.; 
	    }
	    else if (crit3perc[j] > vars->crit[critNum].high)
	    {
		fprintf(results,"no\n");
		score3-= (crit3perc[j] - vars->crit[critNum].high)/100.; 
	    }
	    else
	    {
		fprintf(results,"yes\n");
		score3+=1.0; 
	    }
	}
	critNum++;
    }
    fprintf(results,"Score3 = %lf%s\n",score3,(vars->Crit_mask == 0 || (vars->Crit_mask & 0x4))?"*":"");

    if (vars->Crit_mask == 0 || (vars->Crit_mask & 0x4)) score += score3;

    if (critNum != NUM_CRITERIA) 
    {
	fprintf(results,"Internal error: processed %d criteria out of %d in file\n",
	    critNum,NUM_CRITERIA);
	exit(1);
    }
    if (vars->Match_strategy != 0)
	score = 1.0 / score; /* invert so that low actual score (least squares) is desireable (high)! */

    if (results != stdout)
	fprintf(stdout,"Total Score = %lf ",score);

    fprintf(results,"Total Score = %lf ",score);
    if (vars->Model > 5)
    {
	if (results != stdout)
	    fprintf(stdout,
	    "beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf absorb=%lf gamma=%lf Cmax=%lf IC50=%lf m=%lf time=%lf betae=%e\n",
	    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,absorb,gamma,Cmax,IC50,m,total_stats_time,beta_e);
	fprintf(results,
	    "beta=%e latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf absorb=%lf gamma=%lf Cmax=%lf IC50=%lf m=%lf time=%lf betae=%e\n",
	    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,absorb,gamma,Cmax,IC50,m,total_stats_time,beta_e);
    }
    else
    {
	if (results != stdout)
	    fprintf(stdout,
	    "beta=%e latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf cd8_ic50=%lf kappa=%lf alpha=%lf exp_days=%lf hill=%lf time=%lf betae=%e\n",
	    vars->beta_init,vars->latent_inf_init,vars->log_p_init,vars->c_init,vars->theta_init,vars->delta_init,vars->r_init,vars->inf_init,vars->rinf_init,vars->rho_init,vars->eclipse_init,vars->fpos,vars->an,vars->cd8_ic50_init,vars->kappa_init,vars->alpha_init,vars->exp_days_init,vars->hill,total_stats_time,beta_e);

	fprintf(results,
	    "beta=%e latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf cd8_ic50=%lf kappa=%lf alpha=%lf exp_days=%lf hill=%lf time=%lf betae=%e\n",
	    vars->beta_init,vars->latent_inf_init,vars->log_p_init,vars->c_init,vars->theta_init,vars->delta_init,vars->r_init,vars->inf_init,vars->rinf_init,vars->rho_init,vars->eclipse_init,vars->fpos,vars->an,vars->cd8_ic50_init,vars->kappa_init,vars->alpha_init,vars->exp_days_init,vars->hill,total_stats_time,beta_e);
    }

    fprintf(results,"Final percentage of swabs that were positive= %lf\n",totCrit1Perc);
    fprintf(results,"Shedders percentage of swabs that were positive= %lf\n",subCrit1Perc);

    fflush(results);
    *valid=1;
    gsl_matrix_free(pastTps);
    free(actArray);

    return score;	 /* return score for last patient (if multiple ones) */
}

double fitmodel(globalState *vars)
{
    double score;

    int valid=0;

    if (vars->PDF_on)
    {
	do
	    vars->beta_init = vars->beta_mean+gsl_ran_gaussian (vars->ur, vars->beta_std);
	while (vars->beta_init < 0);
	do
	    vars->betae_init = vars->betae_mean+gsl_ran_gaussian (vars->ur, vars->betae_std);
	while (vars->betae_init < 0);
	do
	    vars->latent_inf_init = vars->latent_inf_mean+gsl_ran_gaussian (vars->ur, vars->latent_inf_std);
	while (vars->latent_inf_init < 0);
	do
	    vars->log_p_init = vars->log_p_mean+gsl_ran_gaussian (vars->ur, vars->log_p_std);
	while (vars->log_p_init < vars->log_p_low || vars->log_p_init > vars->log_p_high);
	do
	    vars->c_init = vars->c_mean+gsl_ran_gaussian (vars->ur, vars->c_std);
	while (vars->c_init < 0);
	do
	    vars->theta_init = vars->theta_mean+gsl_ran_gaussian (vars->ur, vars->theta_std);
	while (vars->theta_init < 0);
	do
	    vars->r_init = vars->r_mean+gsl_ran_gaussian (vars->ur, vars->r_std);
	while (vars->r_init < 0);
	do
	    vars->rinf_init = vars->rinf_mean+gsl_ran_gaussian (vars->ur, vars->rinf_std);
	while (vars->rinf_init < 0);
	do
	    vars->inf_init = vars->inf_mean+gsl_ran_gaussian (vars->ur, vars->inf_std);
	while (vars->inf_init < 0);

	do
	    vars->delta_init = vars->delta_mean+gsl_ran_gaussian (vars->ur, vars->delta_std);
	while (vars->delta_init < 0);

	do
	    vars->eclipse_init = vars->eclipse_mean+gsl_ran_gaussian (vars->ur, vars->eclipse_std);
	while (vars->eclipse_init < 0);

	do
	    vars->fpos = vars->fpos_mean + gsl_ran_gaussian(vars->ur,vars->fpos_std);
	while (vars->fpos < 0);

	do
	    vars->hill = vars->hill_mean + gsl_ran_gaussian(vars->ur,vars->hill_std);
	while (vars->hill < 0);

	do
	    vars->an = vars->an_mean + gsl_ran_gaussian(vars->ur,vars->an_std);
	while (vars->an < 0);

	do
	    vars->alpha_init = vars->alpha_mean + gsl_ran_gaussian(vars->ur,vars->alpha_std);
	while (vars->alpha_init < 0);

	do
	    vars->kappa_init = vars->kappa_mean + gsl_ran_gaussian(vars->ur,vars->kappa_std);
	while (vars->kappa_init < 0);

	do
	    vars->exp_days_init = vars->exp_days_mean + gsl_ran_gaussian(vars->ur,vars->exp_days_std);
	while (vars->exp_days_init < 0);

	do
	    vars->cd8_ic50_init = vars->cd8_ic50_mean + gsl_ran_gaussian(vars->ur,vars->cd8_ic50_std);
	while (vars->cd8_ic50_init < 0);

	do
	    vars->To = vars->To_mean + gsl_ran_gaussian(vars->ur,vars->To_std);
	while (vars->To < 0);

	if (vars->Model > 5 || vars->Model_2 > 5)
	{
	    do
		vars->Cmax_init = vars->Cmax_mean+gsl_ran_gaussian (vars->ur, vars->Cmax_std);
	    while (vars->Cmax_init < 0);

	    do
		vars->IC50_init = vars->IC50_mean+gsl_ran_gaussian (vars->ur, vars->IC50_std);
	    while (vars->IC50_init < 0);

	    do
		if (vars->m_mean != 0)
		    vars->m_init = vars->m_mean+gsl_ran_gaussian (vars->ur, vars->m_std);
		else
		    vars->m_init = vars->m_low +gsl_rng_uniform(vars->ur)* (vars->m_high-vars->m_low);
	    while (vars->m_init < 1);

	    do
		if (vars->gamma_mean != 0)
		    vars->gamma_init = vars->gamma_mean+gsl_ran_gaussian (vars->ur, vars->gamma_std);
		else
		    vars->gamma_init = vars->gamma_low +gsl_rng_uniform(vars->ur)* (vars->gamma_high-vars->gamma_low);
	    while (vars->gamma_init < 0);

	    /* absorb handled with range rather than distribution */
	    do
		if (vars->absorb_mean != 0)
		    vars->absorb_init = vars->absorb_mean+gsl_ran_gaussian (vars->ur, vars->absorb_std);
		else
		    vars->absorb_init = vars->absorb_low +gsl_rng_uniform(vars->ur)* (vars->absorb_high-vars->absorb_low);
	    while (vars->absorb_init < 0);
	}
    }
    else
    {
	/* set mean values for results printout */
	    vars->beta_mean = vars->beta_init;
	    vars->betae_mean = vars->betae_init;
	    vars->latent_inf_mean = vars->latent_inf_init;
	    vars->log_p_mean = vars->log_p_init;
	    vars->c_mean = vars->c_init;
	    vars->theta_mean = vars->theta_init;
	    vars->r_mean = vars->r_init;
	    vars->rinf_mean = vars->rinf_init;
	    vars->inf_mean = vars->inf_init;
	    vars->delta_mean = vars->delta_init;
	    vars->fpos_mean = vars->fpos;
	    vars->hill_mean = vars->hill;
	    vars->rho_mean = vars->rho_init;
	    vars->eclipse_mean = vars->eclipse_init;
	    vars->To_mean = vars->To;
	    vars->an_mean = vars->an;
	    if (vars->Model > 5) /* adjust drug related params */
	    {
		vars->Cmax_mean = vars->Cmax_init;
		vars->IC50_mean = vars->IC50_init;
		vars->m_mean = vars->m_init;
		vars->gamma_mean = vars->gamma_init;
		vars->absorb_mean = vars->absorb_init;
	    }
    }

    vars->hex_time_bias = 0;
    vars->sample_index = 0;

    score=ScoreFunction(&valid, (void *)vars,NULL);

    return score;
}

void usage(char *prog_name)
{

    fprintf(stderr,
	"Usage: %s [-h][-d][-f <input_file>][-c <crit file>][-n <cells>[-r][-s <seed>][-v][-w <write_mask>][-W <crit><weight>]\n",prog_name);
    fprintf(stderr,"\t-h = this help\n");
    fprintf(stderr,"\t-f = optional input file\n");
    fprintf(stderr,"\t-c = optional criteria file\n");
    fprintf(stderr,"\t\tFormat: target lower and upper CIs and mean values for...\n");
    fprintf(stderr,"\t\t\t cumul percent pos swabs (9 bins)\n");
    fprintf(stderr,"\t\t\t peak log VL histograms (8 bins)\n");
    fprintf(stderr,"\t-w <write_mask> = which output (csv) files to generate (Bit-mask 1 shifted by 1-12)\n");
    fprintf(stderr,"\t\t1= cumul & regional ve vs time\n");
    fprintf(stderr,"\t\t2= ve current episode episode\n");
    fprintf(stderr,"\t\t3= cumul & regional vi vs time\n");
    fprintf(stderr,"\t\t4= vi current episode episode\n");
    fprintf(stderr,"\t\t5= cumul & regional inf cells vs time\n");
    fprintf(stderr,"\t\t6= inf cells current episode episode\n");
    fprintf(stderr,"\t\t8= regional R0s vs time\n");
    fprintf(stderr,"\t\t9= run stats (per run)\n");
    fprintf(stderr,"\t\t10= episode stats (per run)\n");
    fprintf(stderr,"\t\t11= more detailed regional info\n");
    fprintf(stderr,"\t\t12= ACV dosing info\n");
    fprintf(stderr,"\t-d = deterministic mode (no distributional draws)\n");
    fprintf(stderr,"\t-n = change the size of the hive (max = 1000)\n");
    fprintf(stderr,"\t-W = change weighting of a scoring criteria\n");
    fprintf(stderr,"\t-v = verbose messages\n");
}
int main(int argc, char *argv[])
{
    char def_file[] = "hhv8_sim.in";
    char def_outp_dir[] = ".";

    char dat_file1[] = "hhv8_sim.dat1";
    char dat_file2[] = "hhv8_sim.dat2";
    char dat_file3[] = "hhv8_sim.dat3";
    char dat_file4[] = "hhv8_sim.dat4";
    char dat_file5[] = "hhv8_sim.dat5";
    char dat_file6[] = "hhv8_sim.dat6";
    char dat_file7[] = "hhv8_sim.dat7";
    char dat_file8[] = "hhv8_sim.dat8";
    char dat_file9[] = "hhv8_sim.dat9";
    char dat_file10[] = "hhv8_sim.dat10";
    char dat_file11[] = "hhv8_sim.dat11";
    char dat_file12[] = "hhv8_sim.dat12";

    char criteria_file[] = "hhv8_sim.crit";
    char *crit_file;

    outDir = &def_outp_dir[0];


    long time_st,time_fin;
    const gsl_rng_type * T;//pointer to gsl rand generator type
  
    long seed;

    /* create a genrator chosen by the environment variable GSL_RNG_TYPE */
    gsl_rng_env_setup();

    T = gsl_rng_default;

    globalState settings;

    theState = &settings;

    settings.inp_file = &def_file[0];
    crit_file = &criteria_file[0];

    settings.ur = gsl_rng_alloc (T);

    int Ncrit = NUM_CRITERIA;

    stoch = 1;
    int writeData = 0;
    int writeMask = 0;

    int crit_num = 0;
    double crit_weight = 0;

    double best;

    hexcell *the_cells[MAX_HEXCELLS];
    int next_available = 1;

    int num_cells;
    hexcell::num_hex_cells = 0;

    for(int i = 0; i < argc; i++) {
	cout << "argv[" << i << "] = " << argv[i] << endl;
	if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"-H")) {
	    usage(argv[0]);
	    exit(0);
	}

	if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"-D")) {
	    stoch = 0;
	}
	if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"-V")) {
	    settings.Verbose = 1;
	}
	if (i +1 < argc && !strcmp(argv[i],"-w") && sscanf(argv[i+1],"%d",&writeMask)==1) {
	    writeData = 1;
	    i++;
	}
	if (!strcmp(argv[i],"-W")) {
	    if (i+2 >=argc ||
		(sscanf(argv[i+1],"%d",&crit_num)!=1) ||
		(sscanf(argv[i+2],"%lf",&crit_weight)!=1))
	    {
		fprintf(stderr,"-W option requires two numeric values! Exiting.");
		exit(1);
	    }
	    if (crit_num < 1 || crit_num > MAX_CRIT_CATEGORIES)
	    {
		fprintf(stderr,"-W option requires criteria index between 1 and %d! Exiting.",MAX_CRIT_CATEGORIES);
		exit(1);
	    }
	    settings.critWeight[crit_num-1]=crit_weight;
	    fprintf(stderr,"Set weight for category %d to %lf\n",crit_num-1,crit_weight);
	    i+=2;
	}
	if ((!strcmp(argv[i],"-c") || !strcmp(argv[i],"-C")) && i +1 < argc) {
	    crit_file = argv[i+1];
	    i++;
	}
	if ((!strcmp(argv[i],"-f") || !strcmp(argv[i],"-F")) && i +1 < argc) {
	    settings.inp_file = argv[i+1];
	    i++;
	}
	if ((!strcmp(argv[i],"-a") || !strcmp(argv[i],"-A")) && i +1 < argc) {
	    settings.alt_inp_file = argv[i+1];
	    i++;
	}
	if ((!strcmp(argv[i],"-o") || !strcmp(argv[i],"-O")) && i +1 < argc) {
	    outDir = argv[i+1];
	    i++;
	}
	if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"-R")) {
	    seed = time (NULL) * getpid();    
	    gsl_rng_set (settings.ur, seed);
	    fprintf(stdout,"Random generator seed is %ld\n",seed);
	}
	if (i +1 < argc && !strcmp(argv[i],"-s") && sscanf(argv[i+1],"%ld",&seed)==1) {
	    gsl_rng_set (settings.ur, seed);
	    fprintf(stdout,"Random generator seed is %ld\n",seed);
	    i++;
	}
	if (i +1 < argc && !strcmp(argv[i],"-n") && sscanf(argv[i+1],"%d",&num_cells)==1) {
	    if (num_cells > MAX_HEXCELLS)
	    {
		fprintf(stderr,"Number of requested cells %d exceeds hardcoded limit %d. Exiting.",
		    num_cells,MAX_HEXCELLS);
		exit(1);
	    }

	    hexcell::num_hex_cells = num_cells;
	    i++;
	}
    }

    time(&time_st);
    settings.time_st = time_st;

    read_input_file(0, settings.inp_file, &settings);

    if (hexcell::num_hex_cells == 0)
	hexcell::num_hex_cells = settings.Regions;

    settings.Regions = hexcell::num_hex_cells;

    settings.cells = the_cells;
  
    for (int i=0; i < hexcell::num_hex_cells; i++)
	the_cells[i] = new hexcell(i);

    if (hexcell::num_hex_cells < MAX_HEXCELLS)
	for (int i=hexcell::num_hex_cells; i < MAX_HEXCELLS; i++)
	    the_cells[i] = NULL;

    for (int i=0; i < hexcell::num_hex_cells; i++)
	if(the_cells[i]->create_neighborhood(the_cells,&next_available) == false)
	{
	    fprintf(stderr,"Encountered a problem creating neighborhood for node %d. Exiting.",i);
	    exit(1);
	}

#ifdef ZERO
    for (int i=0; i < hexcell::num_hex_cells; i++)
	the_cells[i]->print_neighborhood();
#endif


    ////////////////////////////////////////////////////////////////////////
    ///// read criteria CIs through a criteria file that has them //////////////
    ifstream critF;
    critF.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit );

    try {
	critF.open(crit_file);
	for(int j=0; j<=Ncrit-1; j++)
	{	
	    critF >> settings.crit[j].low >> settings.crit[j].high>> settings.crit[j].mean;
	    cout <<j<< " low:"<< settings.crit[j].low << " high:"<< settings.crit[j].high<< " mean:"<< settings.crit[j].mean;
	    cout<<endl;
	}	
	cout<<endl;
	critF.close();
    } catch (ifstream::failure e)
    {
	cerr << "Could not open criteria file "<<crit_file<<"\n";
	exit(1);
    }

    settings.dataF1 = NULL;
    settings.dataF2 = NULL;
    settings.dataF3 = NULL;
    settings.dataF4 = NULL;
    settings.dataF5 = NULL;
    settings.dataF6 = NULL;
    settings.dataF7 = NULL;
    settings.dataF8 = NULL;
    settings.dataF9 = NULL;
    settings.dataF10 = NULL;
    settings.dataF11 = NULL;
    settings.dataF12 = NULL;

    if (writeData) {

	settings.writeOn = 1;

	if(((writeMask & 1) && ((settings.dataF1 = fopen(dat_file1,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file1<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<1)) && ((settings.dataF2 = fopen(dat_file2,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file2<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<2)) && ((settings.dataF3 = fopen(dat_file3,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file3<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<3)) && ((settings.dataF4 = fopen(dat_file4,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file4<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<4)) && ((settings.dataF5 = fopen(dat_file5,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file5<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<5)) && ((settings.dataF6 = fopen(dat_file6,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file6<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<6)) && ((settings.dataF7 = fopen(dat_file7,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file7<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<7)) && ((settings.dataF8 = fopen(dat_file8,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file8<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<8)) && ((settings.dataF9 = fopen(dat_file9,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file9<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<9)) && ((settings.dataF10 = fopen(dat_file10,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file10<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<10)) && ((settings.dataF11 = fopen(dat_file11,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file11<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<11)) && ((settings.dataF12 = fopen(dat_file12,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file12<<"\n";
	    exit(1);
	}
    }

    best =fitmodel(&settings);

    if (settings.Verbose)
    {
	fprintf(stderr,"Final score = %lf for the following values...\n\n",best);
    }

    //free memory from the random number generator
    gsl_rng_free (settings.ur);
    time(&time_fin);

    if (writeData)
    {
	if(settings.dataF1 != NULL) fclose(settings.dataF1);
	if(settings.dataF2 != NULL) fclose(settings.dataF2);
	if(settings.dataF3 != NULL) fclose(settings.dataF3);
	if(settings.dataF4 != NULL) fclose(settings.dataF4);
	if(settings.dataF5 != NULL) fclose(settings.dataF5);
	if(settings.dataF6 != NULL) fclose(settings.dataF6);
	if(settings.dataF7 != NULL) fclose(settings.dataF7);
	if(settings.dataF8 != NULL) fclose(settings.dataF8);
	if(settings.dataF9 != NULL) fclose(settings.dataF9);
	if(settings.dataF10 != NULL) fclose(settings.dataF10);
	if(settings.dataF11 != NULL) fclose(settings.dataF11);
	if(settings.dataF12 != NULL) fclose(settings.dataF12);
    }


}//end of main
