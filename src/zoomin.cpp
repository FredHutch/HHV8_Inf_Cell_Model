/* zoomin: this function traverses a solution space by picking
 * the "best" direction from the current point and shrinking the search 
 * space at each new location.  It counts on there being consistent 
 * gradients towards the ideal solution (not always the case).
 *
 * Zoomin was developed as an alternative to "mem_walk" which chooses several 
 * starting locations and "walks" towards the best solution in the neighborhood
 * for each one, returning the best overall solution.  
 *
 */
#include <pthread.h>
#include<cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

/* The follwoing variables are only accessed from the main thread 
 * (hence they can be static) */
static int total_launched = 0;
static int total_expected = 0;

struct ThreadData {
    void *tid;
    int valid;
    double (*F)(int *valid,gsl_vector *params, void *data,FILE *);
    gsl_vector *dumx;
    double tempval;
    void *data;
    FILE *results;
};

void *ThreadedEval(void *details)
{
    struct ThreadData *tdata = (struct ThreadData *)details;
/*
    fprintf(stderr,"Evaluating (X1=%g,X2=%g,X3=%g,X4=%g)\n", tdata->dumx.vec[1],tdata->dumx.vec[2],tdata->dumx.vec[3],tdata->dumx.vec[4]);
*/
    
    tdata->tempval = (*tdata->F)(&tdata->valid,tdata->dumx,tdata->data,tdata->results);
    pthread_exit(tdata->tid);
}


void run_matrix(double (*F)(int *valid,gsl_vector *params, void *data,FILE *),gsl_rng * r,
	int *valid,gsl_vector *dumx, int thisOne, gsl_vector *lsize,gsl_vector *low_bound,gsl_vector *high_bound,
	gsl_vector *bestloc, double *bestval, int param_mask, int max_params,void *data, 
	int num_threads, int *launched, int searchOrder, pthread_t *threads, struct ThreadData *threadData)  
{
    int rc;	/* thread launch return code */
    double dumfl;
    int oneValid;
    double tempval,uni;

    gsl_vector *newx;

    newx = gsl_vector_alloc(max_params);
    gsl_vector_memcpy(newx,dumx);
    
    int mult;

    /* recurse until we get to the param in mask */
    if (((1<<(thisOne)) & param_mask) == 0) {
	run_matrix(F,r,valid,newx,thisOne+1,lsize,
	    low_bound,high_bound,bestloc,bestval,param_mask,
	    max_params,data,num_threads,launched,searchOrder,threads,threadData);
    }
    else {
	for (int i=0; i < 5; i++) {
	    
	    switch (searchOrder) {
		case 0:  /* low to high (-2,-1,0,1,2)*/
		    mult = i-2; 
		    break;
		case 1: /*middle, low, high (0,-1,1,-2,2) */
		    switch (i) { 
			case 0: mult=0; break;
			case 1: mult=-1; break;
			case 2: mult=1; break;
			case 3: mult=-2; break;
			case 4: mult=2; break;
		    }
		    break;
		default: /* high to low (2,1,0,-1,-2)*/
		    mult = 2-i; 
		    break;
	    }
	    /* reset dumx to params each time */
	    dumfl = gsl_vector_get(dumx,thisOne) + ((double)mult)*gsl_vector_get(lsize,thisOne);
	    gsl_vector_set(newx,thisOne,dumfl);

	    /*  check in bounds: */

	    if( ( gsl_vector_get(newx,thisOne) >= gsl_vector_get(low_bound,thisOne) )&&
	       (gsl_vector_get(newx,thisOne) <= gsl_vector_get(high_bound,thisOne) ) ){

		/* is this the last (or "highest") varied parameter in the list? */
		if (param_mask < 1<<(thisOne+1)) {
		    if (num_threads) {
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

			/* for threaded solution, launch one thread per "neioghbor" */
			threadData[*launched].tid=(void *)*launched;
			threadData[*launched].F=F;
			threadData[*launched].dumx=gsl_vector_alloc(max_params);
			gsl_vector_memcpy(threadData[*launched].dumx,newx);
			threadData[*launched].data = data;
			 
			rc = pthread_create(&threads[*launched], &attr, ThreadedEval, (void *)&threadData[*launched]);
			if (rc){
			    fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
			    printf("ERROR; return code from pthread_create() is %d\n", rc);
			}
			(*launched)=(*launched)+1;
			total_launched++;

			fprintf(stderr,"Launched thread %d (%x)...\n",*launched,
			    threadData[*launched-1].tid);

			/* at this point, if we've launched "threade" threads,
			   we wait for all evaluations to rejoin before proceeding! */

			if (num_threads > 0 && (*launched == num_threads ||
			    total_launched == total_expected)/* last one */) {
			    fprintf(stderr,"Waiting for %d threads...\n",*launched);
			    for (int t=0; t < *launched; t++) {
				void *status;
				int j,k;
				int rc;	/* thread launch return code */
				rc = pthread_join(threads[t], &status);
				if (rc){
				    fprintf(stderr,"ERROR; return code from pthread_join() is %d\n", rc);
				    printf("ERROR; return code from pthread_join() is %d\n", rc);
				}
				/* this code (from Eval) is done here to avoid need for locking between threads */
				if( threadData[t].valid) {

				    /* Compare value returned by thread against current max value */
				    /* if no tie...*/
				    /*fprintf(stderr,"For thread id=%g: value=%lf\n",threadData[t].tid,threadData[t].tempval);*/

				    if( threadData[t].tempval > *bestval ){
					*bestval = threadData[t].tempval;
					uni = 0;
					gsl_vector_memcpy(bestloc,threadData[t].dumx);
				    }

				    else{
				      
				    /* but in case of a tie,
				    choose a uniform r.v. and change if larger than before: */

					if(threadData[t].tempval == *bestval){
					    double r1; r1 = gsl_rng_uniform(r);
					    if( r1 > uni){
					      uni = r1;
					      gsl_vector_memcpy(bestloc,threadData[t].dumx);
					    } 
					} 
				    }
				}
				*valid |= threadData[t].valid;
				gsl_vector_free(threadData[t].dumx);
			    }
			    *launched = 0;
			}
		    }
		    else {
			tempval = (*F)(&oneValid,newx,data,NULL);
			if (oneValid) {
			    /* if no tie...*/
			    if( tempval > *bestval ){
				*bestval = tempval;
				uni = 0;
				gsl_vector_memcpy(bestloc,newx);
			    }

			    else{
				/* but in case of a tie,
				choose a uniform r.v. and change if larger than before: */

				if(tempval == *bestval){
				    double r1; r1 = gsl_rng_uniform(r);
				    if( r1 > uni){
					uni = r1;
					gsl_vector_memcpy(bestloc,newx);
				    } 
				} 

			    }
			} /* valid check */
			*valid |= oneValid;
		    } /* high/low check */
		}
		else
		    run_matrix(F,r,valid,newx,thisOne+1,lsize,low_bound,high_bound,bestloc,bestval,param_mask,max_params,data,num_threads,launched, searchOrder,threads,threadData);
	    }
		
	}
    }
    gsl_vector_free(newx);
}

int zoomin(double (*F)( int *,gsl_vector *,void *,FILE *),
    gsl_rng * r, gsl_vector *params,double *Max,
    gsl_vector *low_bound,gsl_vector *high_bound,
    int param_mask, int max_params,
    int max_steps,int stop_walk,int bvstop_walk,
    int print, double tolerance,void *data, int num_threads, int searchOrder)

{
    int i,step,end,memlength,memvisits,unit,intlsize, dumint;
    double tiny,currentval,bestval,dumfl,dumfl1;
    gsl_vector *bestloc,*lsize,*store_lsize,*dumx;
    int num_params=0;
    int launched=0;

    total_launched=0;

    int valid;    /* flag to tell if return value from F is ok */
   
    pthread_t *threads = NULL;
    struct ThreadData *threadData;

    /* count set bits in param_mask */
    for (int bit=0; bit < max_params; bit++)
	if ((1<<bit) & param_mask)
	    num_params++;


    total_expected = (int)(pow(5.0,num_params)+0.5);

    fprintf(stderr,"Varying %d parameters (%d runs expected)\n",
		num_params,total_expected);

    /* allocate the thread structures (for zoomin max parallelism is 3)*/
    if (num_threads > 0) {
       threads = (pthread_t *) malloc (num_threads*sizeof (pthread_t));
       if( threads == NULL) fprintf(stderr,"malloc error in zoomin: could not allocate threads"); 
       threadData = (struct ThreadData *) malloc (num_threads*sizeof (struct ThreadData));
       if( threadData == NULL) fprintf(stderr,"malloc error in zoomin: could not allocate threads"); 
       for (int thr=0; thr < num_threads;thr++)
        {
	    char filename[20];
	    sprintf(filename,"hhv8_sim.out_%d",thr);
	    if( (threadData[thr].results = fopen(filename,"wt")) == NULL){
		fprintf(stderr,"Could not results output file %s\n",filename);
		exit(1);
	    }
       	}
    }

    dumx=gsl_vector_alloc(max_params);
    lsize=gsl_vector_alloc(max_params);
    bestloc=gsl_vector_alloc(max_params);
    store_lsize=gsl_vector_alloc(max_params);

    unit = 1024; tiny = 1.0e-8;

   /* initialise walker to check mid-values and "neighborhood" +/- 33% */

    for(i=0;i<max_params;i++){
	if ((1<<(i)) & param_mask) {
	    dumfl = ( gsl_vector_get(high_bound,i) + gsl_vector_get(low_bound,i) )/2;
	    gsl_vector_set(dumx,i,dumfl);
	    dumfl = ( gsl_vector_get(high_bound,i) - gsl_vector_get(low_bound,i) )/5;
	    gsl_vector_set(lsize,i,dumfl);
	} else {
	    dumfl = gsl_vector_get(params,i);
	    gsl_vector_set(dumx,i,dumfl);
	    gsl_vector_set(lsize,i,dumfl);
	}
    }

    gsl_vector_memcpy(store_lsize,lsize);

    /* initialize  bestval, bestloc and currentval */

    bestval = -1.0e+30;
    currentval = -1.0e+30;

    step=0;end=0;

    int tstep,quit,bvquit;
    double tempmax,savebestval;

    tstep=0;quit=0;bvquit=0;

    /* 
     * stop the walk if there have been "quit" consective non-improving moves or 
     * the maximum number of zoomin "steps" have been reached
    */
    while( (tstep < max_steps)&&(quit < stop_walk)&&(bvquit < bvstop_walk) ){

	double tempval,uni,sbestval;
	sbestval = bestval;

	  /* propose a move: */

	uni=0.0;

	/* Push evaluation sites onto a list to be performed in parallel by worker threads */
	/* 2*d threads will work in parallel and then report back to see which is best route */
	if (num_threads > 0) {
	    pthread_attr_t attr;
	    pthread_attr_init(&attr);
	    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	    valid=0;
	    run_matrix(F,r,&valid,dumx,0,lsize,low_bound,high_bound,bestloc,&bestval,param_mask,max_params,data,num_threads,&launched,searchOrder,threads,threadData);

	} /*threaded approach */

	else {
	    /* non-threaded approach */

	    valid=0;
	    run_matrix(F,r,&valid,dumx,0,lsize,low_bound,high_bound,bestloc,&bestval,param_mask,max_params,data,num_threads,&launched,searchOrder,threads,threadData);

	} /* end of non-threaded option */

	if (valid) {

	    /* at this point we have updated bestval and bestloc based on nearest neighbors
	    (whether done in parallel or sequentially) */

	    if(bestval >= currentval){
		/* move: */

		gsl_vector_memcpy(params,bestloc);
		gsl_vector_scale(lsize,0.5);


		/* was there improvement? */

		if( bestval - currentval < currentval*tolerance) quit++;

		else quit = 0;

		currentval = bestval;
	    }


	    else{ /* stay here, but shrink lsize by half */
		gsl_vector_scale(lsize,0.5);

		quit++;
	    }

	    /* Did this walk make a global sufficient improvement? */

	    if(bestval - sbestval < sbestval*tolerance) bvquit++;
	    else bvquit = 0;

	    tstep++;
	} else {
	    fprintf(stderr,"Abandoning current walk (all options bogus)!\n");
	    fflush(stderr);
	    break;
	}
   } /* end walk loop */

if(print==1){
  fprintf(stderr,"walk length = %i\n",tstep++); 
  fprintf(stderr,"Best value on current walk = %f\n",currentval);
 } 
  (*Max) = bestval;
   gsl_vector_memcpy(params,bestloc);

  
  gsl_vector_free(lsize);
  gsl_vector_free(store_lsize);
  gsl_vector_free(bestloc);
  gsl_vector_free(dumx);

  /* allocate the thread structures */
  if (num_threads > 0) {
       for (int thr=0; thr < num_threads;thr++)
        {
	    if(threadData[thr].results != NULL)
		fclose(threadData[thr].results);
	}
     if( threads != NULL) free(threads);
     if( threadData != NULL) free (threadData);
  }
  return valid;

}
