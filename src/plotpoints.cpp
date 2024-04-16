#include <cstdlib>
using namespace std;

#include <strings.h>
#include "plotpoints.h"
#include "hexcell.h"

plotPoints::plotPoints(void) 
{
	valid = 0;
	max_points=MAX_POINTS;
	time = (double*) malloc(MAX_POINTS * sizeof(double));
	cd8cells = (long int *)malloc(MAX_POINTS * sizeof (long int));
	vet = (long int *)malloc(MAX_POINTS * sizeof (long int));
	vit = (long int *)malloc(MAX_POINTS * sizeof (long int));
	ve = (long int **)malloc(MAX_HEXCELLS * sizeof (long int*));
	vi = (long int **)malloc(MAX_HEXCELLS * sizeof (long int*));
	color = (int **)malloc(MAX_HEXCELLS * sizeof (int *));
	vet[0]=0;
	vit[0]=0;

	for (int i=0; i < MAX_HEXCELLS; i++)
	{
	    ve[i]=(long int *)malloc(MAX_POINTS * sizeof (long int));
	    ve[i][0]=0;
	    vi[i]=(long int *)malloc(MAX_POINTS * sizeof (long int));
	    vi[i][0]=0;
	    color[i]=(int *)malloc(MAX_POINTS * sizeof (int));
	    color[i][0]=0;
	}

	inf = (long int *)malloc(MAX_POINTS * sizeof (long int));
	inf[0]=0;
	ACV = (double *)malloc(MAX_POINTS * sizeof (double));
	ACV[0]=0;
};

plotPoints::~plotPoints(void) 
{
	free(vet);
	free(vit);
	free(cd8cells);
	free(time);
	free(inf);
	free(ACV);

	for (int i=0; i < MAX_HEXCELLS; i++)
	{
	    free(ve[i]);
	    free(vi[i]);
	    free(color[i]);
	}
	free(ve);
	free(vi);
	free(color);
}

void plotPoints::checkForRealloc()
{
	/* check for need to realloc point array */
	if (valid >= max_points-1) {
	    double *tempTimes = time;
	    time = (double*) malloc(max_points*2 * sizeof(double));
	    bcopy((char *)tempTimes,(char *)time,
		    max_points*sizeof(double));

	    free(tempTimes);

	    long int *tempCd8s = cd8cells;
	    cd8cells = (long int *)malloc(max_points*2 * sizeof (long int));
	    bcopy((char *)tempCd8s,(char *)cd8cells,
		    max_points*sizeof(long int));

	    free(tempCd8s);

	    long int *tempVL = vet;
	    vet = (long int *)malloc(max_points*2 * sizeof (long int));
	    bcopy((char *)tempVL,(char *)vet,
		    max_points*sizeof(long int));

	    free(tempVL);

	    tempVL = vit;
	    vit = (long int *)malloc(max_points*2 * sizeof (long int));
	    bcopy((char *)tempVL,(char *)vit,
		    max_points*sizeof(long int));

	    free(tempVL);

	    long int *tempInf = inf;
	    inf = (long int *)malloc(max_points*2 * sizeof (long int));
	    bcopy((char *)tempInf,(char *)inf,
		    max_points*sizeof(long int));

	    free(tempInf);

	    double *tempACV = ACV;
	    ACV = (double *)malloc(max_points*2 * sizeof (double));
	    bcopy((char *)tempACV,(char *)ACV,
		    max_points*sizeof(double));

	    free(tempACV);

	    for (int i=0; i < MAX_HEXCELLS; i++) {
		long int *tempVe = ve[i];
		ve[i]=(long int *)malloc(max_points*2 * sizeof (long int));
		bcopy((char *)tempVe,(char *)ve[i],
		    max_points*sizeof(long int*));
		free(tempVe);
	    }

	    for (int i=0; i < MAX_HEXCELLS; i++) {
		long int *tempVi = vi[i];
		vi[i]=(long int *)malloc(max_points*2 * sizeof (long int));
		bcopy((char *)tempVi,(char *)vi[i],
		    max_points*sizeof(long int*));
		free(tempVi);
	    }

	    for (int i=0; i < MAX_HEXCELLS; i++) {
		int *tempColor = color[i];
		color[i]=(int *)malloc(max_points*2 * sizeof (int));
		bcopy((char *)tempColor,(char *)color[i],
		    max_points*sizeof(int*));
		free(tempColor);
	    }

	    max_points *= 2;
	}
}
