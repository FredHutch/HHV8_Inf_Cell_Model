#ifndef PLOTPOINTS_H
#define PLOTPOINTS_H

#define MAX_POINTS 1000

class plotPoints {

public:
	double *time;

	long int *vet;
	long int *vit;
	long int **ve;
	long int **vi;

	long int *inf;
	long int *cd8cells;
	double *ACV;
	int **color;

   	int valid;
   	int max_points;

	plotPoints(void);
	~plotPoints(void);

	void checkForRealloc();
};
#endif
