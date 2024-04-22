#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include <cstdlib>
#include<cmath>
using namespace std;

#include <strings.h>

#include "hexcell.h"
#include "global_state.h"

double hexcell::x_offsets[MAX_NEIGHBORS]={0,0,0,0,0,0};
double hexcell::y_offsets[MAX_NEIGHBORS]={0,0,0,0,0,0};

double hexcell::max_x_vertex=0;
double hexcell::max_y_vertex=0;

int hexcell::num_hex_cells=MAX_HEXCELLS;

bool hexcell::offsets_done = false;

hexcell_ref::hexcell_ref(int index) 
{
    cell=index;
    sent=false;
}

hexcell::hexcell(int current) 
{
    /* start with no neighbors */
    if (offsets_done == false)
    {
	for (int i=0; i < MAX_NEIGHBORS; i++)
	{
	/* setup vertx offsets (these apply to all cells) */
	    switch(i)
	    {
		case 0: 
		    x_offsets[0]=0.0;
		    y_offsets[0]=1.0 / pow(3.0,0.5);
		    break;
		case 1: 
		    x_offsets[1]=0.5;
		    y_offsets[1]=0.5 / pow(3.0,0.5);
		    break;
		case 2: 
		    x_offsets[2]=0.5;
		    y_offsets[2]=-0.5 / pow(3.0,0.5);
		    break;
		case 3: 
		    x_offsets[3]=0.0;
		    y_offsets[3]=-1.0 / pow(3.0,0.5);
		    break;
		case 4: 
		    x_offsets[4]=-0.5;
		    y_offsets[4]=-0.5 / pow(3.0,0.5);
		    break;
		case 5: 
		    x_offsets[5]=-0.5;
		    y_offsets[5]=0.5 / pow(3.0,0.5);
		    break;
		default:
		    fprintf(stderr,"Code does not handle more than 6 neighbors!\n");
		    exit(1);
	    }
	}
	offsets_done=true;
    }

    for (int i=0; i < MAX_NEIGHBORS; i++)
    {
	neighbors[i] = NULL;

	x_vertices[i] = x_offsets[i];
	if (abs(x_vertices[i]) > max_x_vertex)
	    max_x_vertex = abs(x_vertices[i]) ;

	y_vertices[i] = y_offsets[i];
	if (abs(y_vertices[i]) > max_y_vertex)
	    max_y_vertex = abs(y_vertices[i]) ;
    }

    myref = new hexcell_ref(current);
    myref->sent=true;

    num_neighbors=0;
    xcoord=0;
    ycoord=0;
}

hexcell::~hexcell(void) 
{
    /* delete neighbors */
    for (int i=0; i < MAX_NEIGHBORS; i++)
	if (neighbors[i] != NULL)
	    delete neighbors[i];

    delete myref;
}

void hexcell::set_coords(double newx, double newy) 
{
    xcoord=newx;
    ycoord=newy;
    for (int i=0; i < MAX_NEIGHBORS; i++)
    {
	x_vertices[i] = xcoord+x_offsets[i];
	if (abs(x_vertices[i]) > max_x_vertex)
	    max_x_vertex = abs(x_vertices[i]) ;

	y_vertices[i] = ycoord+y_offsets[i];
	if (abs(y_vertices[i]) > max_y_vertex)
	    max_y_vertex = abs(y_vertices[i]) ;
    }
}

bool hexcell::create_neighborhood(hexcell *the_cells[],int *start_number)
{
    int left,right;
    float newx,newy;

    double longside=pow(3,0.5)/2; 
    double shortside=0.5;
    double diag = 1.0;

    int current=*start_number;

    /* first fill out our neighbor list */
    for (int i=0; i < MAX_NEIGHBORS; i++)
    {
	if (neighbors[i] == NULL && current < num_hex_cells) 
	{
	    neighbors[i] = new hexcell_ref(current);
	    current++;
	    num_neighbors++;
	}
    }

    /* next, send add_neighbors message to each non-sent neighbor */
    for (int i=0; i < MAX_NEIGHBORS; i++)
    {
	if (neighbors[i] != NULL && !neighbors[i]->sent)
	{
	    if (i==0)
		left = MAX_NEIGHBORS-1;
	    else
		left=i-1;

	    if (i== MAX_NEIGHBORS-1)
		right = 0;
	    else
		right=i+1;

	    switch(i)
	    {
		case 0: 
		    newx=xcoord+shortside;
		    newy=ycoord+longside;
		    break;
		case 1: 
		    newx=xcoord+diag;
		    newy=ycoord;
		    break;
		case 2: 
		    newx=xcoord+shortside;
		    newy=ycoord-longside;
		    break;
		case 3: 
		    newx=xcoord-shortside;
		    newy=ycoord-longside;
		    break;
		case 4: 
		    newx=xcoord-diag;
		    newy=ycoord;
		    break;
		case 5: 
		    newx=xcoord-shortside;
		    newy=ycoord+longside;
		    break;
		default:
		    fprintf(stderr,"Code does not handle more than 6 neighbors!\n");
		    exit(1);
	    }

	    if(the_cells[neighbors[i]->cell]->add_neighbors(myref, neighbors[left], neighbors[right],i,newx,newy) == false)
		return false;
	}
    }
    *start_number= current;
    return true;
}

bool hexcell::add_neighbors(hexcell_ref *from, hexcell_ref *left, hexcell_ref *right, int edge, float newx, float newy)
{
    int found_from=-1;
    int found_left=-1;
    int found_right=-1;

    int new_from=(edge+(MAX_NEIGHBORS/2))%MAX_NEIGHBORS;
    int new_left=(edge+(MAX_NEIGHBORS/2)+1)%MAX_NEIGHBORS;
    int new_right=(edge+(MAX_NEIGHBORS/2)-1)%MAX_NEIGHBORS;

    set_coords(newx,newy);

    /* is "from" al;ready in our neighbor list? */
    for (int i=0; i < MAX_NEIGHBORS; i++)
    {
        if (neighbors[i] != NULL && from != NULL && neighbors[i]->cell == from->cell )
	    found_from=i;
        if (neighbors[i] != NULL && left != NULL && neighbors[i]->cell == left->cell )
	    found_left=i;
        if (neighbors[i] != NULL && right != NULL && neighbors[i]->cell == right->cell )
	    found_right=i;
    }

    if (found_from < 0)
	num_neighbors++;
    else if (new_from != found_from)
    {
	fprintf(stderr, "Mismatch between existing order and new neighbors from %d (from edge=%d vs. existing %d)\n",from->cell,new_from,found_from);
	return false;
    }

    if (found_left < 0)
	num_neighbors++;
    else if (new_left != found_left)
    {
	fprintf(stderr, "Mismatch between existing order and new neighbors from %d (left edge=%d vs. existing %d)\n",from->cell,new_left,found_left);
	return false;
    }

    if (found_right < 0)
	num_neighbors++;
    else if (new_right != found_right)
    {
	fprintf(stderr, "Mismatch between existing order and new neighbors from %d (right edge=%d vs. existing %d)\n",from->cell,new_right,found_right);
	return false;
    }

    if (neighbors[new_from] == NULL && from != NULL)
	neighbors[new_from] = new hexcell_ref(from->cell);

    if (neighbors[new_from] != NULL && from != NULL)
    {
	neighbors[new_from]->cell = from->cell;
	neighbors[new_from]->sent = true;
    }

    if (neighbors[new_left] == NULL && left != NULL)
	neighbors[new_left] = new hexcell_ref(left->cell);

    if (neighbors[new_left] != NULL && left != NULL)
	neighbors[new_left]->cell = left->cell;

    if (neighbors[new_right] == NULL && right != NULL)
	neighbors[new_right] = new hexcell_ref(right->cell);

    if (neighbors[new_right] != NULL && right != NULL)
	neighbors[new_right]->cell = right->cell;

    return true;
}

void hexcell::print_neighborhood(void)
{
	printf("Cell %d:",(myref!=NULL)?myref->cell:-1);
	for (int i=0; i < MAX_NEIGHBORS; i++)
	{
	    printf("edge %d: %d",i,(neighbors[i]!=NULL)?neighbors[i]->cell:-1);
	    if (i < MAX_NEIGHBORS-1)
		printf(", ");
	}
	printf(" coords: %g,%g",xcoord,ycoord);
	printf("\n");
}
