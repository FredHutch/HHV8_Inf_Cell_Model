#ifndef HEXCELL_H
#define HEXCELL_H

#define MAX_HEXCELLS 1000
#define MAX_NEIGHBORS 6
#define MAX_PLAQ_COLORS 7

#define MAX_VL_LOG 9.0
#define MAX_INF_LOG 6.0
#define MAX_CD8_LOG 3.8
#define MIN_CD8_LOG 3.0
#define MAX_CD8 7000.
#define MIN_CD8 1000.0

#define ACTUAL_HEX_DIAM 5.9
#define DTOR 0.0174532925

#define MAX_DIAMETER 3.0
#define MAX_REPRO 16.0

class globalState;

#include <GL/gl.h>
#include <GL/glu.h>

class hexcell_ref {
public:
	bool sent;
	int cell;
	hexcell_ref(int index);
	~hexcell_ref(void){};
};

class hexcell {

public:
	hexcell_ref *neighbors[MAX_NEIGHBORS];
	void set_coords(double newx, double newy);
	bool add_neighbors(hexcell_ref *from, hexcell_ref *left, hexcell_ref *right, int edge, float xcoord, float ycoord); 
	bool create_neighborhood(hexcell *the_cells[],int *start_number);
	void print_neighborhood(void);
	//void create_circle(double r);
	double graph_cell(globalState *the_state, int option, GLfloat *pColors[]);

	hexcell(int index);
	~hexcell(void);

	int num_neighbors;
	hexcell_ref *myref;

	float xcoord;
	float ycoord;

	static double max_x_vertex;
	static double max_y_vertex;

	static int num_hex_cells;

private:
	static double x_offsets[MAX_NEIGHBORS];
	static double y_offsets[MAX_NEIGHBORS];
	static bool offsets_done;

	double x_vertices[MAX_NEIGHBORS];
	double y_vertices[MAX_NEIGHBORS];
};

static GLfloat Black[4] = {0,0,0,1};
static GLfloat Red[4] = {1, 0.2, 0.2, 1};
static GLfloat Pink[4] = {1, 0.5, 0.5, 1};
static GLfloat Green[4] = {0.2, 1, 0.2, 1};
static GLfloat DarkBlue[4] = {0., 0., 0.65, 1};
static GLfloat DarkRed[4] = {0.65, 0., 0., 1};
static GLfloat DarkGreen[4] = {0.0, 0.65, 0.0, 1};
static GLfloat DarkPurple[4] = {0.44, 0., 0.73, 1};
static GLfloat DarkPurple2[4] = {0.17, 0., 0.34, 1};
static GLfloat Purple[4] = {0.57, 0.44, 0.86, 1};
static GLfloat Yellow[4] = {1, 1, 0.2, 1};
static GLfloat White[4] = {1,1,1,1};
static GLfloat Gray[4] = {.75,.75,1,1};
static GLfloat Orange[4] = {0.75, 0.375, 0.0, 1};
static GLfloat Gold[4] = {1.0, 211./255., 32./255., 1};
static GLfloat Brown[4] = {75./255., 25./255., 0.0, 1};
static GLfloat DarkBrown[4] = {37./255., 12./255., 0.0, 1};
static GLfloat Brick[4] = {126./255., 0.0, 33./255., 1};
#endif
