#ifndef HEXCELL_H
#define HEXCELL_H

#define MAX_HEXCELLS 1000
#define MAX_NEIGHBORS 6

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

#endif
