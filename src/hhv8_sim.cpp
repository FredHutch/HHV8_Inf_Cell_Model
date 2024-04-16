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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include "GL/osmesa.h"
#include <png.h>

#ifndef NO_GUI
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#endif

#include "state_var.h"  
#include "global_state.h"  

#include "hexcell.h"

#include "hhv8.xpm"

using namespace std;

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define DEFAULT_WIDTH  800
#define DEFAULT_HEIGHT 600
#define DEFAULT_IMAGE_WIDTH  816
#define DEFAULT_IMAGE_HEIGHT 592
#ifdef __sun
#define DEFAULT_TITLE  "HHV8 simulation (for SunOS)"
#else
#define DEFAULT_TITLE  "HHV8 simulation"
#endif
#define PROGRAM_NAME  "hhv8_sim"
#define PACKAGE_VERSION  "1.0"
#define PACKAGE_BUGREPORT  "dswan@fredhutch.org"

#define ORIGIN_X          0
#define ORIGIN_Y          0
#define MAX_X_COORD          80.
#define MAX_Y_COORD          80.

#define MAX_RGB_VAL 255
#define MAX_SWABS 10000

#define VE_GRAPH 1
#define VI_GRAPH 2
#define INF_GRAPH 3
#define CD8_GRAPH 4
#define ACV_GRAPH 5

static globalState *theState=NULL;
static int batchMode;
static int stoch;
static int drawnTime = 0;

static void *bgr;
static bool draw_routine (GLfloat width, GLfloat height);
void write_png_file( int width, int height,const char* file_name);
static string outDir;

#ifndef NO_GUI
static GLuint font_list_base;

static gchar default_font_string[] = "courier bold 12";

static GLuint default_font_list_base;

static PangoFontDescription *default_font_desc = NULL;

static GdkGLContext *global_glcontext;
static GdkGLDrawable *global_gldrawable;

static GtkWidget *main_window;

static GtkWidget *image;
static GdkPixmap *pixmap;
static GdkGLConfig *glconfig;

static GtkWidget *beta_entry;
static GtkWidget *c_entry;
static GtkWidget *delta_entry;
static GtkWidget *eclipse_entry;
static GtkWidget *inf_entry;
static GtkWidget *p_entry;
static GtkWidget *r_entry;
static GtkWidget *rho_entry;
static GtkWidget *rinf_entry;
static GtkWidget *theta_entry;
static GtkWidget *latent_inf_entry;
static GtkWidget *plot_span_entry;
static GtkWidget *T0_entry;
static GtkWidget *use_rho_entry;

static GtkWidget *backup_slider;
static GtkObject *adjuster;

static GtkWidget *font_name_widget;
static GtkWidget *font_size_widget;

static GtkWidget *background_r;
static GtkWidget *background_g;
static GtkWidget *background_b;
static GtkWidget *rgb_box;

static GtkWidget *zoomin_button;
static GtkWidget *zoomout_button;
static GtkWidget *run_button;
static GtkWidget *stop_button;
static GtkWidget *pause_button;
static GtkWidget *page_left_button;
static GtkWidget *page_right_button;

static GLfloat da_width = 0.0;
static GLfloat da_height = 0.0;

/* Movie parameters */
struct capture_t {
	GtkWidget *movie_name;
	GtkWidget *codec;
	GtkWidget *framerate;
	GtkWidget *bitrate;
	int pipefd;
	int enc_width, enc_height;
};
struct capture_t scapture;


/**************************************************************************
 *  * The following section contains the function prototype declarations.
 *   **************************************************************************/

static GdkGLConfig *configure_gl      (void);
void updateGL(void);
void snap_movie_frame(void);

static GtkWidget   *create_window     (GdkGLConfig *glconfig);
static void stop_cb (GtkWidget* widget, gpointer data);

#ifndef __sun
void wait_cursor( GdkWindow *win)
{
    GdkCursor *cur;
    cur = gdk_cursor_new( GDK_CLOCK );
    gdk_window_set_cursor( win, cur );
    gdk_cursor_unref( cur );
}

void normal_cursor( GdkWindow *win)
{
    gdk_window_set_cursor( win, NULL );
}
#endif

#else
// FreeType Headers
#include <ft2build.h>
#include <freetype/freetype.h>
#include <freetype/ftglyph.h>
#include <freetype/ftoutln.h>
#include <freetype/fttrigon.h>

// Inside Of This Namespace, Give Ourselves The Ability
// To Write Just "vector" Instead Of "std::vector"
using std::vector;
 
// Ditto For String.
using std::string;
 
// This Holds All Of The Information Related To Any
// FreeType Font That We Want To Create. 
class font_data {
public:
    float h;                                        // Holds The Height Of The Font.
    GLuint * textures;                                  // Holds The Texture Id's
    GLuint list_base;                                   // Holds The First Display List Id
 
    // The Init Function Will Create A Font With
    // The Height h From The File fname.
    void init(const char * fname, unsigned int h);
 
    // Free All The Resources Associated With The Font.
        void clean();

	font_data(void){ textures=NULL;}
	~font_data(void){ }
};

// The Flagship Function Of The Library - This Thing Will Print
// Out Text At Window Coordinates X, Y, Using The Font ft_font.
// The Current Modelview Matrix Will Also Be Applied To The Text.
void print(const font_data &ft_font, float x, float y, const char *fmt, ...);
// This Function Gets The First Power Of 2 >= The
// Int That We Pass It.
inline int next_p2 (int a )
{
    int rval=1;
    // rval<<=1 Is A Prettier Way Of Writing rval*=2;
    while(rval<a) rval<<=1;
    return rval;
}
 
// Create A Display List Corresponding To The Given Character.
void make_dlist ( FT_Face face, char ch, GLuint list_base, GLuint * tex_base ) {
 
    // The First Thing We Do Is Get FreeType To Render Our Character
    // Into A Bitmap.  This Actually Requires A Couple Of FreeType Commands:
 
    // Load The Glyph For Our Character.
    if(FT_Load_Glyph( face, FT_Get_Char_Index( face, ch ), FT_LOAD_DEFAULT ))
        throw std::runtime_error("FT_Load_Glyph failed");
 
    // Move The Face's Glyph Into A Glyph Object.
    FT_Glyph glyph;
    if(FT_Get_Glyph( face->glyph, &glyph ))
        throw std::runtime_error("FT_Get_Glyph failed");
 
    // Convert The Glyph To A Bitmap.
    FT_Glyph_To_Bitmap( &glyph, ft_render_mode_normal, 0, 1 );
    FT_BitmapGlyph bitmap_glyph = (FT_BitmapGlyph)glyph;
 
    // This Reference Will Make Accessing The Bitmap Easier.
    FT_Bitmap& bitmap=bitmap_glyph->bitmap;

    // Use Our Helper Function To Get The Widths Of
    // The Bitmap Data That We Will Need In Order To Create
    // Our Texture.
    int width = next_p2( bitmap.width );
    int height = next_p2( bitmap.rows );
     
    // Allocate Memory For The Texture Data.
    GLubyte* expanded_data = new GLubyte[ 2 * width * height];
     
    // Here We Fill In The Data For The Expanded Bitmap.
    // Notice That We Are Using A Two Channel Bitmap (One For
    // Channel Luminosity And One For Alpha), But We Assign
    // Both Luminosity And Alpha To The Value That We
    // Find In The FreeType Bitmap.
    // We Use The ?: Operator To Say That Value Which We Use
    // Will Be 0 If We Are In The Padding Zone, And Whatever
    // Is The FreeType Bitmap Otherwise.
    for(int j=0; j <height;j++) {
	for(int i=0; i < width; i++){
	    expanded_data[2*(i+j*width)]= expanded_data[2*(i+j*width)+1] =
		(i>=bitmap.width || j>=bitmap.rows) ?
		0 : bitmap.buffer[i + bitmap.width*j];
	}
    }
// Now We Just Setup Some Texture Parameters.
glBindTexture( GL_TEXTURE_2D, tex_base[ch]);
glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
 
// Here We Actually Create The Texture Itself, Notice
// That We Are Using GL_LUMINANCE_ALPHA To Indicate That
// We Are Using 2 Channel Data.
glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0,
    GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, expanded_data );
 
// With The Texture Created, We Don't Need The Expanded Data Anymore.
delete [] expanded_data;

   // Now We Create The Display List
    glNewList(list_base+ch,GL_COMPILE);
 
    glBindTexture(GL_TEXTURE_2D,tex_base[ch]);
 
    glPushMatrix();
 
    // First We Need To Move Over A Little So That
    // The Character Has The Right Amount Of Space
    // Between It And The One Before It.
    glTranslatef(bitmap_glyph->left,0,0);
 
    // Now We Move Down A Little In The Case That The
    // Bitmap Extends Past The Bottom Of The Line
    // This Is Only True For Characters Like 'g' Or 'y'.
    glTranslatef(0,bitmap_glyph->top-bitmap.rows,0);
 
    // Now We Need To Account For The Fact That Many Of
    // Our Textures Are Filled With Empty Padding Space.
    // We Figure What Portion Of The Texture Is Used By
    // The Actual Character And Store That Information In
    // The x And y Variables, Then When We Draw The
    // Quad, We Will Only Reference The Parts Of The Texture
    // That Contains The Character Itself.
    float   x=(float)bitmap.width / (float)width,
    y=(float)bitmap.rows / (float)height;
 
    // Here We Draw The Texturemapped Quads.
    // The Bitmap That We Got From FreeType Was Not
    // Oriented Quite Like We Would Like It To Be,
    // But We Link The Texture To The Quad
    // In Such A Way That The Result Will Be Properly Aligned.
    glBegin(GL_QUADS);
    glTexCoord2d(0,0); glVertex2f(0,bitmap.rows);
    glTexCoord2d(0,y); glVertex2f(0,0);
    glTexCoord2d(x,y); glVertex2f(bitmap.width,0);
    glTexCoord2d(x,0); glVertex2f(bitmap.width,bitmap.rows);
    glEnd();
    glPopMatrix();
    glTranslatef(face->glyph->advance.x >> 6 ,0,0);
 
    // Increment The Raster Position As If We Were A Bitmap Font.
    // (Only Needed If You Want To Calculate Text Length)
    // glBitmap(0,0,0,0,face->glyph->advance.x >> 6,0,NULL);
 
    // Finish The Display List
    glEndList();
}

font_data our_font;

void font_data::init(const char * fname, unsigned int h) {
    // Allocate Some Memory To Store The Texture Ids.
    textures = new GLuint[128];
 
    this->h=h;
 
    // Create And Initilize A FreeType Font Library.
    FT_Library library;
    if (FT_Init_FreeType( &library ))
        throw std::runtime_error("FT_Init_FreeType failed");
 
    // The Object In Which FreeType Holds Information On A Given
    // Font Is Called A "face".
    FT_Face face;
 
    // This Is Where We Load In The Font Information From The File.
    // Of All The Places Where The Code Might Die, This Is The Most Likely,
    // As FT_New_Face Will Fail If The Font File Does Not Exist Or Is Somehow Broken.
    if (FT_New_Face( library, fname, 0, &face ))
        throw std::runtime_error("FT_New_Face failed (there is probably a problem with your font file)");
 
    // For Some Twisted Reason, FreeType Measures Font Size
    // In Terms Of 1/64ths Of Pixels.  Thus, To Make A Font
    // h Pixels High, We Need To Request A Size Of h*64.
    // (h << 6 Is Just A Prettier Way Of Writing h*64)
    FT_Set_Char_Size( face, h << 6, h << 6, 96, 96);
 
    // Here We Ask OpenGL To Allocate Resources For
    // All The Textures And Display Lists Which We
    // Are About To Create. 
    list_base=glGenLists(128);
    glGenTextures( 128, textures );
 
    // This Is Where We Actually Create Each Of The Fonts Display Lists.
    for(unsigned char i=0;i<128;i++)
        make_dlist(face,i,list_base,textures);
 
    // We Don't Need The Face Information Now That The Display
    // Lists Have Been Created, So We Free The Assosiated Resources.
    FT_Done_Face(face);
 
    // Ditto For The Font Library.
    FT_Done_FreeType(library);
}

void font_data::clean() {
    glDeleteLists(list_base,128);
    glDeleteTextures(128,textures);
    delete [] textures;
}

// A Fairly Straightforward Function That Pushes
// A Projection Matrix That Will Make Object World
// Coordinates Identical To Window Coordinates.
inline void pushScreenCoordinateMatrix() {
    glPushAttrib(GL_TRANSFORM_BIT);
    GLint   viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(viewport[0],viewport[2],viewport[1],viewport[3],1.0,-1.0);
    glPopAttrib();
}
 
// Pops The Projection Matrix Without Changing The Current
// MatrixMode.
inline void pop_projection_matrix() {
    glPushAttrib(GL_TRANSFORM_BIT);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glPopAttrib();
}

// Much Like NeHe's glPrint Function, But Modified To Work
// With FreeType Fonts.
void print(const font_data &ft_font, float x, float y, const char *fmt, ...)  {
         
    // We Want A Coordinate System Where Distance Is Measured In Window Pixels.
    pushScreenCoordinateMatrix();                                  
         
    GLuint font=ft_font.list_base;
    // We Make The Height A Little Bigger.  There Will Be Some Space Between Lines.
    float h=ft_font.h/.63f;                                                
    char    text[256];                                  // Holds Our String
    va_list ap;                                     // Pointer To List Of Arguments
 
    if (fmt == NULL)                                    // If There's No Text
        *text=0;                                    // Do Nothing
    else {
        va_start(ap, fmt);                              // Parses The String For Variables
        vsprintf(text, fmt, ap);                            // And Converts Symbols To Actual Numbers
        va_end(ap);                                 // Results Are Stored In Text
    }
 
    // Here Is Some Code To Split The Text That We Have Been
    // Given Into A Set Of Lines. 
    // This Could Be Made Much Neater By Using
    // A Regular Expression Library Such As The One Available From
    // boost.org (I've Only Done It Out By Hand To Avoid Complicating
    // This Tutorial With Unnecessary Library Dependencies).
    const char *start_line=text;
    vector<string> lines;
    const char *c;
    for(c=text;*c;c++) {
        if(*c=='\n') {
            string line;
            for(const char *n=start_line;n<c;n++) line.append(1,*n);
            lines.push_back(line);
            start_line=c+1;
        }
    }
    if(start_line) {
        string line;
        for(const char *n=start_line;n<c;n++) line.append(1,*n);
        lines.push_back(line);
    }
 
    glPushAttrib(GL_LIST_BIT | GL_CURRENT_BIT  | GL_ENABLE_BIT | GL_TRANSFORM_BIT);
    glMatrixMode(GL_MODELVIEW);
    glDisable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     
 
    glListBase(font);

    float modelview_matrix[16];    
    glGetFloatv(GL_MODELVIEW_MATRIX, modelview_matrix);
 
    // This Is Where The Text Display Actually Happens.
    // For Each Line Of Text We Reset The Modelview Matrix
    // So That The Line's Text Will Start In The Correct Position.
    // Notice That We Need To Reset The Matrix, Rather Than Just Translating
    // Down By h. This Is Because When Each Character Is
    // Drawn It Modifies The Current Matrix So That The Next Character
    // Will Be Drawn Immediately After It. 
    for(int i=0;i<lines.size();i++) {
        glPushMatrix();
        glLoadIdentity();
        glTranslatef(x,y-h*i,0);
	  //glRasterPos2f ((float)x, (float)y);
        glMultMatrixf(modelview_matrix);
 
    // The Commented Out Raster Position Stuff Can Be Useful If You Need To
    // Know The Length Of The Text That You Are Creating.
    // If You Decide To Use It Make Sure To Also Uncomment The glBitmap Command
    // In make_dlist().
        // glRasterPos2f(0,0);
        glCallLists(lines[i].length(), GL_UNSIGNED_BYTE, lines[i].c_str());
        //float rpos[4];
        //glGetFloatv(GL_CURRENT_RASTER_POSITION ,rpos);
        // float len=x-rpos[0]; (Assuming No Rotations Have Happend)
 
        glPopMatrix();
	break;
    }
 
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glPopAttrib();         
 
    pop_projection_matrix();
}
static void *
initialize_hidden(GLfloat width, GLfloat height)
{
	OSMesaContext ctx;
	void *buffer;

	/* Create an RGB-mode context */
	ctx = OSMesaCreateContext( GL_RGBA, NULL );

	/* Allocate the image buffer */
	buffer = malloc( width * height * 4 );

	/* Bind the buffer to the context and make it current */
	OSMesaMakeCurrent( ctx, buffer, GL_UNSIGNED_BYTE, width, height );

	glClearColor(1., 1., 1., 1.0);
	glClearDepth(1.0);
	glDisable(GL_DEPTH_TEST);

	glClear( GL_COLOR_BUFFER_BIT );
	glViewport(0, 0, width, height);
	glMatrixMode( GL_PROJECTION );
	 glDisable ( GL_LIGHTING ) ;
	glLoadIdentity();
	  glOrtho(-60.0, 60.0,-60., 60.,  1.0, -1.0);
	glMatrixMode(GL_MODELVIEW);
	glFlush();
	return buffer;
}
#endif

void abort_(const char * s, ...)
{
        va_list args;
        va_start(args, s);
        vfprintf(stderr, s, args);
        fprintf(stderr, "\n");
        va_end(args);
        abort();
}
static GLfloat* pColors[MAX_PLAQ_COLORS+1];
static int simLock = 0;
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

void resetColorMask(int *colorMask, int regions, int plaqColors[])
{
    *colorMask = 0;
    for (int i=0; i < regions; i++)
	if (plaqColors[i] > 0)
	    *colorMask |= 1 << (plaqColors[i] -1);
}

int nextColor(int *colorMask, int *last)
{
    int retColor;
    for (int i=0; i < MAX_PLAQ_COLORS-1; i++)
	if ((*colorMask & (1 << i)) == 0)
	{
	    *colorMask |= 1 << i;
	    *last = i+1;
	    return i+1;
	}
    // if all are taken, then return one beyond last one or 1st*/
    retColor= (*last<MAX_PLAQ_COLORS-1)?(*last+1):1;
    *last = retColor;

    return retColor;
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
bool model5 (int inT0, double time, 
	state_var* Ec[], state_var* Ip[], 
	state_var* Tp[], gsl_matrix *pastTps, 
	state_var* Vi[], state_var*  Ve[], 
	state_var*  Veadj[], unsigned int *Iplaquebirths, 
	int *totalPlaques, int plaqColors[], 
	gsl_vector *params, globalState *vars, 
	double Repro[], double logRepro[], double p, double latent_inf, int *lastColor, 
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
	unsigned long int *newVirons, 
	unsigned long int *IthisTot, unsigned long int *IdthisTot, unsigned long int *VithisTot, 
	unsigned long int *VethisTot, int selectedRegions[]) 
{
        double tstep;   
	double db, beta, beta_un, hill, fpos,  theta, FIp, FInf, delta, c, 
		rho, r, an, inf, rinf, eclipse, cd8_ic50, alpha, kappa, expand_days; 
	double decline;

	unsigned long int Inf_e[MAX_HEXCELLS], Inf_neu;

	long int Vedelta;
	unsigned long int Vesaved;

	long int Videlta;
	unsigned long int Visaved;

	unsigned long int Ve_Death, Ve_Inf, Tcell_Death;

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

	/* set parameters using passed param vector */

	/* NOTE: latent_inf and p passed in to allow for drug effects (in models 6 & 7) */

	beta=gsl_vector_get(params,0);
	c=gsl_vector_get(params,3);
	theta=gsl_vector_get(params,4);
	inf=gsl_vector_get(params,5);
	r=gsl_vector_get(params,6);
	rinf=gsl_vector_get(params,7);
	delta=gsl_vector_get(params,8);
	rho=gsl_vector_get(params,10);
	eclipse=gsl_vector_get(params,11);
	beta_un=pow(10.,gsl_vector_get(params,12));

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

	int colorMask=0;

	// clear color mask based on colors in use
	resetColorMask(&colorMask,vars->Regions, plaqColors);
	    
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
		    "beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
		    beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse);
		exit(1);
	    }
	    if (theState->Verbose > 1 && Vetot > 4e9 && Velast < 4e9)
	    {
		fprintf(stderr,"Vetot Over 4 billion (%lu) at t=%lf! (was %lu)\n",Vetot, time, Velast);
		fprintf(stderr,
		    "beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
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
		    "beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
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
		    "beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf\n",
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
	*newVirons=0; 

	for (int i =0; i < vars->Regions; i++)
	{
	    Sinf[i] = (stoch == 1) ? gsl_ran_binomial (vars->ur, beta*tstep, Ip[i]->val()): (beta*tstep* Ip[i]->val());

	    Sinftot += Sinf[i];
        }

	Inf_neu = (stoch == 1) ? gsl_ran_poisson (vars->ur, latent_inf*tstep):latent_inf*tstep;

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

		if (plaqColors[index] == 0)
		    plaqColors[index] = nextColor(&colorMask,lastColor);

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
		if (plaqColors[i] == 0)
		    plaqColors[i] = nextColor(&colorMask,lastColor);
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
			"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf\n",
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
	unsigned int AdjDeath, AdjInf;

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
		    "beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf\n",
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
			"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf kappa=%lf alpha=%lf cd8_ic50=%lf fpos=%lf\n",
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

int zoomin(double (*F)( int *,gsl_vector *,void *,FILE *),gsl_rng *r,
    gsl_vector *params,double *Max,
    gsl_vector *low_bound,gsl_vector *high_bound,
    int param_mask, int max_params,
    int max_steps,int stop_walk,int bvstop_walk,
    int print, double tolerance,void *data, int threaded, int searchOrder);

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

    SET_PARAMETER("TSim",TSim);
    SET_PARAMETER("T0_Sim",T0_Sim);
    SET_OPT_PARAMETER("T2_Sim",T2_Sim);
    SET_OPT_PARAMETER("T3_Sim",T3_Sim);
    SET_PARAMETER("tstep",tstep);
    SET_PARAMETER("Tmindel2",Tmindel2);

    SET_PARAMETER("Io",Io);
    SET_PARAMETER("To",To);

    SET_OPT_PARAMETER("alpha_init",alpha_init);
    SET_OPT_PARAMETER("alpha_low",alpha_low);
    SET_OPT_PARAMETER("alpha_high",alpha_high);
    SET_OPT_PARAMETER("alpha_mean",alpha_mean);
    SET_OPT_PARAMETER("alpha_std",alpha_std);

    SET_PARAMETER("beta_init",beta_init);
    SET_PARAMETER("beta_low",beta_low);
    SET_PARAMETER("beta_high",beta_high);
    SET_PARAMETER("log_p_init",log_p_init);
    SET_PARAMETER("log_p_low",log_p_low);
    SET_PARAMETER("log_p_high",log_p_high);
    SET_PARAMETER("latent_inf_init",latent_inf_init);
    SET_PARAMETER("latent_inf_low",latent_inf_low);
    SET_PARAMETER("latent_inf_high",latent_inf_high);
    SET_PARAMETER("theta_init",theta_init);
    SET_PARAMETER("theta_low",theta_low);
    SET_PARAMETER("theta_high",theta_high);
    SET_PARAMETER("r_init",r_init);
    SET_PARAMETER("r_low",r_low);
    SET_PARAMETER("r_high",r_high);
    SET_PARAMETER("delta_init",delta_init);
    SET_PARAMETER("delta_low",delta_low);
    SET_PARAMETER("delta_high",delta_high);
    SET_PARAMETER("c_init",c_init);
    SET_PARAMETER("c_low",c_low);
    SET_PARAMETER("c_high",c_high);
    SET_PARAMETER("inf_init",inf_init);
    SET_PARAMETER("inf_low",inf_low);
    SET_PARAMETER("inf_high",inf_high);
    SET_PARAMETER("inf_mean",inf_mean);
    SET_PARAMETER("inf_std",inf_std);
    SET_PARAMETER("rinf_init",rinf_init);
    SET_PARAMETER("rinf_low",rinf_low);
    SET_PARAMETER("rinf_high",rinf_high);
    SET_PARAMETER("rinf_mean",rinf_mean);
    SET_PARAMETER("rinf_std",rinf_std);
    SET_PARAMETER("rho_init",rho_init);
    SET_PARAMETER("rho_low",rho_low);
    SET_PARAMETER("rho_high",rho_high);
    SET_PARAMETER("db",db);
    SET_PARAMETER("To_mean",To_mean);
    SET_PARAMETER("To_std",To_std);
    SET_PARAMETER("an",an);
    SET_PARAMETER("fpos",fpos);
    SET_OPT_INT_PARAMETER("density_killing",density_killing);
    SET_OPT_PARAMETER("hill",hill);
    SET_PARAMETER("swabInterval",swabInterval);
    SET_OPT_PARAMETER("statInterval",statInterval);
    SET_OPT_PARAMETER("Tolerance",Tolerance);
    SET_PARAMETER("refresh",refresh);

    SET_OPT_INT_PARAMETER("del1",del1);
    SET_OPT_INT_PARAMETER("Tdelay_on",Tdelay_on);
    SET_INT_PARAMETER("Tmin",Tmin);
    SET_INT_PARAMETER("N_runs",N_runs);
    SET_OPT_INT_PARAMETER("Max_steps",Max_steps);
    SET_OPT_INT_PARAMETER("Param_mask",Param_mask);
    SET_OPT_INT_PARAMETER("Stop_walk",Stop_walk);
    SET_OPT_INT_PARAMETER("Bvstop_walk",Bvstop_walk);
    SET_OPT_INT_PARAMETER("Printmax",Printmax);
    SET_OPT_INT_PARAMETER("Threading",Threading);

    SET_OPT_INT_PARAMETER("Fit_model",Fit_model);
    SET_OPT_INT_PARAMETER("Rand_start",Rand_start);
    SET_INT_PARAMETER("Calc_T0",Calc_T0);
    SET_INT_PARAMETER("writeOn",writeOn);
    SET_INT_PARAMETER("Regions",Regions);
    SET_INT_PARAMETER("Crit_mask",Crit_mask);
    SET_INT_PARAMETER("Match_strategy",Match_strategy);
    SET_OPT_INT_PARAMETER("Search_order",Search_order);
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
	SET_PARAMETER("log_betaun_init",log_betaun_init);
	SET_PARAMETER("log_betaun_low",log_betaun_low);
	SET_PARAMETER("log_betaun_high",log_betaun_high);
	SET_PARAMETER("eclipse_init",eclipse_init);
	SET_PARAMETER("eclipse_low",eclipse_low);
	SET_PARAMETER("eclipse_high",eclipse_high);
    }
    if (vars->Model >= 5) 
    {
	SET_OPT_INT_PARAMETER("Model_0",Model_0);
	SET_PARAMETER("gamma_init",gamma_init);
	SET_PARAMETER("gamma_high",gamma_high);
	SET_PARAMETER("gamma_low",gamma_low);
	SET_PARAMETER("gamma_mean",gamma_mean);
	SET_PARAMETER("Cmax_init",Cmax_init);
	SET_PARAMETER("Cmax_high",Cmax_high);
	SET_PARAMETER("Cmax_low",Cmax_low);
	SET_PARAMETER("Cmax_mean",Cmax_mean);
	SET_OPT_PARAMETER("Cmax_0",Cmax_0);
	SET_PARAMETER("IC50_init",IC50_init);
	SET_PARAMETER("IC50_high",IC50_high);
	SET_PARAMETER("IC50_low",IC50_low);
	SET_PARAMETER("IC50_mean",IC50_mean);
	SET_PARAMETER("m_init",m_init);
	SET_PARAMETER("m_high",m_high);
	SET_PARAMETER("m_low",m_low);
	SET_PARAMETER("m_mean",m_mean);
	SET_PARAMETER("bolus",bolus);
	SET_PARAMETER("absorb_init",absorb_init);
	SET_PARAMETER("absorb_high",absorb_high);
	SET_PARAMETER("absorb_low",absorb_low);
	SET_PARAMETER("absorb_mean",absorb_mean);

	SET_OPT_INT_PARAMETER("absorb_hrs",absorb_hrs);
	SET_OPT_INT_PARAMETER("gamma_hrs",gamma_hrs);

	SET_PARAMETER("beta_mean",beta_mean);
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
	SET_PARAMETER("cd8_ic50_low",cd8_ic50_low);
	SET_PARAMETER("cd8_ic50_high",cd8_ic50_high);
	SET_PARAMETER("cd8_ic50_std",cd8_ic50_std);

	SET_OPT_PARAMETER("exp_days_init",exp_days_init);
	SET_OPT_PARAMETER("exp_days_mean",exp_days_mean);
	SET_OPT_PARAMETER("exp_days_std",exp_days_std);

	SET_PARAMETER("kappa_init",kappa_init);
	SET_PARAMETER("kappa_low",kappa_low);
	SET_PARAMETER("kappa_high",kappa_high);
	SET_PARAMETER("kappa_mean",kappa_mean);
	SET_PARAMETER("kappa_std",kappa_std);

	SET_PARAMETER("beta_std",beta_std);
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
    SET_OPT_INT_PARAMETER("AutoSnapshot",AutoSnapshot);
    SET_OPT_PARAMETER("SnapshotInterval",SnapshotInterval);
    SET_OPT_PARAMETER("plot_span",plot_span);
    SET_OPT_INT_PARAMETER("plotStyle1",plotStyle1);
    SET_OPT_INT_PARAMETER("plotOpt1",plotOpt1);
    SET_OPT_INT_PARAMETER("plotOpt2",plotOpt2);
    SET_OPT_INT_PARAMETER("plotOpt3",plotOpt3);
    SET_OPT_INT_PARAMETER("plotOpt4",plotOpt4);
    SET_OPT_INT_PARAMETER("plotOpt5",plotOpt5);
    SET_OPT_INT_PARAMETER("plotOpt6",plotOpt6);

    SET_OPT_INT_PARAMETER("plotVe",plotVe);
    SET_OPT_INT_PARAMETER("plotVi",plotVi);
    SET_OPT_INT_PARAMETER("plotInf",plotInf);
    SET_OPT_INT_PARAMETER("plotCd8",plotCd8);
    SET_OPT_INT_PARAMETER("plotACV",plotACV);

    SET_OPT_INT_PARAMETER("plotRegions",plotRegions);
    SET_OPT_INT_PARAMETER("plotColor",plotColor);
    SET_OPT_INT_PARAMETER("plotLogs",plotLogs);
}
#ifndef NO_GUI
void take_screenshot( GtkWidget *ok_button, char *filename );
#endif

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
double ScoreFunction(int *valid, gsl_vector *ParamVector,
	void *data, FILE *results)
{ 
    double score =0.;

    double score1=0.,score2=0.,score3=0.;

    *valid=0;

    globalState *vars = (globalState *)data;

    static int printCount = vars->Printmax;

    int Nepis = 10000;

    double shed_thresh=100.; /* level of detection for virus shedding episode */

    int model;

    double Size_limit, Episode_limit, T1Sim, T0_Sim;
    double TSim, tstep, timep, period_p, refreshp, nextsnap=0;

    int snapnum=0;
    int snap_this_frame=0;

    double swabT=0., max_T=0., First_T=0., Last_T=0.,cont_First_T=0.;
    bool inContEpisode=false;

    double start_inf_time[MAX_HEXCELLS];
    double start_R0[MAX_HEXCELLS];

    double totPlaques[Nepis]; 

    double episodDur[Nepis], 
	    First_VL[Nepis], 
	    Last_VL[Nepis], 
	    maxFirstReg[Nepis];

    double cont_episodDur[Nepis]; 
    double period_episodDur[Nepis]; 

    double Repro[MAX_HEXCELLS];
    double logRepro[MAX_HEXCELLS];
    double ReproMax[MAX_HEXCELLS];
    double ReproMin[MAX_HEXCELLS];
	    
    double ViOnset[MAX_HEXCELLS];
    double ViLifeMax[MAX_HEXCELLS];
    double ViLifeMin[MAX_HEXCELLS];
    double ViLifeSpan;

    int swabs[VL_BINS];
    int criteria1[VL_BINS];
    double crit1perc[VL_BINS];

    double minCrit1Perc=100.;
    double maxCrit1Perc=0.;
    double avgCrit1Perc=0.;

    int peaks[EPI_BINS];
    int criteria2[EPI_BINS];
    double crit2perc[EPI_BINS];

    int dur_cats[DUR_BINS];
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
    double timeAtGlobalStart=-1;


    int VL_bin, epiD_bin, maxVL_bin;

    double beta, theta, beta_un, delta, c, r, latent_inf, p, avgVL; 
    double rho, inf, rinf, eclipse;

    unsigned long int pastVL;
    unsigned long int currentVL=0;
    unsigned long int measuredVL;
    unsigned long int past_measuredVL;

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
    int period_totalEpis;
    int totalSwabs;
    int totalPosSwabs;
    int posSwabs;
    int withPosSwabs=0;
    int thisSubjSwabs=0;
    int shedderSwabs=0;

    double patSwabMeds[Nepis];	// Track median pos swab VL for each patient (report avg)
    double patSwabVars[Nepis];	// Track variance of swabs for each patient (report avg)

    double totSwabs[MAX_SWABS];
    double totPosSwabs[MAX_SWABS];
    double totPeaks[MAX_SWABS];
    double totDurations[MAX_SWABS];
    double thisPosSwabs[MAX_SWABS];

    double avgEpisYr;

    double time_in_days;
    double total_stats_time = 0;
    double stats_time = 0;

    unsigned int Iplaquebirths =0; /* used for model 5 plaques */
    int plaque_births =0; /* used for model 5 plaques */

    /* model 1 and 5 state variables */
    state_var *Ip[MAX_HEXCELLS], *Tp[MAX_HEXCELLS];
    state_var *pre_Tp[MAX_HEXCELLS];

    /* model 5 only state variables */
    state_var *Ve[MAX_HEXCELLS];
    state_var *Ec[MAX_HEXCELLS], *Veadj[MAX_HEXCELLS], *Vi[MAX_HEXCELLS];
    state_var *Ithis[MAX_HEXCELLS], *Vithis[MAX_HEXCELLS], *Vethis[MAX_HEXCELLS];
    state_var *Id_this[MAX_HEXCELLS];

    unsigned long int Iun_p = 0;

    double *actArray;

    int transmissions=0;
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

    double swab_over_1 = 0;
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

    int lastColor; /* tracks last color used when color tiling plaques */
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
    unsigned long int newInfcells = 0;	/* newly infected cells */
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
    Size_limit = vars->Size_limit; 
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
    beta=gsl_vector_get(ParamVector,0);
    latent_inf=gsl_vector_get(ParamVector,1);
    p=pow(10.,gsl_vector_get(ParamVector,2));
    c=gsl_vector_get(ParamVector,3);
    theta=gsl_vector_get(ParamVector,4);
    inf=gsl_vector_get(ParamVector,5);
    r=gsl_vector_get(ParamVector,6);
    rinf=gsl_vector_get(ParamVector,7);
    delta=gsl_vector_get(ParamVector,8);
    rho=gsl_vector_get(ParamVector,10);
    eclipse=gsl_vector_get(ParamVector,11);
    beta_un=pow(10.,gsl_vector_get(ParamVector,12));

    if (vars->Model >= 5) /* get drug related params */
    {
	Cmax=gsl_vector_get(ParamVector,13);
	IC50=gsl_vector_get(ParamVector,14);
	m=gsl_vector_get(ParamVector,15);
	if (vars->gamma_hrs)
	    gamma=1.0 / (gsl_vector_get(ParamVector,16)/ 24.);
	else
	    gamma=gsl_vector_get(ParamVector,16); 

	if (vars->absorb_hrs)
	    absorb=(gsl_vector_get(ParamVector,17)/ 24.);
	else
	    absorb=gsl_vector_get(ParamVector,17); 

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
	
    for (int j=0; j <vars->Regions && j < MAX_HEXCELLS; j++)
    {
	// initialize all state variables
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

	    /* reset parameter vectors incase any changed */
	    gsl_vector_set(ParamVector,0,vars->beta_init); 
	    gsl_vector_set(ParamVector,1,vars->latent_inf_init); 
	    gsl_vector_set(ParamVector,2,vars->log_p_init); 
	    gsl_vector_set(ParamVector,3,vars->c_init); 
	    gsl_vector_set(ParamVector,4,vars->theta_init); 
	    gsl_vector_set(ParamVector,5,vars->inf_init); 
	    gsl_vector_set(ParamVector,6,vars->r_init); 
	    gsl_vector_set(ParamVector,7,vars->rinf_init); 
	    gsl_vector_set(ParamVector,8,vars->delta_init); 
	    gsl_vector_set(ParamVector,10,vars->rho_init); 
	    gsl_vector_set(ParamVector,11,vars->eclipse_init); 


	    /* set parameters using passed param vector */
	    beta=gsl_vector_get(ParamVector,0);
	    latent_inf=gsl_vector_get(ParamVector,1);
	    p=pow(10.,gsl_vector_get(ParamVector,2));
	    c=gsl_vector_get(ParamVector,3);
	    theta=gsl_vector_get(ParamVector,4);
	    inf=gsl_vector_get(ParamVector,5);
	    r=gsl_vector_get(ParamVector,6);
	    rinf=gsl_vector_get(ParamVector,7);
	    delta=gsl_vector_get(ParamVector,8);
	    rho=gsl_vector_get(ParamVector,10);
	    eclipse=gsl_vector_get(ParamVector,11);

	    if (vars->Model > 5) /* adjust drug related params */
	    {
		gsl_vector_set(ParamVector,13,vars->Cmax_init); 
		gsl_vector_set(ParamVector,14,vars->IC50_init); 
		gsl_vector_set(ParamVector,15,vars->m_init); 
		gsl_vector_set(ParamVector,16,vars->gamma_init); 
		gsl_vector_set(ParamVector,17,vars->absorb_init); 

		Cmax=gsl_vector_get(ParamVector,13);
		IC50=gsl_vector_get(ParamVector,14);
		m=gsl_vector_get(ParamVector,15);
		if (vars->gamma_hrs)
		    gamma=1.0 / (gsl_vector_get(ParamVector,16) / 24.);
		else
		    gamma=gsl_vector_get(ParamVector,16); 

		if (vars->absorb_hrs)
		    absorb=(gsl_vector_get(ParamVector,17)/ 24.);
		else
		    absorb=gsl_vector_get(ParamVector,17); 

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
	    if (vars->Calc_T0 == -1)
		sprintf(t0_file, "%s.%x.%d",t0_default,pthread_self(),randfile);
	    else
		sprintf(t0_file, "%s.%d",t0_default,randfile);
	}

	if ((!batchMode || vars->Calc_T0 <= 0 || vars->Calc_T0==3) && 
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
	    fprintf(stderr,"Total T0 over %d regions read in is %u\n",vars->Regions,Tptot);
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
		    fprintf(stderr,"Total T0 over %d regions randomly generated is %u\n",vars->Regions,Tptot);
	    }
	    else
		Tp[0]->initialz_with(vars->To);

	    if (vars->Calc_T0 != 0)
	    {
		// initialize all state variables
		for (int j=0; j <vars->Regions; j++)
		{
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
		    if (vars->points != NULL)
		    {
			vars->vet[j][vars->sample_index] = 0;
			vars->vit[j][vars->sample_index] = 0;
			vars->inf[j][vars->sample_index]=Ip[j]->val();
			vars->cd8[j][vars->sample_index]=Tp[j]->val();
			vars->repro[j][vars->sample_index]=0;
			vars->color[j][vars->sample_index] = 0;
		    }
		}
		time_in_days=0.; 
		stats_time=0.; 
		swab_over_1=0.; 
		time_over_1=0.; 
		time_over_2=0.; 
		swabT=0.;
		lastColor=0;
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
		newVirons=0;
		Itot_p = 0;
		Iun_p = 0;

		if (vars->Model >= 4)
		{
		    ACV = 0;
		    lastBolus = -vars->bolus;
		    firstBolus = 1;
		    infuse = 0;
		    doses=0;
		}
    #ifndef NO_GUI
    #ifndef __sun
		if (!batchMode)
		{
		    GdkWindow *win = gtk_widget_get_window( main_window );
		    wait_cursor(win);
		}
    #endif
    #endif
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
			    "To=%lf, beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf absorb=%lf gamma=%lf Cmax=%lf IC50=%lf m=%lf regions=%d swabInt=%lf\n",
			    vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,absorb,gamma,Cmax,IC50,m,vars->Regions,vars->swabInterval);
		    else
			fprintf(stderr,
			    "To=%lf, beta=%lf latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf time=%lf\n",
			    vars->To,beta,latent_inf,vars->log_p_init,c,theta,delta,r,inf,rinf,rho,eclipse,vars->fpos,vars->an,time_in_days);
		}

		if (vars->Calc_T0 < 2)
		{
		    while (time_in_days<= T0_Sim) 
		    {
			if (vars->Model_0 == 4 || vars->Model_0 == 5)
			{
			    if (!model5(1,time_in_days,Ec, Ip, 
				    Tp, pastTps, Vi, Ve, 
				    Veadj, &Iplaquebirths, &totalPlaques,
				    plaqColors, ParamVector, vars, 
				    Repro, logRepro, p, latent_inf,&lastColor,
				    Vethis, Vithis, Ithis, Id_this, start_inf_time,
				    &repro_max, &repro_mean, &log_repro_mean,
				    &repro_std, &log_repro_std,
				    &Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newVirons,
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

			    if (!model5(1,time_in_days,Ec, Ip, 
				    Tp, pastTps, Vi, Ve, Veadj, &Iplaquebirths, 
				    &totalPlaques, plaqColors, ParamVector, vars, 
				    Repro,logRepro,p_ACV,latent_inf_ACV,&lastColor,
				    Vethis, Vithis, Ithis, Id_this, start_inf_time,
				    &repro_max, &repro_mean, &log_repro_mean,
				    &repro_std, &log_repro_std,
				    &Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newVirons,
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

			    if (!model5(1,time_in_days,Ec, Ip, 
				    Tp, pastTps, Vi, Ve, Veadj, &Iplaquebirths, 
				    &totalPlaques, plaqColors, ParamVector, vars, 
				    Repro,logRepro,p_ACV,latent_inf_ACV,&lastColor,
				    Vethis, Vithis, Ithis, Id_this, start_inf_time,
				    &repro_max, &repro_mean, &log_repro_mean,
				    &repro_std, &log_repro_std,
				    &Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newVirons,
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
    #ifndef NO_GUI
			if (!batchMode)
			{
			  if (gtk_events_pending ())
			    gtk_main_iteration();
			}
    #endif
		    }
		    if (vars->Verbose > 1)
			fprintf(stderr,
			    "[%d]: Calculated T0=%lu (after %lf days)\n",
			    Nr_count,T0,time_in_days);
    #ifndef NO_GUI
	#ifndef __sun
		    if (!batchMode)
		    {
			GdkWindow *win = gtk_widget_get_window( main_window );
			normal_cursor(win);
		    }
	#endif
    #endif
		    if (vars->Fit_model == 0 || vars->Calc_T0 < 0)
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
				fprintf(stderr,"(thread : wrote T0 file %s (elapsed=%ld)\n",
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

	// initialize all state variables
	// (Leave Ve, Ip and Tp at burn-in levels)
	for (int j=0; j <vars->Regions; j++)
	{
	    Ec[j]->initialz_with(0); 
	    Ip[j]->initialz_with(0); 
	    Ve[j]->initialz_with(0); 
	    Vi[j]->initialz_with(0); 
	    Vethis[j]->initialz_with(Ve[j]->val()); 
	    Vithis[j]->initialz_with(Vi[j]->val()); 
	    Ithis[j]->initialz_with(Ip[j]->val()); 
	    Id_this[j]->initialz_with(0); 
	    Veadj[j]->initialz_with(0); 

	    if (vars->points != NULL)
	    {
		vars->vet[j][vars->sample_index] = Ve[j]->val();
		vars->vit[j][vars->sample_index] = Vi[j]->val();
		vars->inf[j][vars->sample_index]=Ip[j]->val();
		vars->cd8[j][vars->sample_index]=Tp[j]->val();
		vars->repro[j][vars->sample_index]=0;
		vars->color[j][vars->sample_index]=0;
	    }
	}
	ACV = 0;
	lastColor=0;
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
	newVirons=0;
	cd8size=0;
	vet=0;
	vet_p=0;
	vit=0;
	totalPlaques=0;
	infTot=0;
	infGlobalMax=0;
	timeAtGlobalStart=0;
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
	refreshp=0.; 
	for(int i=0;i<Nepis;i++)
	{
	    maxVL[i] = 0;
	    maxFirstReg[i] = 0.;
	    First_VL[i] = 0.;
	    Last_VL[i] = 0.;
	    episodDur[i] = 0.;
	    cont_episodDur[i] = 0.;
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
	past_measuredVL = 0.;
	swabsThisEpi=0;

	bool inEpisode=false;
	bool inEpisode_p=false;
	bool any_neg=false;
	double max_dur = 0;
	Iplaquebirths=0;
	plaque_births=0;
	avgVL=0.;

	for (int j=0; j <vars->Regions; j++) {
	    ReproMax[j] = 0;
	    ReproMin[j] = 100.;
	    ViOnset[j] = 0;
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
	    if (vars->dataF11 != NULL) 
	    {
		fprintf(vars->dataF11,"episode,");
		fprintf(vars->dataF11,"duration,");
		fprintf(vars->dataF11,"overall peak,");
		fprintf(vars->dataF11,"plaques,");
		fprintf(vars->dataF11,"first region peak\n");
	    }
	    if (vars->dataF12 != NULL) 
	    {
		fprintf(vars->dataF12,"time,");
		fprintf(vars->dataF12,"T_max");
		fprintf(vars->dataF12,",T_tot");
		fprintf(vars->dataF12,",Ve_max");
		fprintf(vars->dataF12,",Ve_tot");
		fprintf(vars->dataF12,",R_max");
		fprintf(vars->dataF12,",R_mean");
		fprintf(vars->dataF12,",R_stddev");
		fprintf(vars->dataF12,"\n");
	    }
	    if (vars->dataF13 != NULL) 
	    {
		fprintf(vars->dataF13,"time");
		fprintf(vars->dataF13,",model");
		fprintf(vars->dataF13,",bolus");
		fprintf(vars->dataF13,",doses");
		fprintf(vars->dataF13,",ACV");
		fprintf(vars->dataF13,",Infus");
		fprintf(vars->dataF13,",lastBolus");
		fprintf(vars->dataF13,",p");
		fprintf(vars->dataF13,",p_ACV");
		fprintf(vars->dataF13,",latent_inf");
		fprintf(vars->dataF13,",latent_inf_ACV");
		fprintf(vars->dataF13,",vet");
		fprintf(vars->dataF13,"\n");
	    }
	}

	int sig_episodes=0;
	infTot=0;
	infGlobalMax=0;
	bool episodeAtStart = false;

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
	T_to_date = 0;
	T_inf_tot = 0;
	T_bump_tot = 0;
	T_dump_tot = 0;
	T_kills_tot = 0;
	Ideaths_tot = 0;

	while (((Total_epis > 0 && totalEpis + counter < Total_epis) || 
			Total_epis==0) &&
		(((Episode_limit > 0 && counter < Episode_limit) || 
		  (Episode_limit == 0 && time_in_days<=TSim )) && counter<Nepis && !episodeStop && !vars->stopFlag))
	{
	    if (vars->alt_inp_file != NULL &&
		vars->Input_refresh > 0 && time_in_days>=vars->Input_refresh &&
		time_in_days<=vars->Input_refresh+tstep)
	    {
		fprintf(results,"Reading in %s at t=%lf\n",
			vars->alt_inp_file,time_in_days);
		if (vars->Model > 5)
		    fprintf(results,
			"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 ,absorb=%lf, gamma=%lf,Cmax=%lf,IC50=%lf,m=%lf,time=%lf\n",
			beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11,absorb,gamma,Cmax,IC50,m,stats_time);
		else
		    fprintf(results,
			"To=%lf,beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 ,time=%lf\n",
			vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11,stats_time);

		read_input_file(1, vars->alt_inp_file,vars);

		// if we are to use mean/stddev then redo parameters after file load in case they changed
		// (no change allowed if PDF_on==1)
		if (vars->PDF_on==2)
		{
		    do
			vars->beta_init = vars->beta_mean+gsl_ran_gaussian (vars->ur, vars->beta_std);
		    while (vars->beta_init < 0);
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
		gsl_vector_set(ParamVector,0,vars->beta_init); 
		gsl_vector_set(ParamVector,1,vars->latent_inf_init); 
		gsl_vector_set(ParamVector,2,vars->log_p_init); 
		gsl_vector_set(ParamVector,3,vars->c_init); 
		gsl_vector_set(ParamVector,4,vars->theta_init); 
		gsl_vector_set(ParamVector,5,vars->inf_init); 
		gsl_vector_set(ParamVector,6,vars->r_init); 
		gsl_vector_set(ParamVector,7,vars->rinf_init); 
		gsl_vector_set(ParamVector,8,vars->delta_init); 
		gsl_vector_set(ParamVector,10,vars->rho_init); 
		gsl_vector_set(ParamVector,11,vars->eclipse_init); 
		gsl_vector_set(ParamVector,12,vars->log_betaun_init); 
		/* set parameters using passed param vector */
		beta=gsl_vector_get(ParamVector,0);
		latent_inf=gsl_vector_get(ParamVector,1);
		p = pow(10,gsl_vector_get(ParamVector,2));
		c=gsl_vector_get(ParamVector,3);
		theta=gsl_vector_get(ParamVector,4);
		inf=gsl_vector_get(ParamVector,5);
		r=gsl_vector_get(ParamVector,6);
		rinf=gsl_vector_get(ParamVector,7);
		delta=gsl_vector_get(ParamVector,8);
		rho=gsl_vector_get(ParamVector,10);
		eclipse=gsl_vector_get(ParamVector,11);
		beta_un=pow(10,gsl_vector_get(ParamVector,12));

		if (vars->Model > 5) /* get drug related params */
		{
		    gsl_vector_set(ParamVector,13,vars->Cmax_init); 
		    gsl_vector_set(ParamVector,14,vars->IC50_init); 
		    gsl_vector_set(ParamVector,15,vars->m_init); 
		    gsl_vector_set(ParamVector,16,vars->gamma_init); 
		    gsl_vector_set(ParamVector,17,vars->absorb_init); 

		    Cmax=gsl_vector_get(ParamVector,13);
		    IC50=gsl_vector_get(ParamVector,14);
		    m=gsl_vector_get(ParamVector,15);
		    if (vars->gamma_hrs)
			gamma=1.0 / (gsl_vector_get(ParamVector,16) / 24.);
		    else
			gamma=gsl_vector_get(ParamVector,16); 

		    if (vars->absorb_hrs)
			absorb=(gsl_vector_get(ParamVector,17)/ 24.);
		    else
			absorb=gsl_vector_get(ParamVector,17); 

		    bolus = vars->bolus;
		    if (vars->Cmax_0 == 0.0)
			Cmax_0 = 0.0;
		    else
			Cmax_0 = vars->Cmax_0;
		}
		if (vars->Model > 5)
		    fprintf(results,
			"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 ,absorb=%lf, gamma=%lf,Cmax=%lf,IC50=%lf,m=%lf\n",
			beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11,absorb,gamma,Cmax,IC50,m);
		else
		    fprintf(results,
			"To=%lf,beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 \n",
			vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11);
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

		if (!model5(0,time_in_days,Ec,Ip,Tp,pastTps,  
			Vi, Ve, Veadj, &Iplaquebirths, &totalPlaques, plaqColors, 
			ParamVector, vars, Repro,logRepro,p,
			latent_inf,&lastColor,
			Vethis, Vithis, Ithis, Id_this, start_inf_time,
			&repro_max, &repro_mean, &log_repro_mean,
			&repro_std, &log_repro_std,
			&Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newVirons,
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

		if (!model5(0,time_in_days,Ec,Ip,Tp,pastTps,  
			Vi, Ve, Veadj, &Iplaquebirths, &totalPlaques, plaqColors, 
			ParamVector, vars, Repro,logRepro,p_ACV,
			latent_inf_ACV,&lastColor,
			Vethis, Vithis, Ithis, Id_this, start_inf_time,
			&repro_max, &repro_mean, &log_repro_mean,
			&repro_std, &log_repro_std,
			&Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newVirons,
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

		if (!model5(0,time_in_days,Ec,Ip,Tp,pastTps,  
			Vi, Ve, Veadj, &Iplaquebirths, &totalPlaques, plaqColors, 
			ParamVector, vars, Repro,logRepro,p_ACV,
			latent_inf_ACV,&lastColor,
			Vethis, Vithis, Ithis, Id_this, start_inf_time,
			&repro_max, &repro_mean, &log_repro_mean,
			&repro_std, &log_repro_std,
			&Tinf,&Tbump,&Tdump,&Tkills,&Ideaths,&newTcells, &newInfcells, &newVirons,
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
    #ifdef OLDWAY
		    if (ViOnset[j] > 0 && Vi[j]->val() == 0) {
			if (ViLifeMax[j] < time_in_days - ViOnset[j])
			    ViLifeMax[j] = time_in_days - ViOnset[j];
			if (ViLifeMin[j] > time_in_days - ViOnset[j])
			    ViLifeMin[j] = time_in_days - ViOnset[j];
		    }
		    if (ViOnset[j] == 0 && Vi[j]->val() > 0)
			ViOnset[j] = time_in_days;
		    if ( Vi[j]->val() == 0)
		    {
			ViOnset[j] = 0;
		    }
    #endif

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
#ifndef NO_GUI
				if (!batchMode)
				    stop_cb(main_window,NULL);
				else
#endif
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
		timeAtGlobalStart = time_in_days;
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
		    First_VL[counter] = measuredVL;
		    First_T = time_in_days;

		    inEpisode=true;

		    if (vars->Calc_T0 == 2 && 
			(First_T - 1.5*vars->swabInterval < T0_Sim || First_T - 1.5*tstep <  T0_Sim))
			episodeAtStart = true;

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
			if (vars->dataF11 != NULL) {
			    fprintf(vars->dataF11,"%d,",counter);
			    fprintf(vars->dataF11,"%lf,",episodDur[counter]);
			    fprintf(vars->dataF11,"%lu,",maxVL[counter]);
			    fprintf(vars->dataF11,"%lf,",totPlaques[counter]);
			    fprintf(vars->dataF11,"%lf,",maxFirstReg[counter]);
			    fprintf(vars->dataF11,"\n");
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
			    fprintf(stderr,"(thread : Patient %d hit %d episodes at t=%lf (elapsed=%ld)\n",
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
		past_measuredVL = measuredVL;
		if (vars->CritOn && time_in_days >= vars->crit_start)
		{
		    if (measuredVL > 0)
		    {
			totSwabs[totalSwabs]=log10(measuredVL);
		    }
		    else
			totSwabs[totalSwabs]=0;

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
		cont_episodDur[cont_counter]=time_in_days - cont_First_T;
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
		    if (vars->dataF12 != NULL) {
			unsigned long int tmax = 0;
			unsigned long int ttot = 0;
			unsigned long int vemax = 0;
			unsigned long int vetot = 0;
			fprintf(vars->dataF12,"%lf,",time_in_days);
			for (int j=0; j <vars->Regions; j++)
			{
			    if (Tp[j]->val() > tmax)
				tmax = Tp[j]->val();
			    ttot +=Tp[j]->val();
			    if (Ve[j]->val() > vemax)
				vemax = Ve[j]->val();
			    vetot +=Ve[j]->val();
			}
			fprintf(vars->dataF12,"%lu",tmax);
			fprintf(vars->dataF12,",%lu",ttot);
			fprintf(vars->dataF12,",%lu",vemax);
			fprintf(vars->dataF12,",%lu",vetot);
			fprintf(vars->dataF12,",%lf",repro_max);
			fprintf(vars->dataF12,",%lf",repro_mean);
			fprintf(vars->dataF12,",%lf",repro_std);
			fprintf(vars->dataF12,"\n");
		    }
		    if (vars->dataF13 != NULL) 
		    {
			fprintf(vars->dataF13,"%lf,",time_in_days);
			fprintf(vars->dataF13,"%d,",model);
			fprintf(vars->dataF13,"%lf,",bolus);
			fprintf(vars->dataF13,"%d,",doses);
			fprintf(vars->dataF13,"%lf,",ACV);
			fprintf(vars->dataF13,"%lf,",infuse);
			fprintf(vars->dataF13,"%lf,",lastBolus);
			fprintf(vars->dataF13,"%lf,",p);
			fprintf(vars->dataF13,"%lf,",p_ACV);
			fprintf(vars->dataF13,"%lf,",latent_inf);
			fprintf(vars->dataF13,"%lf,",latent_inf_ACV);
			fprintf(vars->dataF13,"%lu",vet);
			fprintf(vars->dataF13,"\n");
		    }
		}

		vars->time = time_in_days;

		if (vars->points != NULL && abs(model) >= 4) 
	        {
		    if (++vars->sample_index == MAX_SAMPLES)
			vars->sample_index = 0;

		    for (int j=0; j <vars->Regions; j++)
		    {
			vars->vet[j][vars->sample_index] = Ve[j]->val();
			vars->vit[j][vars->sample_index] = Vi[j]->val();
			vars->inf[j][vars->sample_index]=Ip[j]->val();
			vars->cd8[j][vars->sample_index]=Tp[j]->val();
			vars->repro[j][vars->sample_index]=Repro[j];
			vars->color[j][vars->sample_index]=plaqColors[j];
		    }
		    if (vars->points->valid < vars->points->max_points-1) {
			vars->points->valid++;

			vars->points->time[vars->points->valid] = time_in_days;
			vars->points->cd8cells[vars->points->valid] = cd8size;
			vars->points->vet[vars->points->valid] = vet;
			vars->points->vit[vars->points->valid] = vit;
			vars->points->inf[vars->points->valid] = infTot;
			vars->points->ACV[vars->points->valid] = ACV;
			for (int j=0; j <vars->Regions; j++)
			{
			    vars->points->ve[j][vars->points->valid] = Ve[j]->val();
			    vars->points->vi[j][vars->points->valid] = Vi[j]->val();
			    vars->points->color[j][vars->points->valid] = plaqColors[j];
			}

			if (time_in_days > vars->max_time)
			    vars->max_time = 2* vars->max_time;
			if (vet > vars->max_vl || vit > vars->max_vl)
			    vars->max_vl = 2* vars->max_vl;
			if (infTot > vars->max_inf)
			    vars->max_inf = 2* vars->max_inf;
			if (cd8size > vars->max_cd8s)
			    vars->max_cd8s = 2* vars->max_cd8s;
#ifndef NO_GUI
			if (image != NULL)
			{
			    if (time_in_days-refreshp > vars->refresh) { 
				updateGL();
#ifdef OLD_WAY
				snap_movie_frame();
#endif
				refreshp=time_in_days;
			    }
			    if (vars->AutoSnapshot &&  
				((snapnum == 0) ||
				(time_in_days>=nextsnap))) { 
				updateGL();
				snap_this_frame=1;  // done after graph updates!
				nextsnap+=vars->SnapshotInterval;
			    }
			}
			else
#endif
			{
			    if (vars->AutoSnapshot &&  
				((snapnum == 0) ||
				(time_in_days>=nextsnap))) { 
				snap_this_frame=1;  // done after graph updates!
				nextsnap+=vars->SnapshotInterval;
			    }
			}
		    }
		    /* check for need to realloc point array */
		    vars->points->checkForRealloc();
	      }
	      timep = time_in_days;
	  }

	  /* check for graph updates atleast once (and whenever paused) */
#ifndef NO_GUI
	  /* check for graph updates atleast once (and whenever paused) */
	  if (!batchMode)
	  {
	      do {
		  if (gtk_events_pending ())
		    gtk_main_iteration();
		  if (snap_this_frame)
		  {
		    char snapfile[100];
		    sprintf(snapfile,
			"%s/snapshot_%d.png",outDir.c_str(),snapnum++);
		    take_screenshot( NULL, snapfile);
		    snap_this_frame=0;
		  }
	      } while (vars->pauseFlag == 1);
	  }
#else
	  if (snap_this_frame)
	  {
	    char snapfile[100];
	    sprintf(snapfile,
		"%s/snapshot_%d.png",outDir.c_str(),snapnum++);

	    draw_routine(DEFAULT_IMAGE_WIDTH,DEFAULT_IMAGE_HEIGHT);

	    write_png_file(DEFAULT_IMAGE_WIDTH,DEFAULT_IMAGE_HEIGHT,snapfile);
	    snap_this_frame=0;
	  }
#endif
	  time_in_days = time_in_days + tstep;

	  if (vars->CritOn && time_in_days >= vars->crit_start &&
	      (vars->Calc_T0 != 2 || time_in_days > T0_Sim))
	  {
	      total_stats_time += tstep;
	      stats_time += tstep;
	  }
	  inEpisode_p = inEpisode;

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
	    cont_episodDur[cont_counter]=time_in_days - cont_First_T;
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
	fprintf(results, "Run %d completed in %lf days with %g positive swabs (max episode = %lf days)\n", 
		Nr_count, time_in_days, totalPosSwabs/(float)totalSwabs,max_dur);
	fprintf(results,
		"beta=%lf latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf cd8_ic50=%lf kappa=%lf alpha=%lf exp_days=%lf hill=%lf\n",
		vars->beta_init,vars->latent_inf_init,vars->log_p_init,vars->c_init,vars->theta_init,vars->delta_init,vars->r_init,vars->inf_init,vars->rinf_init,vars->rho_init,vars->eclipse_init,vars->fpos,vars->an,vars->cd8_ic50_init,vars->kappa_init,vars->alpha_init,vars->exp_days_init,vars->hill);

	fprintf(results,"Total infect cell births %lu\n",Inf_to_date);
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

	if (vars->AutoSnapshot)
	    break;

    }//end of runs loop

    if (totalEpis == 0)
    {
	if (printCount == vars->Printmax) 
	{
		fprintf(results,"No episodes!\n");
		if (vars->Model > 5)
		    fprintf(results,
			"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 ,absorb=%lf, gamma=%lf,Cmax=%lf,IC50=%lf,m=%lf,time=%lf\n",
			beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11,absorb,gamma,Cmax,IC50,m,stats_time);
		else
		    fprintf(results,
			"To=%lf,beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 ,time=%lf\n",
			vars->To,beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11,stats_time);
	}
	score = 0.;
	minCrit1Perc = 0;
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

    if (printCount == vars->Printmax) fprintf(results,"Total percentage of swabs that were positive= %lf (%d of %d)\n",totCrit1Perc,totalPosSwabs,totalSwabs);
    if (printCount == vars->Printmax) fprintf(results,"Total percentage of shedder swabs that were positive= %lf (%d of %d)\n",subCrit1Perc,totalPosSwabs,shedderSwabs);
    if (printCount == vars->Printmax) fprintf(results,"Rate of episodes with >1mm plaques= %lf\n",epi_1mm* (365./total_stats_time));
    if (printCount == vars->Printmax) fprintf(results,"Percent Time in episodes with >1mm plaques= %lf\n",100.*time_over_1/total_stats_time);
    if (printCount == vars->Printmax) fprintf(results,"Rate of episodes with >2mm plaques= %lf\n",epi_2mm* (365./total_stats_time));
    if (printCount == vars->Printmax) fprintf(results,"Percent Time in episodes with >2mm plaques= %lf\n",100.*time_over_2/total_stats_time);

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
	    if (printCount == vars->Printmax) fprintf(results,"criteria 1[%d]: %lf vs. mean of %lf - err = %lf * %lf * %lf\n",
		j+1,crit1perc[j], vars->crit[critNum].mean,err,weight,vars->critWeight[0]);
	    score1+= err*weight;
	}
	else 
	{
	    if (printCount == vars->Printmax) fprintf(results,"criteria 1[%d]: %lf < %lf < %lf - ",
		j+1,vars->crit[critNum].low, crit1perc[j], vars->crit[critNum].high);
	    
	    if (crit1perc[j] < vars->crit[critNum].low )
	    {
		if (printCount == vars->Printmax) fprintf(results,"no\n");
		score1-= (vars->crit[critNum].low - crit1perc[j])/100.; 
	    }
	    else if (crit1perc[j] > vars->crit[critNum].high)
	    {
		if (printCount == vars->Printmax) fprintf(results,"no\n");
		score1-= (crit1perc[j] - vars->crit[critNum].high)/100.; 
	    }
	    else
	    {
		if (printCount == vars->Printmax) fprintf(results,"yes\n");
		score1+=1.0; 
	    }
	}
	critNum++;
    }
    if (printCount == vars->Printmax) fprintf(results,"Score1 = %lf%s\n",score1,(vars->Crit_mask == 0 || (vars->Crit_mask & 1))?"*":"");

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
	    if (printCount == vars->Printmax) fprintf(results,"criteria 2[%d]: %lf < %lf < %lf - ",
		j+1,vars->crit[critNum].low, crit2perc[j], vars->crit[critNum].high);
	    
	    if (crit2perc[j] < vars->crit[critNum].low )
	    {
		if (printCount == vars->Printmax) fprintf(results,"no\n");
		score2-= (vars->crit[critNum].low - crit2perc[j])/100.; 
	    }
	    else if (crit2perc[j] > vars->crit[critNum].high)
	    {
		if (printCount == vars->Printmax) fprintf(results,"no\n");
		score2-= (crit2perc[j] - vars->crit[critNum].high)/100.; 
	    }
	    else
	    {
		if (printCount == vars->Printmax) fprintf(results,"yes\n");
		score2+=1.0; 
	    }
	}
	critNum++;
    }
    if (printCount == vars->Printmax) fprintf(results,"Score2 = %lf%s\n",score2,(vars->Crit_mask == 0 || (vars->Crit_mask & 0x2))?"*":"");

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

    int crit3start = critNum;
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
	    if (printCount == vars->Printmax) fprintf(results,"criteria 3[%d]: %lf < %lf < %lf - ",
		j+1,vars->crit[critNum].low, crit3perc[j], vars->crit[critNum].high);
	    
	    if (crit3perc[j] < vars->crit[critNum].low )
	    {
		if (printCount == vars->Printmax) fprintf(results,"no\n");
		score3-= (vars->crit[critNum].low - crit3perc[j])/100.; 
	    }
	    else if (crit3perc[j] > vars->crit[critNum].high)
	    {
		if (printCount == vars->Printmax) fprintf(results,"no\n");
		score3-= (crit3perc[j] - vars->crit[critNum].high)/100.; 
	    }
	    else
	    {
		if (printCount == vars->Printmax) fprintf(results,"yes\n");
		score3+=1.0; 
	    }
	}
	critNum++;
    }
    if (printCount == vars->Printmax) fprintf(results,"Score3 = %lf%s\n",score3,(vars->Crit_mask == 0 || (vars->Crit_mask & 0x4))?"*":"");

    if (vars->Crit_mask == 0 || (vars->Crit_mask & 0x4)) score += score3;

    if (critNum != NUM_CRITERIA) 
    {
	fprintf(results,"Internal error: processed %d criteria out of %d in file\n",
	    critNum,NUM_CRITERIA);
	exit(1);
    }
    if (vars->Match_strategy != 0)
	score = 1.0 / score; /* invert so that low actual score (least squares) is desireable (high)! */

    if (printCount == vars->Printmax) {
	if (results != stdout)
	    fprintf(stdout,"Total Score = %lf ",score);

	fprintf(results,"Total Score = %lf ",score);
	if (vars->Model > 5)
	{
	    if (results != stdout)
		fprintf(stdout,
		"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 absorb=%lf gamma=%lf Cmax=%lf IC50=%lf m=%lf time=%lf\n",
		beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11,absorb,gamma,Cmax,IC50,m,total_stats_time);
	    fprintf(results,
		"beta=%lf latent_inf=%lf p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf beta_un=%lfe-11 absorb=%lf gamma=%lf Cmax=%lf IC50=%lf m=%lf time=%lf\n",
		beta,latent_inf,p,c,theta,delta,r,inf,rinf,rho,eclipse,beta_un*1.0e11,absorb,gamma,Cmax,IC50,m,total_stats_time);
	}
	else
	{
	    if (results != stdout)
		fprintf(stdout,
		"beta=%lf latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf cd8_ic50=%lf kappa=%lf alpha=%lf exp_days=%lf hill=%lf time=%lf\n",
		vars->beta_init,vars->latent_inf_init,vars->log_p_init,vars->c_init,vars->theta_init,vars->delta_init,vars->r_init,vars->inf_init,vars->rinf_init,vars->rho_init,vars->eclipse_init,vars->fpos,vars->an,vars->cd8_ic50_init,vars->kappa_init,vars->alpha_init,vars->exp_days_init,vars->hill,total_stats_time);

	    fprintf(results,
		"beta=%lf latent_inf=%lf log_p=%lf c=%lf theta=%lf delta=%lf r=%lf inf=%lf rinf=%lf rho=%lf eclipse=%lf fpos=%lf an=%lf cd8_ic50=%lf kappa=%lf alpha=%lf exp_days=%lf hill=%lf time=%lf\n",
		vars->beta_init,vars->latent_inf_init,vars->log_p_init,vars->c_init,vars->theta_init,vars->delta_init,vars->r_init,vars->inf_init,vars->rinf_init,vars->rho_init,vars->eclipse_init,vars->fpos,vars->an,vars->cd8_ic50_init,vars->kappa_init,vars->alpha_init,vars->exp_days_init,vars->hill,total_stats_time);
	}
	/* reset printCount */
	printCount = 1;
    }
    else
	printCount++;

	fprintf(results,"Final percentage of swabs that were positive= %lf\n",totCrit1Perc);
	fprintf(results,"Shedders percentage of swabs that were positive= %lf\n",subCrit1Perc);

    fflush(results);
    *valid=1;
    gsl_matrix_free(pastTps);
    free(actArray);

    return score;	 /* return score for last patient (if multiple ones) */
}

/* This function runs one simulation.  It is the entry point when 
 * the program is run in the non-batch mode */
int runSimulation ( globalState *vars)
{
    gsl_vector *params;

    int valid=0;

    double score=0;

    params = gsl_vector_alloc(NUM_PARAMS);

    gsl_vector_set(params,0,vars->beta_init); 
    gsl_vector_set(params,1,vars->latent_inf_init); 
    gsl_vector_set(params,2,vars->log_p_init); 
    gsl_vector_set(params,3,vars->c_init); 
    gsl_vector_set(params,4,vars->theta_init); 
    gsl_vector_set(params,5,vars->inf_init); 
    gsl_vector_set(params,6,vars->r_init); 
    gsl_vector_set(params,7,vars->rinf_init); 
    gsl_vector_set(params,8,vars->delta_init); 
    gsl_vector_set(params,10,vars->rho_init); 
    gsl_vector_set(params,11,vars->eclipse_init); 
    gsl_vector_set(params,12,vars->log_betaun_init); 

    if (vars->Model >= 5)
    {
	gsl_vector_set(params,13,vars->Cmax_init); 
	gsl_vector_set(params,14,vars->IC50_init); 
	gsl_vector_set(params,15,vars->m_init); 
	gsl_vector_set(params,16,vars->gamma_init); 
	gsl_vector_set(params,17,vars->absorb_init); 
    }

    score=ScoreFunction(&valid, params, (void *)vars,NULL);

    gsl_vector_free(params);

    return valid;
}

double fitmodel(gsl_vector *params, globalState *vars)
{
    gsl_vector *low_bound,*high_bound;

    double score;

    int valid=0;
    int param_count;

    param_count=NUM_PARAMS;
    low_bound=gsl_vector_alloc(param_count);
    high_bound=gsl_vector_alloc(param_count);

    if (vars->PDF_on)
    {
	do
	    vars->beta_init = vars->beta_mean+gsl_ran_gaussian (vars->ur, vars->beta_std);
	while (vars->beta_init < 0);
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

    /* bounds for beta */
    gsl_vector_set(params,0,vars->beta_init); 
    gsl_vector_set(low_bound,0,vars->beta_low); 
    gsl_vector_set(high_bound,0,vars->beta_high); 

    /* bounds for latent_inf */
    gsl_vector_set(params,1,vars->latent_inf_init); 
    gsl_vector_set(low_bound,1,vars->latent_inf_low); 
    gsl_vector_set(high_bound,1,vars->latent_inf_high); 

    /* bounds for p */
    gsl_vector_set(params,2,vars->log_p_init); 
    gsl_vector_set(low_bound,2,vars->log_p_low); 
    gsl_vector_set(high_bound,2,vars->log_p_high); 

    /* bounds for c */
    gsl_vector_set(params,3,vars->c_init); 
    gsl_vector_set(low_bound,3,vars->c_low); 
    gsl_vector_set(high_bound,3,vars->c_high); 

    /* bounds for theta */
    gsl_vector_set(params,4,vars->theta_init); 
    gsl_vector_set(low_bound,4,vars->theta_low); 
    gsl_vector_set(high_bound,4,vars->theta_high); 

    /* bounds for inf */
    gsl_vector_set(params,5,vars->inf_init); 
    gsl_vector_set(low_bound,5,vars->inf_low); 
    gsl_vector_set(high_bound,5,vars->inf_high); 

    /* bounds for r */
    gsl_vector_set(params,6,vars->r_init); 
    gsl_vector_set(low_bound,6,vars->r_low); 
    gsl_vector_set(high_bound,6,vars->r_high); 

    /* bounds for rinf */
    gsl_vector_set(params,7,vars->rinf_init); 
    gsl_vector_set(low_bound,7,vars->rinf_low); 
    gsl_vector_set(high_bound,7,vars->rinf_high); 

    /* bounds for delta */
    gsl_vector_set(params,8,vars->delta_init); 
    gsl_vector_set(low_bound,8,vars->delta_low); 
    gsl_vector_set(high_bound,8,vars->delta_high); 

    /* bounds for rho */
    gsl_vector_set(params,10,vars->rho_init); 
    gsl_vector_set(low_bound,10,vars->rho_low); 
    gsl_vector_set(high_bound,10,vars->rho_high); 

    /* bounds for eclipse */
    gsl_vector_set(params,11,vars->eclipse_init); 
    gsl_vector_set(low_bound,11,vars->eclipse_low); 
    gsl_vector_set(high_bound,11,vars->eclipse_high); 

    /* bounds for betaun */
    gsl_vector_set(params,12,vars->log_betaun_init); 
    gsl_vector_set(low_bound,12,vars->log_betaun_low); 
    gsl_vector_set(high_bound,12,vars->log_betaun_high); 

    if (vars->Model >= 5)
    {
	/* bounds for Cmax */
	gsl_vector_set(params,13,vars->Cmax_init); 
	gsl_vector_set(low_bound,13,vars->Cmax_low); 
	gsl_vector_set(high_bound,13,vars->Cmax_high); 

	/* bounds for IC50 */
	gsl_vector_set(params,14,vars->IC50_init); 
	gsl_vector_set(low_bound,14,vars->IC50_low); 
	gsl_vector_set(high_bound,14,vars->IC50_high); 

	/* bounds for m */
	gsl_vector_set(params,15,vars->m_init); 
	gsl_vector_set(low_bound,15,vars->m_low); 
	gsl_vector_set(high_bound,15,vars->m_high); 

	/* bounds for gamma */
	gsl_vector_set(params,16,vars->gamma_init); 
	gsl_vector_set(low_bound,16,vars->gamma_low); 
	gsl_vector_set(high_bound,16,vars->gamma_high); 

	/* bounds for absorb */
	gsl_vector_set(params,17,vars->absorb_init); 
	gsl_vector_set(low_bound,17,vars->absorb_low); 
	gsl_vector_set(high_bound,17,vars->absorb_high); 
    }

    if (vars->Fit_model > 0)
    {
	zoomin(ScoreFunction,vars->ur,params,&score,
	  low_bound,high_bound, vars->Param_mask,param_count,
	  vars->Max_steps,vars->Stop_walk,vars->Bvstop_walk,
	  (int)(vars->Printmax>0),vars->Tolerance,(void *)vars,vars->Threading,vars->Search_order);
    }
    else
    {
	if ((batchMode == 0 || vars->AutoSnapshot != 0) && vars->points == NULL)
	    vars->points = new plotPoints();
    
	vars->stopFlag = 0;
    
	vars->plot_bias = 0.;
	vars->hex_time_bias = 0;
	vars->sample_index = 0;

	score=ScoreFunction(&valid, params, (void *)vars,NULL);
    }

    gsl_vector_free(low_bound);
    gsl_vector_free(high_bound);

    return score;
}

/**************************************************************************
 * The following section contains all the callback function definitions.
 *  **************************************************************************/

png_byte color_type;
png_byte bit_depth;

png_structp png_ptr;
png_infop info_ptr;
int number_of_passes;
png_bytep * row_pointers;

void write_png_file( int width, int height,const char* file_name)
{
        /* create file */
        FILE *fp = fopen(file_name, "wb");
        if (!fp)
                abort_("[write_png_file] File %s could not be opened for writing", file_name);

        /* initialize stuff */
        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[write_png_file] png_create_write_struct failed");

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[write_png_file] png_create_info_struct failed");

	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);

	char *currentP=(char *)bgr;
	for (int j=0; j < height; j++)
	{
	    row_pointers[height-j-1]=(png_bytep)currentP;
	    currentP+=width*4;
	}

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during init_io");

        png_init_io(png_ptr, fp);


        /* write header */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing header");

        png_set_IHDR(png_ptr, info_ptr, width, height,
                     8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        png_write_info(png_ptr, info_ptr);


        /* write bytes */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing bytes");

        png_write_image(png_ptr, row_pointers);


        /* end write */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during end of write");

        png_write_end(png_ptr, NULL);

        /* cleanup heap allocation */
        free(row_pointers);

        fclose(fp);
}
#ifndef NO_GUI

bool ShowMessageBox(std::string title, std::string message)
{
    GtkWidget*dialog = gtk_message_dialog_new (GTK_WINDOW(main_window),
                                 GTK_DIALOG_DESTROY_WITH_PARENT,
                                 GTK_MESSAGE_ERROR,
                                 GTK_BUTTONS_CLOSE,
                                 message.c_str());

    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);

    return true;
}

void take_screenshot( GtkWidget *ok_button, char *filename )
{
	static GdkPixbuf *screenshot = NULL;

	screenshot = gdk_pixbuf_get_from_drawable( screenshot, GDK_DRAWABLE(GTK_WIDGET(image)->window), gdk_colormap_get_system(), 0, 0, 0, 0, scapture.enc_width, scapture.enc_height );
	gdk_pixbuf_save( GDK_PIXBUF(screenshot), filename, "png", NULL, NULL );
	//g_free( filename );
}


void get_movie_name( GtkWidget *ok_button, char *filename )
{
	gtk_entry_set_text( GTK_ENTRY(scapture.movie_name), filename );
	g_free(filename);
}


void file_name_changed( GtkWidget *widget, GtkWidget *check_button )
{
	if( strlen( gtk_entry_get_text(GTK_ENTRY(widget)) ) < 1 )
		gtk_widget_set_sensitive( GTK_WIDGET(check_button), FALSE );
	else
		gtk_widget_set_sensitive( GTK_WIDGET(check_button), TRUE );
}

void font_changed( GtkWidget *widget )
{
    PangoFontDescription *font_desc;
    PangoFont *font;
    string font_name;
    string font_size;


    font_name = gtk_entry_get_text(GTK_ENTRY(font_name_widget));
    font_size = gtk_entry_get_text(GTK_ENTRY(font_size_widget));

    font_name += " "+font_size;

    if( font_name.length() > 0 )
    {
	font_desc = pango_font_description_from_string (font_name.c_str());
        font_list_base = glGenLists (128);

	font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	if (font == NULL)
	{
	    string err_str ="*** Can't load font '"+font_name+"'\n";
	    ShowMessageBox("Warning",err_str);

	    g_print (err_str.c_str());
	    glDeleteLists(font_list_base,128);

	    font_list_base = default_font_list_base;
	    font_desc = default_font_desc;
	    font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	}
	else
	{
	    g_print ("*** Loaded font '%s'\n", font_name.c_str());
	    glDeleteLists(default_font_list_base,128);
	    default_font_list_base = font_list_base;
	    default_font_desc = font_desc;
	}
    }
}
void background_changed(GtkWidget *widget)
{
    double back_r = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_r) );
    double back_g = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_g) );
    double back_b = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_b) );

    glClearColor(back_r/255., back_g/255., back_b/255., 1.0);

    GdkColor init_color;

    init_color.red = 65535 * (back_r/255.);;
    init_color.green = 65535 * (back_g/255.);
    init_color.blue = 65535 * (back_b/255.);

    gtk_widget_modify_bg(rgb_box,GTK_STATE_NORMAL, &init_color);
}


void check_button_clicked( GtkWidget *widget, gpointer which_check_button )
{
    /* If it is the checkbox for saving animation */
    if( which_check_button == "movie" ) {
	if( GTK_TOGGLE_BUTTON( widget )->active ) 
	    theState->AutoSnapshot=1;
	else 
	    theState->AutoSnapshot=0;
    }
#ifdef OLD_WAY
	typedef enum {
		OGGM,
		Raw,
		WindowsMedia,
		QuickTime
	} encode_type_t;
	char movie_command[1000];
	FILE *pipe_file;
	int iter_codec = gtk_combo_box_get_active( GTK_COMBO_BOX(scapture.codec) );
	int nfps = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(scapture.framerate) );
	int nkbps = gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(scapture.bitrate) );
	int bytes_per_frame = scapture.enc_width*scapture.enc_height*3; /* 3 bytes = RGB */
	int length = 0;

	/* If it is the checkbox for saving animation */
	if( which_check_button == "movie" ) {
		if( GTK_TOGGLE_BUTTON( widget )->active ) {
			length += sprintf( movie_command, "gst-launch fdsrc fd=0 blocksize=%d ! video/x-raw-rgb, width=%d, height=%d, bpp=24, depth=24, red_mask=16711680, green_mask=65280, blue_mask=255, endianness=4321, framerate=%d/1 ! ffmpegcolorspace ! ", bytes_per_frame, scapture.enc_width, scapture.enc_height, nfps );
			switch( iter_codec ) {
				case Raw: /* Just dump it in an avi container: default extension = .avi */
					length += sprintf(movie_command+length, "avimux");
				break;
				case WindowsMedia: /* An msmpeg4 stream in an asf container: default extension = .wmv */
					length += sprintf(movie_command+length, "ffenc_msmpeg4 bitrate=%d ! ffmux_asf", nkbps);
				break;
				case QuickTime: /* An mpeg4 stream in a mov container: default extension = .mov */
					length += sprintf(movie_command+length, "ffenc_mpeg4 bitrate=%d ! ffmux_mov", nkbps);
				break;
				case OGGM:
				default: /* default extension = .ogm */
					length += sprintf(movie_command+length, "theoraenc bitrate=%d ! oggmux", nkbps);
				break;
			}
			sprintf(movie_command+length, " ! filesink location=\"%s\"", gtk_entry_get_text(GTK_ENTRY(scapture.movie_name)) );
			fprintf( stderr, "\n%s\n\n", movie_command );

			pipe_file = popen( movie_command, "w" );
			scapture.pipefd = fileno( pipe_file );
			if( !pipe_file ) {
				gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON(widget), FALSE );
				fprintf( stderr, "%s: cannot open movie pipe", PROGRAM_NAME );
				return;
			}
		}
		else {
			close( scapture.pipefd );
			scapture.pipefd = -1;
			//gtk_widget_set_sensitive( GTK_WIDGET(widget), FALSE );
		}
	}
#endif
}


void open_file_selection( GtkWidget *widget, gpointer data )
{
	GtkWidget *dialog;
	char *filename = NULL;
	void (*filename_handler)(GtkWidget *, char *) = (void (*)(GtkWidget*, char*))data;
	/* Function to call if a filename selected */

	dialog = gtk_file_chooser_dialog_new("Select File...", NULL, GTK_FILE_CHOOSER_ACTION_SAVE, GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT, NULL);
	if( gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT ) {
		filename = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(dialog) );
	}
	gtk_widget_destroy( dialog );
	if( filename != NULL )
		filename_handler( widget, filename );
}


static gboolean
initialize_image(GLfloat width, GLfloat height)
{
	guint count;

	float aspectRatio=(float)width/(float)height;
	GdkGLPixmap *glpixmap;

	//gtk_widget_set_size_request( GTK_WIDGET(image), width, height );
	global_gldrawable = gtk_widget_get_gl_drawable( image );
	global_glcontext = gtk_widget_get_gl_context( image );

	if( gdk_gl_drawable_gl_begin(global_gldrawable, global_glcontext) ) {

		glClear( GL_COLOR_BUFFER_BIT );
		glViewport(0, 0, width, height);
		glMatrixMode( GL_PROJECTION );
		 glDisable ( GL_LIGHTING ) ;
		glLoadIdentity();
		//if (width <= height)
		  //glOrtho(-60., 60., -60.0/aspectRatio, 60.0/aspectRatio, 1.0, -1.0);
		//else
		  //glOrtho(-60.0*aspectRatio, 60.0*aspectRatio,-60., 60.,  1.0, -1.0);
		  glOrtho(-60.0, 60.0,-60., 60.,  1.0, -1.0);
		//glOrtho( 0.0, width/arena->rescale_factor, height/arena->rescale_factor, 0.0, -1.0, 0.0 );
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glFlush();

		if (gdk_gl_drawable_is_double_buffered (global_gldrawable))
		  gdk_gl_drawable_swap_buffers (global_gldrawable);
		gdk_gl_drawable_gl_end( global_gldrawable );
	}
	gtk_widget_queue_draw( image );

	/* Encoders encode better with heights or widths that are multiples of 16 */
	/* so the area to be encoded is clipped to satisfy this. */
	scapture.enc_width = (int)width - (int)width%16;
	scapture.enc_height = (int)height - (int)height%16;
	return true;
}
/***
 *  *** The "realize" signal handler. All the OpenGL initialization
 *   *** should be performed here, such as default background colour,
 *    *** certain states etc.
 *     ***/
static void
realize (GtkWidget *widget,
         gpointer   data)
{
  GdkGLContext *glcontext = gtk_widget_get_gl_context (widget);
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);

  GLfloat width = widget->allocation.width;
  GLfloat height = widget->allocation.height;

  da_width = widget->allocation.width;
  da_height = widget->allocation.height;

  PangoFontDescription *font_desc;
  PangoFont *font;

  /*** OpenGL BEGIN ***/
  if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext))
    return;

  /*
   * Generate font display lists.
   */
  if (default_font_desc == NULL)
  {
        default_font_desc = pango_font_description_from_string (default_font_string);
        default_font_list_base = glGenLists (128);
  }

  font_list_base = default_font_list_base;

  font_desc = default_font_desc;

  font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
  if (font == NULL)
    {
      g_print ("*** Can't load font '%s'\n", default_font_string);
      exit (1);
    }

  /*** Fill in the details here. ***/
  //glClearColor(170./255., 180./255., 190./255., 1.0);

  double back_r = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_r) );
  double back_g = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_g) );
  double back_b = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_b) );

  glClearColor(back_r/255., back_g/255., back_b/255., 1.0);
  glClearDepth(1.0);
  /*glShadeModel(GL_FLAT);*/
  glDisable(GL_DEPTH_TEST);

  gdk_gl_drawable_gl_end (gldrawable);
  /*** OpenGL END ***/

  initialize_image(width, height);

  return;
}

static gboolean
configure_event (GtkWidget         *widget,
                 GdkEventConfigure *event)
{
	GLfloat width = widget->allocation.width;
	GLfloat height = widget->allocation.height;
	GLfloat w = event->width;
	GLfloat h = event->height;
	gboolean value = TRUE;

	if (width != da_width || height != da_height)
		realize(widget,NULL);
	else
		value =  initialize_image(width, height);

	return value;
}
void setCheckVi(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotVi = value;
        updateGL();
    }
}

void setCheckVe(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotVe = value;
        updateGL();
    }
}

void setCheckInf(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotInf = value;
        updateGL();
    }
}

void setCheckCd8(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotCd8 = value;
        updateGL();
    }
}

void setCheckACV(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotACV = value;
        updateGL();
    }
}

void setCheckLogs(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotLogs = value;
        updateGL();
    }
}

void setCheckColor(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotColor = value;
        updateGL();
    }
}

void setCheckRegions(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotRegions = value;
        updateGL();
    }
}

void setCheckWriteOn(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->writeOn = value;
        updateGL();
    }
}

void setCheckScrollAxes(GtkWidget *widget, gpointer dummy)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->scrollAxes = value;
        updateGL();
    }
}

static void
setStyle1 (GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotStyle1 = value;
    }
    gtk_widget_show(zoomin_button);
    gtk_widget_show(zoomout_button);
    gtk_widget_show(page_left_button);
    gtk_widget_show(page_right_button);
    updateGL();
}

void setOpt1(GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotOpt1 = value;
    }
    updateGL();
}
void setOpt2(GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotOpt2 = value;
    }
    updateGL();
}
void setOpt3(GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotOpt3 = value;
    }
    updateGL();
}
void setOpt4(GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotOpt4 = value;
    }
    updateGL();
}
void setOpt5(GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotOpt5 = value;
    }
    updateGL();
}
void setOpt6(GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotOpt6 = value;
    }
    updateGL();
}

void setOpt7(GtkWidget* widget, gpointer data)
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotOpt7 = value;
    }
    updateGL();
}

static double plot_span_saved=10.;
static double plot_bias_saved=0.;
static void
run_cb (GtkWidget* widget, gpointer data)
{
    gtk_widget_set_sensitive( GTK_WIDGET(run_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(pause_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(stop_button), TRUE );

    if (theState->pauseFlag == 1) {
	theState->pauseFlag = 0;
	//theState->plot_span = plot_span_saved;
	theState->plot_bias = plot_bias_saved;
	//theState->hex_time_bias = 0.0;
	//gtk_widget_set_sensitive( GTK_WIDGET(zoomin_button), FALSE );
	//gtk_widget_set_sensitive( GTK_WIDGET(zoomout_button), FALSE );
	gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), FALSE );
	gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), FALSE );
	gtk_widget_hide( GTK_WIDGET(backup_slider) );
	updateGL();
	return;
    }

    // try to prevent double clicking or premature sim launches
    if (simLock != 0) {
	//theState->app->beep();
	return;
    }


    simLock = 1;

    if (theState->points != NULL)
	delete theState->points;

    theState->points = new plotPoints();

    theState->stopFlag = 0;

    theState->plot_bias = 0.;
    theState->hex_time_bias = 0;
    theState->sample_index = 0;

    /* reset scroll bar if necessary */
    ((GtkAdjustment *)adjuster)->value = 20;
    gtk_signal_emit_by_name(GTK_OBJECT(adjuster), "value_changed");

    gtk_widget_hide( GTK_WIDGET(backup_slider) );

    char value_str[100];
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(beta_entry)));
    sscanf(value_str,"%lf", &theState->beta_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(c_entry)));
    sscanf(value_str,"%lf", &theState->c_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(delta_entry)));
    sscanf(value_str,"%lf", &theState->delta_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(eclipse_entry)));
    sscanf(value_str,"%lf", &theState->eclipse_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(inf_entry)));
    sscanf(value_str,"%lf", &theState->inf_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(p_entry)));
    sscanf(value_str,"%lf", &theState->log_p_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(r_entry)));
    sscanf(value_str,"%lf", &theState->r_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(rho_entry)));
    sscanf(value_str,"%lf", &theState->rho_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(rinf_entry)));
    sscanf(value_str,"%lf", &theState->rinf_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(theta_entry)));
    sscanf(value_str,"%lf", &theState->theta_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(latent_inf_entry)));
    sscanf(value_str,"%lf", &theState->latent_inf_init);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(T0_entry)));
    sscanf(value_str,"%lf", &theState->To);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(plot_span_entry)));
    sscanf(value_str,"%lf", &theState->plot_span);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(use_rho_entry)));
    sscanf(value_str,"%d", &theState->use_rho);

    runSimulation(theState);  

    simLock = 0;

    updateGL();
}
static void
pause_cb (GtkWidget* widget, gpointer data)
{
    // disallow double clicking or premature sim stopping
    if (theState->pauseFlag ==1) {
	//theState->app->beep();
	return;
    }

    theState->pauseFlag = 1;
    //plot_span_saved = theState->plot_span;
    plot_bias_saved = theState->plot_bias;
    //gtk_widget_set_sensitive( GTK_WIDGET(zoomin_button), TRUE );
    //gtk_widget_set_sensitive( GTK_WIDGET(zoomout_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(pause_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(run_button), TRUE );
    gtk_widget_show( GTK_WIDGET(backup_slider) );
    updateGL();
}

static void
stop_cb (GtkWidget* widget, gpointer data)
{
    // disallow double clicking or premature sim stopping
    if (simLock == 0) {
	//theState->app->beep();
	return;
    }

    theState->stopFlag = 1;
    theState->pauseFlag = 0;

    //gtk_widget_set_sensitive( GTK_WIDGET(zoomin_button), TRUE );
    //gtk_widget_set_sensitive( GTK_WIDGET(zoomout_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(stop_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(pause_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(run_button), TRUE );
    gtk_widget_hide( GTK_WIDGET(backup_slider) );
    updateGL();
}

static void
zoomin_cb (GtkWidget* widget, gpointer data)
{
    char value_str[100];
    if (theState->pauseFlag ==0) {
	//theState->app->beep();
	//return;
    }
    theState->plot_span = theState->plot_span/2.0;
    sprintf(value_str,"%lf", theState->plot_span);
    gtk_entry_set_text (GTK_ENTRY(plot_span_entry), value_str);

    //emit param8Changed(QString::number(theState->plot_span));
    updateGL();
}
static void
zoomout_cb (GtkWidget* widget, gpointer data)
{
    char value_str[100];
    if (theState->pauseFlag ==0) {
	//theState->app->beep();
	//return;
    }
    theState->plot_span = theState->plot_span*2.0;
    sprintf(value_str,"%lf", theState->plot_span);
    gtk_entry_set_text (GTK_ENTRY(plot_span_entry), value_str);

    //emit param8Changed(QString::number(theState->plot_span));
    updateGL();
}

static void
page_left_cb (GtkWidget* widget, gpointer data)
{
    // disallow scrolling during run mode
    if (theState->pauseFlag !=1 && theState->stopFlag != 1) {
	//theState->app->beep();
	return;
    }
    if (theState->plotStyle1 != 0)
    {
    if (theState->time+(theState->plot_bias- theState->plot_span) >= 0.0)
	theState->plot_bias -= theState->plot_span;
    } 
    updateGL();
}
static void
page_right_cb (GtkWidget* widget, gpointer data)
{
    // disallow scrolling during run mode
    if (theState->pauseFlag !=1 && theState->stopFlag != 1) {
	//theState->app->beep();
	return;
    }
    if (theState->plotStyle1 != 0)
    {
    if (theState->plot_bias < 0.0)
	theState->plot_bias += theState->plot_span;

    }
    updateGL();
}

void scrollTime(GtkWidget* widget, gpointer data)
{
    double value;

    // disallow scrolling during run mode
    if (theState->pauseFlag !=1 && theState->stopFlag != 1) {
	//theState->app->beep();
	return;
    }
    value = gtk_adjustment_get_value(GTK_ADJUSTMENT(adjuster));

    theState->hex_time_bias = value - 20.0;

    if (theState->hex_time_bias < 0 && theState->time+theState->hex_time_bias < 0)
	theState->hex_time_bias = -theState->time;

    updateGL();
}

void snap_movie_frame(void)
{
	static GdkPixbuf *screenshot = NULL;

	if( !(scapture.pipefd < 0) ) {
		screenshot = gdk_pixbuf_get_from_drawable( screenshot, GDK_DRAWABLE(GTK_WIDGET(image)->window), gdk_colormap_get_system(), 0, 0, 0, 0, scapture.enc_width, scapture.enc_height );
		write( scapture.pipefd, gdk_pixbuf_get_pixels(GDK_PIXBUF(screenshot)), gdk_pixbuf_get_rowstride(GDK_PIXBUF(screenshot))*gdk_pixbuf_get_height(GDK_PIXBUF(screenshot)) );
	}
}

void updateGL(void)
{
    //advance_handler();
    gdk_window_invalidate_rect (image->window, &image->allocation, FALSE);
}
//! [8-]
void resizeGL(int width, int height)
{
//    int side = qMin(width, height);
//    glViewport((width - side) / 2, (height - side) / 2, side, side);

    float aspectRatio=(float)width/(float)height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
		 glDisable ( GL_LIGHTING ) ;
    glLoadIdentity();
if (width <= height)
    glOrtho(-60., 60., -60.0/aspectRatio, 60.0/aspectRatio, 1.0, -1.0);
else
    glOrtho(-60.0*aspectRatio, 60.0*aspectRatio,-60., 60.,  1.0, -1.0);
    glMatrixMode(GL_MODELVIEW);
    glDisable(GL_DEPTH_TEST);
    glLoadIdentity();
}
/***
 *  *** The "expose_event" signal handler. All the OpenGL re-drawing should
 *   *** be done here. This is repeatedly called as the painting routine
 *    *** every time the 'expose'/'draw' event is signalled.
 *     ***/
static gboolean
expose_event (GtkWidget      *widget,
              GdkEventExpose *event,
              gpointer        data)
{
  GdkGLContext *glcontext = gtk_widget_get_gl_context (widget);
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);

GLfloat width = widget->allocation.width;
GLfloat height = widget->allocation.height;
static GdkPixbuf *screenshot = NULL;
  char val_str[1000];

  /*** OpenGL BEGIN ***/
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext))
	return FALSE;

  if (width != da_width || height != da_height)
	realize(widget,NULL);

    draw_routine (width,height);
  /* Swap buffers */
    gdk_gl_drawable_swap_buffers (gldrawable);

  gdk_gl_drawable_gl_end (gldrawable);
  /*** OpenGL END ***/

  return true;
}
/***
 *  *** The "button_press_event" signal handler. Any processing required when
 *   *** mouse buttons (only left and middle buttons) are pressed on the OpenGL-
 *    *** capable drawing area should be done here.
 *     ***/
static gboolean
button_press_event (GtkWidget      *widget,
                    GdkEventButton *event,
                    gpointer        data)
{
  int x,y;
  if (event->button == 1)
    {
      x = event->x;
      y = widget->allocation.height - event->y;

      gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);
      return TRUE;
    }

  return FALSE;
}

static void
quit_cb (GtkWidget* widget, gpointer data)
{
  gtk_main_quit();
}

void destroy_widget( GtkWidget *widget, gpointer data)
{
	GtkWidget *widget_to_kill = (GtkWidget *) data;

	gtk_widget_destroy( widget_to_kill );
}

void destroy( void )
{
	if (theState != NULL)
	{
	    theState->pauseFlag = 0;
	    theState->stopFlag = 1;
	}

	gtk_main_quit();
}

static void attr_list_insert( PangoAttrList *attrlist, PangoAttribute *attr )
{
	attr->start_index = 0;
	attr->end_index = G_MAXINT;
	pango_attr_list_insert( attrlist, attr );
}

void select_font(GtkWidget *widget, gpointer label)
{

  GtkResponseType result;

  GtkWidget *dialog = gtk_font_selection_dialog_new("Select Font");
  do
  {
      result = (GtkResponseType)gtk_dialog_run(GTK_DIALOG(dialog));

      if (result == GTK_RESPONSE_OK || result == GTK_RESPONSE_APPLY)
      {

	PangoFontDescription *font_desc;
	PangoFont *font;
	gchar *fontname = gtk_font_selection_dialog_get_font_name(
				GTK_FONT_SELECTION_DIALOG(dialog));

	font_desc = pango_font_description_from_string(fontname);
	font_list_base = glGenLists (128);

	font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	if (font == NULL)
	{
	    string err_str ="*** Can't use font '";
	    err_str +=fontname;
	    err_str +="' for OpenGL drawing\n";

	    ShowMessageBox("Warning",err_str);

	    g_print (err_str.c_str());
	    glDeleteLists(font_list_base,128);

	    font_list_base = default_font_list_base;
	    font_desc = default_font_desc;
	    font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	}
	else
	{
	    g_print ("*** Loaded font '%s'\n", fontname);
	    glDeleteLists(default_font_list_base,128);
	    default_font_list_base = font_list_base;
	    default_font_desc = font_desc;
	    result = GTK_RESPONSE_CANCEL;
	}

	gtk_widget_modify_font(GTK_WIDGET(label), font_desc);

	g_free(fontname);
      }
  }
  while(result == GTK_RESPONSE_OK || result == GTK_RESPONSE_APPLY);

  gtk_widget_destroy(dialog);
}



/**************************************************************************
 *  * The following section contains the GUI building function definitions.
 *   **************************************************************************/

void about_window( void )
{
	GtkWidget *about_window;
	GtkWidget *button;
	GtkWidget *frame;
	GtkWidget *hseparator;
	GtkWidget *label;
	GtkWidget *vbox, *hbox;
	GdkPixmap *apixmap;
	GdkBitmap *mask;
	GtkWidget *aimage;
	PangoAttrList *attrlist;

	about_window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
	gtk_window_set_title( GTK_WINDOW(about_window), "About hhv8_sim");
	g_signal_connect( GTK_OBJECT(about_window), "destroy", G_CALLBACK(destroy_widget), about_window );

	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 15);
	gtk_container_add( GTK_CONTAINER( about_window ), vbox );

	button = gtk_button_new_with_label( "OK" );
	g_signal_connect( GTK_OBJECT( button ), "clicked", G_CALLBACK( destroy_widget ), about_window );
	gtk_box_pack_end( GTK_BOX( vbox ), button, FALSE, FALSE, 0 );

	hseparator = gtk_hseparator_new();
	gtk_box_pack_end( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( vbox ), hbox, FALSE, FALSE, 0 );

	apixmap = gdk_pixmap_create_from_xpm_d( main_window->window, &mask, NULL, hhv8_xpm );
	aimage = gtk_image_new_from_pixmap( apixmap, mask );
	gtk_container_add( GTK_CONTAINER( hbox ), aimage );

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_box_pack_start( GTK_BOX( hbox ), vbox, FALSE, FALSE, 0 );

	attrlist = pango_attr_list_new();
	attr_list_insert( attrlist, pango_attr_size_new( 20*PANGO_SCALE ) );
	attr_list_insert( attrlist, pango_attr_weight_new( PANGO_WEIGHT_BOLD ) );
	attr_list_insert( attrlist, pango_attr_foreground_new( 8738, 39835, 5911 ) );
	label = gtk_label_new( DEFAULT_TITLE " - Version " PACKAGE_VERSION );
	gtk_label_set_attributes( GTK_LABEL(label), attrlist );
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_CENTER );
	gtk_label_set_line_wrap( GTK_LABEL(label), FALSE );
	gtk_box_pack_start( GTK_BOX( vbox ), label, FALSE, FALSE, 0 );
	pango_attr_list_unref( attrlist );

	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type( GTK_FRAME( frame ), GTK_SHADOW_IN );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0 );	

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 10);
	gtk_container_add( GTK_CONTAINER( frame ), vbox );

	label = gtk_label_new( "A Visualization Program for recurrent HHV8 infection\nDistributed under the GNU GPL v3.0" );
	gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	gtk_box_pack_start( GTK_BOX( vbox ), label, FALSE, FALSE, 0 );	

	label = gtk_label_new( "Authors:\n  David Swan <"PACKAGE_BUGREPORT">\n  Dr. Joshua Schiffer" );
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	gtk_box_pack_start( GTK_BOX( vbox ), label, FALSE, FALSE, 0 );

	gtk_widget_show_all(about_window);
}


void make_gui(void)
{
#ifdef __sun
	static GtkActionEntry menu_entries[] = {
		{ "FileMenu", NULL, "_File" },
		{ "HelpMenu", NULL, "_Help" },
		{ "Quit", GTK_STOCK_QUIT, "_Quit", "<control>Q", "Exit the program", destroy },
		{ "HelpAbout", NULL, "_About", NULL, "About hhv8_sim...", about_window },
	};
#else
	static GtkActionEntry menu_entries[] = {
		{ "FileMenu", NULL, "_File" },
		{ "HelpMenu", NULL, "_Help" },
		{ "Quit", GTK_STOCK_QUIT, "_Quit", "<control>Q", "Exit the program", destroy },
		{ "HelpAbout", GTK_STOCK_ABOUT, "_About", NULL, "About hhv8_sim...", about_window },
	};
#endif
	static const char *ui_description =
	"<ui>"
	"  <menubar name='MainMenu'>"
	"    <menu action='FileMenu'>"
	"      <menuitem action='Quit'/>"
	"    </menu>"
	"    <menu action='HelpMenu'>"
	"      <menuitem action='HelpAbout'/>"
	"    </menu>"
	"  </menubar>"
	"</ui>";
	GtkActionGroup *action_group;
	GtkUIManager *ui_manager;
	GtkAccelGroup *accel_group;
	GError *error;
	GtkWidget *menubar;
	GtkWidget *moviebutton = NULL;

	struct rt_t {
		GtkWidget *time;
		GtkWidget *units;
	} real_time;

	GtkObject *spinbutton_adj;
	GdkPixmap *apixmap;
	GdkBitmap *mask;

	GtkWidget *button;
	GtkWidget *frame;
	GtkWidget *combo;
	GtkWidget *label;
	GtkWidget *notebook;
	GtkWidget *scrolledwindow;
	GtkWidget *vbox,*svbox,*hbox;
	GtkWidget *optBox;

	GtkWidget *plotLines;
	GtkWidget *drawVe;
	GtkWidget *drawVi;
	GtkWidget *drawInf;
	GtkWidget *drawCD8s;
	GtkWidget *drawR0s;
	GtkWidget *drawNums;

	GtkWidget *plotColors;
	GtkWidget *plotLogs;
	GtkWidget *plotRegions;
	GtkWidget *writeEnable;
	GtkWidget *scrollAxes;

	GtkWidget *plot_ve_button;
	GtkWidget *plot_vi_button;
	GtkWidget *plot_inf_button;
	GtkWidget *plot_cd8_button;
	GtkWidget *plot_ACV_button;


	/* Creates the main Window ; sets its title ; attaches the closing behavior */
	main_window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
	gtk_window_set_default_size( GTK_WINDOW(main_window), DEFAULT_WIDTH,DEFAULT_HEIGHT);
	gtk_window_set_title( GTK_WINDOW(main_window), DEFAULT_TITLE);
	g_signal_connect( GTK_OBJECT(main_window), "destroy", G_CALLBACK(destroy), NULL );
	gtk_quit_add_destroy(1, GTK_OBJECT(main_window));
        /* Connect signal handlers to the window */
        //g_signal_connect (G_OBJECT (main_window), "delete_event", G_CALLBACK (quit_cb), NULL);
	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_container_add( GTK_CONTAINER( main_window ), vbox );

	/* Create Menu */
	action_group = gtk_action_group_new ("MenuActions");
	gtk_action_group_add_actions( action_group, menu_entries, G_N_ELEMENTS(menu_entries), main_window);
	ui_manager = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group( ui_manager, action_group, 0 );
	accel_group = gtk_ui_manager_get_accel_group( ui_manager );
	gtk_window_add_accel_group( GTK_WINDOW(main_window), accel_group );
	error = NULL;
	if( !gtk_ui_manager_add_ui_from_string( ui_manager, ui_description, -1, &error ) ) {
		g_message( "Building menus failed: %s", error->message );
		g_error_free( error );
		exit( EXIT_FAILURE );
	}
	menubar = gtk_ui_manager_get_widget( ui_manager, "/MainMenu" );
	gtk_box_pack_start( GTK_BOX( vbox ), menubar, FALSE, FALSE, 0 );

	/* Create the Notebook */
	notebook = gtk_notebook_new();
	gtk_notebook_set_tab_pos( GTK_NOTEBOOK( notebook ), GTK_POS_TOP );
	gtk_box_pack_start( GTK_BOX( vbox ), notebook, TRUE, TRUE, 0 );


	/* Create Simulation TAB */
	vbox = gtk_vbox_new( FALSE, 10 );
	gtk_container_set_border_width( GTK_CONTAINER( vbox ), 10 );
	gtk_notebook_append_page( GTK_NOTEBOOK( notebook ), vbox, gtk_label_new( "Simulation" ) );
	/* Create frame to hold the graph */
	hbox = gtk_hbox_new( FALSE, 0 );
	gtk_box_pack_start( GTK_BOX( vbox ), hbox, TRUE, TRUE, 0 );
	//frame = gtk_frame_new(NULL);
	//gtk_frame_set_shadow_type( GTK_FRAME(frame), GTK_SHADOW_IN );
	//gtk_box_pack_start( GTK_BOX( hbox ), frame, TRUE, TRUE, 0 );
	gtk_widget_show (hbox);

	/* Create glconfig attributes for gtkglext */
	glconfig = gdk_gl_config_new_by_mode( (GdkGLConfigMode)(GDK_GL_MODE_RGB | GDK_GL_MODE_DOUBLE ));
	if( glconfig == NULL ) {
		g_print ("*** Problem with gtkglext.\n");
		exit(1);
	}
	image = gtk_drawing_area_new();
	gtk_widget_set_size_request( GTK_WIDGET(image), DEFAULT_WIDTH, DEFAULT_HEIGHT);
	gtk_widget_set_gl_capability( image, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE );

	//gtk_widget_add_events(GTK_WIDGET(image), GDK_CONFIGURE);

	/* Connect signal handlers to the drawing area */
	g_signal_connect_after (GTK_OBJECT (image), "realize",
				G_CALLBACK (realize), NULL);
	g_signal_connect_after (GTK_OBJECT (image), "configure_event",
				G_CALLBACK (configure_event), NULL);
	g_signal_connect(GTK_OBJECT(image), "expose_event", G_CALLBACK(expose_event), NULL);
	//gtk_container_add( GTK_CONTAINER(frame), image );

	gtk_box_pack_start( GTK_BOX( hbox ), image, TRUE, TRUE, 0 );

	/* Create vbox to hold plot settings */
	optBox = gtk_vbox_new( FALSE, 0 );
	gtk_box_pack_end( GTK_BOX(hbox), optBox, FALSE, FALSE, 0 );
	gtk_widget_show (optBox);

	/* Create a check box for lines plots */
	plotLines = gtk_check_button_new_with_label("Plot vs. time");
	g_signal_connect (G_OBJECT (plotLines), "toggled",
			  G_CALLBACK (setStyle1), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plotLines, TRUE, FALSE, 0 );
	gtk_widget_show (plotLines);


	/* Create check box for cell plot options */
	drawVe = gtk_check_button_new_with_label("Show Ve by Region");
	g_signal_connect (G_OBJECT (drawVe), "toggled",
			  G_CALLBACK (setOpt1), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), drawVe, TRUE, FALSE, 0 );
	gtk_widget_show (drawVe);

	drawVi = gtk_check_button_new_with_label("Show Vi by Region");
	g_signal_connect (G_OBJECT (drawVi), "toggled",
			  G_CALLBACK (setOpt6), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), drawVi, TRUE, FALSE, 0 );
	gtk_widget_show (drawVi);

	drawInf = gtk_check_button_new_with_label("Show Inf Cells by Region");
	g_signal_connect (G_OBJECT (drawInf), "toggled",
			  G_CALLBACK (setOpt2), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), drawInf, TRUE, FALSE, 0 );

	drawCD8s = gtk_check_button_new_with_label("Show CD8 Counts by Region");
	g_signal_connect (G_OBJECT (drawCD8s), "toggled",
			  G_CALLBACK (setOpt3), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), drawCD8s, TRUE, FALSE, 0 );

	drawR0s = gtk_check_button_new_with_label("Show Repro Num by Region");
	g_signal_connect (G_OBJECT (drawR0s), "toggled",
			  G_CALLBACK (setOpt4), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), drawR0s, TRUE, FALSE, 0 );

	drawNums = gtk_check_button_new_with_label("Show Region Numbering");
	g_signal_connect (G_OBJECT (drawNums), "toggled",
			  G_CALLBACK (setOpt7), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), drawNums, TRUE, FALSE, 0 );

	plotColors = gtk_check_button_new_with_label( "Color Cells" );
	g_signal_connect( GTK_OBJECT( plotColors ), "toggled", G_CALLBACK(setCheckColor), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plotColors, TRUE, FALSE, 0 );

	plotLogs = gtk_check_button_new_with_label( "Use Log values" );
	g_signal_connect( GTK_OBJECT( plotLogs ), "toggled", G_CALLBACK(setCheckLogs), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plotLogs, TRUE, FALSE, 0 );

	plotRegions = gtk_check_button_new_with_label( "Plot ve/vi by region" );
	g_signal_connect( GTK_OBJECT( plotRegions ), "toggled", G_CALLBACK(setCheckRegions), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plotRegions, TRUE, FALSE, 0 );

	writeEnable = gtk_check_button_new_with_label( "Enable data writing" );
	g_signal_connect( GTK_OBJECT( writeEnable ), "toggled", G_CALLBACK(setCheckWriteOn), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), writeEnable, TRUE, FALSE, 0 );

	scrollAxes = gtk_check_button_new_with_label( "Scroll time axes" );
	g_signal_connect( GTK_OBJECT( scrollAxes ), "toggled", G_CALLBACK(setCheckScrollAxes), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), scrollAxes, TRUE, FALSE, 0 );

	/* Create the bottom box for pause, play, etc. buttons */
	hbox = gtk_hbox_new( FALSE, 10 );
	/*
       *    * Run simulation button.
       *       */

	run_button = gtk_button_new_from_stock( GTK_STOCK_GO_FORWARD );

	g_signal_connect (G_OBJECT (run_button), "clicked",
			  G_CALLBACK (run_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), run_button, TRUE, TRUE, 5);
	
	gtk_widget_show (run_button);
	/*
       *    * Stop simulation button.
       *       */

	stop_button = gtk_button_new_from_stock( GTK_STOCK_STOP );

	g_signal_connect (G_OBJECT (stop_button), "clicked",
			  G_CALLBACK (stop_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), stop_button, TRUE, TRUE, 5);
	
	gtk_widget_show (stop_button);
	gtk_widget_set_sensitive( GTK_WIDGET(stop_button), FALSE );
	/*
       *    * Pause simulation button.
       *       */

	pause_button = gtk_button_new_with_label( "Pause" );

	g_signal_connect (G_OBJECT (pause_button), "clicked",
			  G_CALLBACK (pause_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), pause_button, TRUE, TRUE, 5);
	
	gtk_widget_show (pause_button);
	gtk_widget_set_sensitive( GTK_WIDGET(pause_button), FALSE );

	/* Create cell time scroll bar (hidden for now) */
        gdouble lower = 0.0;
        gdouble upper = 20.01;
        gdouble step_increment = 0.01;
        gdouble page_increment = 0.1;
        gdouble page_size = 0.01;

	adjuster = gtk_adjustment_new (20.0, lower, upper, step_increment,
                                                         page_increment, page_size);


	backup_slider = gtk_hscrollbar_new(GTK_ADJUSTMENT(adjuster));

	gtk_box_pack_start (GTK_BOX (hbox), backup_slider, TRUE, TRUE, 5);

	g_signal_connect (G_OBJECT (adjuster), "value-changed",
			  G_CALLBACK (scrollTime), NULL);
	/*
       *    * Zoom in simulation button.
       *       */

	zoomin_button = gtk_button_new_with_label ("Zoom in");

	g_signal_connect (G_OBJECT (zoomin_button), "clicked",
			  G_CALLBACK (zoomin_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), zoomin_button, TRUE, TRUE, 5);
	
	gtk_widget_show (zoomin_button);

	/*
       *    * Zoom out simulation button.
       *       */

	zoomout_button = gtk_button_new_with_label ("Zoom out");

	g_signal_connect (G_OBJECT (zoomout_button), "clicked",
			  G_CALLBACK (zoomout_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), zoomout_button, TRUE, TRUE, 5);
	
	gtk_widget_show (zoomout_button);

	gtk_widget_set_sensitive( GTK_WIDGET(zoomin_button), TRUE );
	gtk_widget_set_sensitive( GTK_WIDGET(zoomout_button), TRUE );
	/*
       *    * Page left simulation button.
       *       */

	page_left_button = gtk_button_new_with_label ("<<");

	g_signal_connect (G_OBJECT (page_left_button), "clicked",
			  G_CALLBACK (page_left_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), page_left_button, TRUE, TRUE, 5);
	
	gtk_widget_show (page_left_button);

	/*
       *    * Page right simulation button.
       *       */

	page_right_button = gtk_button_new_with_label (">>");

	g_signal_connect (G_OBJECT (page_right_button), "clicked",
			  G_CALLBACK (page_right_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), page_right_button, TRUE, TRUE, 5);
	
	gtk_widget_show (page_right_button);
	gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), FALSE );
	gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), FALSE );

	/* Button to take screenshot */
	button = gtk_button_new_with_label( "Screenshot" );
	g_signal_connect( GTK_OBJECT( button ), "clicked",G_CALLBACK( open_file_selection ), (void *)&take_screenshot );
	gtk_box_pack_start( GTK_BOX( hbox ), button, TRUE, TRUE, 5 );
	gtk_widget_show (button);
	/* Record to movie */
	moviebutton = gtk_check_button_new_with_label( "Save animation." );
	g_signal_connect( GTK_OBJECT( moviebutton ), "clicked", G_CALLBACK(check_button_clicked), (void *)("movie"));
	//gtk_widget_set_sensitive( GTK_WIDGET(moviebutton), FALSE );
	gtk_box_pack_start( GTK_BOX( hbox ), moviebutton, TRUE, TRUE, 5 );
	gtk_widget_show (moviebutton);
	gtk_widget_show (hbox);
	gtk_box_pack_end( GTK_BOX( vbox ), hbox, FALSE, FALSE, 0 );

	/* Create Properties TAB */
	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_notebook_append_page( GTK_NOTEBOOK( notebook ), vbox, gtk_label_new( "Properties" ) );
	gtk_container_set_border_width( GTK_CONTAINER( vbox ), 10 );

	/* Parameters frame */
	frame = gtk_frame_new( "Simulation Parameters" );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0);

	GtkWidget *table = gtk_table_new(5, 6, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
	gtk_widget_show (table);

  	GtkWidget *label1 = gtk_label_new("Beta:");
  	GtkWidget *label2 = gtk_label_new("Beta_e:");
  	GtkWidget *label3 = gtk_label_new("C:");
  	GtkWidget *label4 = gtk_label_new("Delta:");
  	GtkWidget *label5 = gtk_label_new("Inf:");
  	GtkWidget *label6 = gtk_label_new("Log P:");
  	GtkWidget *label7 = gtk_label_new("R:");
  	GtkWidget *label8 = gtk_label_new("Rho:");
  	GtkWidget *label9 = gtk_label_new("RInf:");
  	GtkWidget *label10 = gtk_label_new("Theta:");
  	GtkWidget *label11 = gtk_label_new("Vburstrate:");
  	GtkWidget *label12 = gtk_label_new("Eclipse:");
  	GtkWidget *label13 = gtk_label_new("T0:");
  	GtkWidget *label14 = gtk_label_new("Plot span:");
  	GtkWidget *label15 = gtk_label_new("Use rho:");

	beta_entry = gtk_entry_new();
	c_entry = gtk_entry_new();
	delta_entry = gtk_entry_new();
	eclipse_entry = gtk_entry_new();
	inf_entry = gtk_entry_new();
	p_entry = gtk_entry_new();
	r_entry = gtk_entry_new();
	rho_entry = gtk_entry_new();
	rinf_entry = gtk_entry_new();
	theta_entry = gtk_entry_new();
	latent_inf_entry = gtk_entry_new();
	eclipse_entry = gtk_entry_new();
	T0_entry = gtk_entry_new();
	plot_span_entry = gtk_entry_new();
	use_rho_entry = gtk_entry_new();

	gtk_table_attach_defaults(GTK_TABLE(table), label1, 0, 1, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), label2, 0, 1, 1, 2);
	gtk_table_attach_defaults(GTK_TABLE(table), label3, 0, 1, 2, 3);
	gtk_table_attach_defaults(GTK_TABLE(table), label4, 0, 1, 3, 4);
	gtk_table_attach_defaults(GTK_TABLE(table), label5, 0, 1, 4, 5);
	gtk_table_attach_defaults(GTK_TABLE(table), label6, 2, 3, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), label7, 2, 3, 1, 2);
	gtk_table_attach_defaults(GTK_TABLE(table), label8, 2, 3, 2, 3);
	gtk_table_attach_defaults(GTK_TABLE(table), label9, 2, 3, 3, 4);
	gtk_table_attach_defaults(GTK_TABLE(table), label10, 2, 3, 4, 5);
	gtk_table_attach_defaults(GTK_TABLE(table), label11, 4, 5, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), label12, 4, 5, 1, 2);
	gtk_table_attach_defaults(GTK_TABLE(table), label13, 4, 5, 2, 3);
	gtk_table_attach_defaults(GTK_TABLE(table), label14, 4, 5, 3, 4);
	gtk_table_attach_defaults(GTK_TABLE(table), label15, 4, 5, 4, 5);


	gtk_table_attach_defaults(GTK_TABLE(table), beta_entry, 1, 2, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), c_entry, 1, 2, 2, 3);
	gtk_table_attach_defaults(GTK_TABLE(table), delta_entry, 1, 2, 3, 4);
	gtk_table_attach_defaults(GTK_TABLE(table), inf_entry, 1, 2, 4, 5);
	gtk_table_attach_defaults(GTK_TABLE(table), p_entry, 3, 4, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), r_entry, 3, 4, 1, 2);
	gtk_table_attach_defaults(GTK_TABLE(table), rho_entry, 3, 4, 2, 3);
	gtk_table_attach_defaults(GTK_TABLE(table), rinf_entry, 3, 4, 3, 4);
	gtk_table_attach_defaults(GTK_TABLE(table), theta_entry, 3, 4, 4, 5);
	gtk_table_attach_defaults(GTK_TABLE(table), latent_inf_entry, 5, 6, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), eclipse_entry, 5, 6, 1, 2);
	gtk_table_attach_defaults(GTK_TABLE(table), T0_entry, 5, 6, 2, 3);
	gtk_table_attach_defaults(GTK_TABLE(table), plot_span_entry, 5, 6, 3, 4);
	gtk_table_attach_defaults(GTK_TABLE(table), use_rho_entry, 5, 6, 4, 5);

	char value_str[100];
	sprintf(value_str,"%lf", theState->beta_init);
	gtk_entry_set_text (GTK_ENTRY(beta_entry), value_str);
	sprintf(value_str,"%lf", theState->c_init);
	gtk_entry_set_text (GTK_ENTRY(c_entry), value_str);
	sprintf(value_str,"%lf", theState->delta_init);
	gtk_entry_set_text (GTK_ENTRY(delta_entry), value_str);
	sprintf(value_str,"%lf", theState->eclipse_init);
	gtk_entry_set_text (GTK_ENTRY(eclipse_entry), value_str);
	sprintf(value_str,"%lf", theState->inf_init);
	gtk_entry_set_text (GTK_ENTRY(inf_entry), value_str);
	sprintf(value_str,"%lf", theState->log_p_init);
	gtk_entry_set_text (GTK_ENTRY(p_entry), value_str);
	sprintf(value_str,"%lf", theState->r_init);
	gtk_entry_set_text (GTK_ENTRY(r_entry), value_str);
	sprintf(value_str,"%lf", theState->rho_init);
	gtk_entry_set_text (GTK_ENTRY(rho_entry), value_str);
	sprintf(value_str,"%lf", theState->rinf_init);
	gtk_entry_set_text (GTK_ENTRY(rinf_entry), value_str);
	sprintf(value_str,"%lf", theState->theta_init);
	gtk_entry_set_text (GTK_ENTRY(theta_entry), value_str);
	sprintf(value_str,"%lf", theState->latent_inf_init);
	gtk_entry_set_text (GTK_ENTRY(latent_inf_entry), value_str);
	sprintf(value_str,"%lf", theState->To);
	gtk_entry_set_text (GTK_ENTRY(T0_entry), value_str);
	sprintf(value_str,"%lf", theState->plot_span);
	gtk_entry_set_text (GTK_ENTRY(plot_span_entry), value_str);
	sprintf(value_str,"%d", theState->use_rho);
	gtk_entry_set_text (GTK_ENTRY(use_rho_entry), value_str);

	gtk_widget_show (label1);
	gtk_widget_show (label2);
	gtk_widget_show (label3);
	gtk_widget_show (label4);
	gtk_widget_show (label5);
	gtk_widget_show (label6);
	gtk_widget_show (label7);
	gtk_widget_show (label8);
	gtk_widget_show (label9);
	gtk_widget_show (label10);
	gtk_widget_show (label11);
	gtk_widget_show (label12);
	gtk_widget_show (label13);
	gtk_widget_show (label14);
	gtk_widget_show (label15);

	gtk_widget_show (beta_entry);
	gtk_widget_show (c_entry);
	gtk_widget_show (delta_entry);
	gtk_widget_show (eclipse_entry);
	gtk_widget_show (inf_entry);
	gtk_widget_show (p_entry);
	gtk_widget_show (r_entry);
	gtk_widget_show (rho_entry);
	gtk_widget_show (rinf_entry);
	gtk_widget_show (theta_entry);
	gtk_widget_show (latent_inf_entry);
	gtk_widget_show (eclipse_entry);

	/* Parameters frame */
	frame = gtk_frame_new( "Time Plot Variables (select 1 or 2)" );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0);

	table = gtk_table_new(3, 2, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
	gtk_widget_show (table);

	plot_ve_button = gtk_check_button_new_with_label( "Plot Ve vs time" );
	g_signal_connect( GTK_OBJECT( plot_ve_button ), "toggled", G_CALLBACK(setCheckVe), NULL);

	plot_vi_button = gtk_check_button_new_with_label( "Plot Vi vs time" );
	g_signal_connect( GTK_OBJECT( plot_vi_button ), "toggled", G_CALLBACK(setCheckVi), NULL);

	plot_inf_button = gtk_check_button_new_with_label( "Plot Inf vs time" );
	g_signal_connect( GTK_OBJECT( plot_inf_button ), "toggled", G_CALLBACK(setCheckInf), NULL);

	plot_cd8_button = gtk_check_button_new_with_label( "Plot CD8s vs time" );
	g_signal_connect( GTK_OBJECT( plot_cd8_button ), "toggled", G_CALLBACK(setCheckCd8), NULL);

	plot_ACV_button = gtk_check_button_new_with_label( "Plot ACV vs time" );
	g_signal_connect( GTK_OBJECT( plot_ACV_button ), "toggled", G_CALLBACK(setCheckACV), NULL);

	gtk_table_attach_defaults(GTK_TABLE(table), plot_ve_button, 0, 1, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), plot_vi_button, 1, 2, 0, 1);
	gtk_table_attach_defaults(GTK_TABLE(table), plot_inf_button, 0, 1, 1, 2);
	gtk_table_attach_defaults(GTK_TABLE(table), plot_cd8_button, 1, 2, 1, 2);
	gtk_table_attach_defaults(GTK_TABLE(table), plot_ACV_button, 2, 3, 1, 2);

	/* Capture Option frame */
	frame = gtk_frame_new( "Capture Options" );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0);
	svbox = gtk_vbox_new( FALSE, 10 );
	gtk_container_set_border_width( GTK_CONTAINER( svbox ), 5 );
	gtk_container_add( GTK_CONTAINER( frame ), svbox );

	/* File name */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Filename:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	scapture.movie_name = gtk_entry_new();
	g_signal_connect( G_OBJECT(scapture.movie_name),"changed",G_CALLBACK(file_name_changed), (void *)moviebutton );
	gtk_box_pack_start( GTK_BOX( hbox ), scapture.movie_name, TRUE, TRUE, 0 );
	button = gtk_button_new_with_label( "Browse..." );
	g_signal_connect( G_OBJECT( button ), "clicked",G_CALLBACK(open_file_selection), (void *)&get_movie_name);
	gtk_box_pack_start( GTK_BOX( hbox ), button, FALSE, TRUE, 0 );

	/* Capture output format */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Video CODEC:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	scapture.codec = gtk_combo_box_new_text();
	gtk_box_pack_start( GTK_BOX( hbox ), scapture.codec, FALSE, FALSE, 0);

	gtk_combo_box_append_text( GTK_COMBO_BOX(scapture.codec), "OGG Media" );
	gtk_combo_box_append_text( GTK_COMBO_BOX(scapture.codec), "Uncompressed" );
	gtk_combo_box_append_text( GTK_COMBO_BOX(scapture.codec), "WindowsMedia" );
	gtk_combo_box_append_text( GTK_COMBO_BOX(scapture.codec), "QuickTime" );
	gtk_combo_box_set_active( GTK_COMBO_BOX(scapture.codec), 0 );

	/* Frame rate */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Playback frame rate:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	spinbutton_adj = gtk_adjustment_new( 10, 1, INT_MAX, 1, 10, 0 );
	scapture.framerate = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), scapture.framerate, FALSE, FALSE, 0);
	label = gtk_label_new( "fps" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);

	/* Bitrate */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Encoder bitrate:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	spinbutton_adj = gtk_adjustment_new( 600, 1, INT_MAX, 1, 10, 0 );
	scapture.bitrate = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), scapture.bitrate, FALSE, FALSE, 0);
	label = gtk_label_new( "kbps" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);

	/* Capture Option frame */
	frame = gtk_frame_new( "GUI Options" );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0);
	svbox = gtk_vbox_new( FALSE, 10 );
	gtk_container_set_border_width( GTK_CONTAINER( svbox ), 5 );
	gtk_container_add( GTK_CONTAINER( frame ), svbox );

	/* Font name */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Font:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	font_name_widget = gtk_entry_new();
	gtk_box_pack_start( GTK_BOX( hbox ), font_name_widget, TRUE, TRUE, 0 );
	gtk_entry_set_text (GTK_ENTRY(font_name_widget), "Courier bold");

	/* Font size */
	label = gtk_label_new( "Font size:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	spinbutton_adj = gtk_adjustment_new( 12, 8, 24, 2, 2, 0 );



	font_size_widget = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), font_size_widget, FALSE, FALSE, 0);

	button = gtk_button_new_with_label( "Load..." );
	gtk_box_pack_start( GTK_BOX( hbox ), button, FALSE, TRUE, 0 );
	g_signal_connect( G_OBJECT( button ), "clicked",G_CALLBACK(font_changed), (void *)button);

	/* Background RGB */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	rgb_box = gtk_event_box_new( );
	label = gtk_label_new( "Background RGB:" );
	gtk_container_add (GTK_CONTAINER (rgb_box), label);
	gtk_box_pack_start( GTK_BOX( hbox ), rgb_box, FALSE, FALSE, 0);

	spinbutton_adj = gtk_adjustment_new( 255, 0, MAX_RGB_VAL, 1, 10, 0 );
	background_r = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), background_r, FALSE, FALSE, 0);
	label = gtk_label_new( "R" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	g_signal_connect( G_OBJECT(background_r),"changed",G_CALLBACK(background_changed), (void *)NULL );

	spinbutton_adj = gtk_adjustment_new( 255, 0, MAX_RGB_VAL, 1, 10, 0 );
	background_g = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), background_g, FALSE, FALSE, 0);
	label = gtk_label_new( "G" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	g_signal_connect( G_OBJECT(background_g),"changed",G_CALLBACK(background_changed), (void *)NULL );

	spinbutton_adj = gtk_adjustment_new( 255, 0, MAX_RGB_VAL, 1, 10, 0 );
	background_b = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), background_b, FALSE, FALSE, 0);
	label = gtk_label_new( "B" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	g_signal_connect( G_OBJECT(background_b),"changed",G_CALLBACK(background_changed), (void *)NULL );

	GdkColor init_color;

	init_color.red = 65535;
	init_color.green = 65535;
	init_color.blue = 65535;

	gtk_widget_modify_bg(rgb_box,GTK_STATE_NORMAL, &init_color);

	gtk_widget_show_all( main_window );
	gtk_widget_set_sensitive( GTK_WIDGET(backup_slider), TRUE );
	gtk_widget_hide( GTK_WIDGET(backup_slider) );

	/* since the button is present, start this off! */
	theState->writeOn = 0;

	/* these set here since widgets realized by now */
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotLines), theState->plotStyle1 != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (drawVe), theState->plotOpt1 != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (drawInf), theState->plotOpt2 != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (drawCD8s), theState->plotOpt3 != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (drawR0s), theState->plotOpt4 != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (drawVi), theState->plotOpt6 != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (drawNums), theState->plotOpt7 != 0);

	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotColors), theState->plotColor != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotLogs), theState->plotLogs != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotRegions), theState->plotRegions != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (writeEnable), theState->writeOn != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (scrollAxes), theState->scrollAxes != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_ve_button), theState->plotVe != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_vi_button), theState->plotVi != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_inf_button), theState->plotInf != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_cd8_button), theState->plotCd8 != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_ACV_button), theState->plotACV != 0);
}

/**************************************************************************
 *  * The following section contains utility function definitions.
 *   **************************************************************************/

/***
 *  *** Configure the OpenGL framebuffer.
 *   ***/
static GdkGLConfig *
configure_gl (void)
{
  GdkGLConfig *glconfig;

  /* Try double-buffered visual */
  glconfig = gdk_gl_config_new_by_mode ((GdkGLConfigMode)(GDK_GL_MODE_RGB    |
                                        GDK_GL_MODE_DOUBLE));
  if (glconfig == NULL)
    {
      g_print ("\n*** Cannot find the double-buffered visual.\n");
      g_print ("\n*** Trying single-buffered visual.\n");

      /* Try single-buffered visual */
      glconfig = gdk_gl_config_new_by_mode ((GdkGLConfigMode)(GDK_GL_MODE_RGB));
      if (glconfig == NULL)
        {
          g_print ("*** No appropriate OpenGL-capable visual found.\n");
          exit (1);
        }
    }

  return glconfig;
}


#endif

void renderText(double x, double y, string instring)
{
  /*
   * Show font description string.
   */
#ifndef NO_GUI
  glRasterPos2f ((float)x, (float)y);
  glListBase (font_list_base);
  glCallLists (instring.length(), GL_UNSIGNED_BYTE, instring.c_str());
#else
  glRasterPos2f ((float)x, (float)y);
  float rpos[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION ,rpos);

  print (our_font, rpos[0], rpos[1], instring.c_str());
#endif
}

void renderGraphText(double x, double y, string instring)
{
  /*
   * Show font description string.
   */
#ifndef NO_GUI
  glRasterPos2f ((float)x, (float)y);
  glListBase (font_list_base);
  glCallLists (instring.length(), GL_UNSIGNED_BYTE, instring.c_str());
#else
  glRasterPos2f ((float)x, (float)y);
  float rpos[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION ,rpos);

  print (our_font, rpos[0]+40., rpos[1]+40., instring.c_str());
#endif
}

void draw_filled_box(GLfloat theColor[],double x, double y, double w, double h)
{

  glBegin(GL_POLYGON);
  glColor4fv(theColor);
  glVertex2d(x,y);
  glVertex2d(x,y+h);
  glVertex2d(x,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y);
  glVertex2d(x+w,y);
  glVertex2d(x,y);
  glEnd();
  // draw outline
  glBegin(GL_LINES);
  glColor4fv(Black);
  glVertex2d(x,y);
  glVertex2d(x,y+h);
  glVertex2d(x,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y);
  glVertex2d(x+w,y);
  glVertex2d(x,y);
  glEnd();
}
void draw_outline()
{

  // draw outline
  glBegin(GL_LINES);
  glColor4fv(Black);
  glVertex2d(-60.,-60.0);
  glVertex2d(-60.,60.0);
  glVertex2d(-60.,60.0);
  glVertex2d(60.,60.0);
  glVertex2d(60.,60.0);
  glVertex2d(60.,-60.0);
  glVertex2d(60.,-60.0);
  glVertex2d(-60.,-60.0);
  glEnd();
}
void draw_axes()
{
  int i;

  int axis1=0;

  // displayed interval
  double start_time;
  double end_time;


  // draw X axis line
  glBegin(GL_LINES);
  glColor4fv(Black);
  glVertex2d(ORIGIN_X,ORIGIN_Y);
  glVertex2d(MAX_X_COORD,ORIGIN_Y);

  // draw Y axis line - 2
  glVertex2d(ORIGIN_X,ORIGIN_Y);
  glVertex2d(ORIGIN_X,MAX_Y_COORD);
  glVertex2d(MAX_X_COORD,ORIGIN_Y);
  glVertex2d(MAX_X_COORD,MAX_Y_COORD);
  glEnd();

  // draw 9 horiz division lines (1/10th of X max apart)
  for (i=1; i < 10; i++) 
  {
      glBegin(GL_LINES);
      glColor4fv(Gray);
      glVertex2d(ORIGIN_X,ORIGIN_Y+(int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)));
      glVertex2d(MAX_X_COORD,ORIGIN_Y+(int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)));
      glEnd();
  }

  glColor4fv(Black);
  renderGraphText(10,-10, "Time (days)");

  // are we showing the whole run?
  if (theState->plot_span == 0.0)
  {
	start_time = 0.0;
	end_time = (double)theState->max_time;
  }
  // if not, which interval are we showing?
  else
  {
	if (theState->scrollAxes)
	{
	    //always show last n days where n=theState->plot_span
	    //
	    if (theState->time+theState->hex_time_bias <= theState->plot_span)
	    {
		start_time = 0.0;
	    }
	    else
	    {
		start_time = (theState->time+theState->hex_time_bias)-theState->plot_span;
		start_time=MAX(0.0, start_time+theState->plot_bias);
	    }

	    end_time = start_time + theState->plot_span;
	}
	else
	{
	    //divide current time by interval to get start time
	    //
	    if (theState->time+theState->hex_time_bias > theState->plot_span)
	    {
		start_time = theState->plot_span * (double)((int)((theState->time+theState->hex_time_bias)/theState->plot_span));
		start_time=MAX(0.0, start_time+theState->plot_bias);
	    }
	    else
	    {
		start_time = 0.0;
	    }

	    end_time = start_time + theState->plot_span;
	}
  }
  double tick_interval = (theState->plot_span)/5.0;
  double first_tick = (double)(((int)(start_time / tick_interval))*tick_interval);
  if (start_time > first_tick)
	first_tick+=tick_interval;

  for (int i=0; i < 6; i++) {
      char val_str[100];
      double xcoord;

      /* determine the next 'tick' based on start time */
      // sprintf(val_str,"%g",(double)((int)(start_time+0.5)) + ((end_time - start_time)/5.0)*(i));
      sprintf(val_str,"%g",first_tick+i*tick_interval);
      xcoord = 16*(first_tick - start_time)/tick_interval+16*i;
      if (xcoord >= ORIGIN_X && xcoord <= MAX_X_COORD)
      {
	  glColor4fv(Black);
	  renderGraphText(xcoord, -5.0, val_str);
	  if (xcoord > ORIGIN_X && xcoord < MAX_X_COORD)
	  {
	      glBegin(GL_LINES);
	      glColor4fv(Gray);
	      glVertex2d(xcoord,ORIGIN_Y);
	      glVertex2d(xcoord,ORIGIN_Y+2);
	      glEnd();
	  }
      }
  }
	
  // Take care of left axis labeling and title
  if (theState->plotVe)
  {
      glColor4fv(DarkRed);
      axis1=VE_GRAPH;
      renderGraphText(ORIGIN_X,90, string("log 10 cell free HHV8 DNA copies/mL"));
  }
  else if (theState->plotVi)
  {
      glColor4fv(Orange);
      axis1=VI_GRAPH;
      renderGraphText(ORIGIN_X,90, string("log 10 cell assoc HHV8 DNA copies/mL"));
  }
  else if (theState->plotInf)
  {
      glColor4fv(DarkGreen);
      axis1=INF_GRAPH;
      renderGraphText(ORIGIN_X,90, string("log 10 Infected Cells/mL"));
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",((double)((int)(log10((double)theState->max_inf)+0.5))/10.0)*i);
	  renderGraphText(-18, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }

  else if (theState->plotCd8)
  {
      glColor4fv(Black);
      axis1=CD8_GRAPH;
      renderGraphText(ORIGIN_X, 90.0, string("CD8s - Total TCells"));
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",(theState->max_cd8s/10.0)*i);
	  renderGraphText(-18, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }

  else if (theState->plotACV)
  {
      glColor4fv(DarkBlue);
      axis1=ACV_GRAPH;
      renderGraphText(ORIGIN_X, 90.0, string("Drug / IC50"));
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",((theState->max_ACV/theState->IC50_init)/10.0)*i);
	  renderGraphText(-18, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }
  if (theState->plotVi || theState->plotVe)
  {
      if (theState->plotVe)
	  glColor4fv(DarkRed);
      else
	  glColor4fv(Orange);
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",((double)((int)(log10((double)theState->max_vl)+0.5))/10.0)*i);
	  renderGraphText(-18, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }
  // Take care of right axis labeling and title
  if (theState->plotVi && axis1 !=VI_GRAPH)
  {
      glColor4fv(Orange);
      renderGraphText(ORIGIN_X,82, string("log 10 cell assoc HHV8 DNA copies/mL"));
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",((double)((int)(log10((double)theState->max_vl)+0.5))/10.0)*i);
	  renderGraphText(MAX_X_COORD, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }
  else if (theState->plotInf && axis1 !=INF_GRAPH)
  {
      glColor4fv(DarkGreen);
      renderGraphText(ORIGIN_X,82, string("log 10 Infected Cells/mL"));
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",((double)((int)(log10((double)theState->max_inf)+0.5))/10.0)*i);
	  renderGraphText(MAX_X_COORD, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }

  else if (theState->plotCd8 && axis1 !=CD8_GRAPH)
  {
      glColor4fv(Black);
      renderGraphText(ORIGIN_X, 82.0, string("CD8s - Total TCells"));
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",(theState->max_cd8s/10.0)*i);
	  renderGraphText(MAX_X_COORD, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }

  else if (theState->plotACV && axis1 !=ACV_GRAPH)
  {
      glColor4fv(DarkBlue);
      renderGraphText(ORIGIN_X, 82.0, string("Drug / IC50"));
      for (i=1; i < 10; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",((theState->max_ACV/theState->IC50_init)/10.0)*i);
	  renderGraphText(MAX_X_COORD, (int)(i*(0.1)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
  }
}

//! [7]
void draw_graph()
{
  double timeCoord, timeCoordp;

  double vetCoord;
  double vitCoord;
  double veCoords[MAX_HEXCELLS];
  double viCoords[MAX_HEXCELLS];

  double infCoord;
  double cd8Coord;
  double ACVCoord;

  // displayed interval (based on index of samples)
  double start_time;
  double end_time;
  double time_spread;

  int start_sample;
  int end_sample;
  int curr_sample;

  int k;
  int traces=0;

  // No points yet? why bother...
  if (theState->points == NULL) 
	return;

  start_sample=0;
  end_sample=theState->points->valid;
  curr_sample=theState->points->valid;
  time_spread = theState->max_time;

  // are we showing the whole run?
  // if not, which interval are we showing?
  if (theState->plot_span > 0.0)
  {
        time_spread = theState->plot_span;
	if (theState->scrollAxes)
	{
	    //always show last n days where n=theState->plot_span
	    //
	    if (theState->time+theState->hex_time_bias <= theState->plot_span)
	    {
		start_time = 0.0;
		start_sample = 0;
	    }
	    else
	    {
		start_time = (theState->time+theState->hex_time_bias)-theState->plot_span;
		start_time=MAX(0.0, start_time+theState->plot_bias);
		for (int j=1; j < theState->points->valid; j++) 
		{
		    if (theState->points->time[j] >= start_time)
		    {
			start_sample = j;
			break;
		    }
		}
	    }
	}
	else
	{
	    //divide current time by interval to get start time
	    //
	    if (theState->time+theState->hex_time_bias > theState->plot_span)
	    {
		start_time = theState->plot_span * (double)((int)((theState->time+theState->hex_time_bias)/theState->plot_span));
		start_time=MAX(0.0, start_time+theState->plot_bias);
		for (int j=1; j < theState->points->valid; j++) 
		{
		    if (theState->points->time[j] >= start_time)
		    {
			start_sample = j;
			break;
		    }
		}
	    }
	    else
	    {
		start_time = 0.0;
		start_sample = 0;
	    }
	}

	end_time = MIN(start_time + theState->plot_span,theState->time);

	for (int j=1; j < theState->points->valid; j++) 
	{
	    if (theState->points->time[j] > theState->time+theState->hex_time_bias)
	    {
		curr_sample = j;
		break;
	    }
	}

	for (int j=1; j < theState->points->valid; j++) 
	{
	    if (theState->points->time[j] >= end_time)
	    {
		end_sample = j;
		break;
	    }
	}
  }
  glBegin(GL_LINES);
  if (theState->points != NULL) {
	double val;
	val = (theState->points->vet[start_sample] > 0)?log10((double)(theState->points->vet[start_sample])):0;
	vetCoord = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;

	for (k=0; k < theState->Regions;k++)
	{
	    val = (theState->points->ve[k][start_sample] > 0)?log10((double)(theState->points->ve[k][start_sample])):0;
	    veCoords[k] = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;
	}

	val = (theState->points->vit[start_sample] > 0)?log10((double)(theState->points->vit[start_sample])):0;
	vitCoord = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;

	for (k=0; k < theState->Regions;k++)
	{
	    val = (theState->points->vi[k][start_sample] > 0)?log10((double)(theState->points->vi[k][start_sample])):0;
	    viCoords[k] = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;
	}

	val = (theState->points->inf[start_sample] > 0)?log10((double)(theState->points->inf[start_sample])):0;
	infCoord = (val/log10((double)theState->max_inf)) *MAX_Y_COORD;

	val = (double)(theState->points->cd8cells[start_sample]);
	cd8Coord = (val/(double)theState->max_cd8s) *MAX_Y_COORD;

	val = (double)(theState->points->ACV[start_sample]);
	ACVCoord = (val/(double)theState->max_ACV) *MAX_Y_COORD;

	/* plot trace 1 first */
	//for (int graphs=0; graphs < 2; graphs++) {
	    timeCoordp = ((double)(theState->points->time[start_sample]-start_time)/time_spread)
			*MAX_X_COORD;
	for (int i = start_sample+1; i < end_sample; i++) {
          timeCoord = ((double)(theState->points->time[i]-start_time)/time_spread) *MAX_X_COORD;

	  if (theState->plotVe)
	  {
	      traces=VE_GRAPH;
	      if (theState->plotRegions)
	      {
		  for (k=0; k < theState->Regions;k++)
		  {
		      if (theState->points->color[k][i] > 0 && theState->points->color[k][i]  < MAX_PLAQ_COLORS)
			  glColor4fv(pColors[theState->points->color[k][i]]);
		      else
			  glColor4fv(Black);
		      glVertex2d(timeCoordp,veCoords[k]);
		      val = (theState->points->ve[k][i] > 0)?log10((double)(theState->points->ve[k][i])):0;
		      veCoords[k] = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;
		      glVertex2d(timeCoord,veCoords[k]);
		  }
	      }
	      else
	      {
		  /* write dark red line from last viral load reading to current one */
		  glColor4fv(DarkRed);
		  glVertex2d(timeCoordp,vetCoord);
		  val = (theState->points->vet[i] > 0)?log10((double)(theState->points->vet[i])):0;
		  vetCoord = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;
		  glVertex2d(timeCoord,vetCoord);
	      }

	  }
	  if (theState->plotVi && traces != VI_GRAPH)
	  {
	      traces=VI_GRAPH;
	      if (theState->plotRegions)
	      {
		  for (k=0; k < theState->Regions;k++)
		  {
		      if (theState->points->color[k][i] > 0 && theState->points->color[k][i]  < MAX_PLAQ_COLORS)
			  glColor4fv(pColors[theState->points->color[k][i]]);
		      else
			  glColor4fv(Black);
		      glVertex2d(timeCoordp,viCoords[k]);
		      val = (theState->points->vi[k][i] > 0)?log10((double)(theState->points->vi[k][i])):0;
		      viCoords[k] = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;
		      glVertex2d(timeCoord,viCoords[k]);
		  }
	      }
	      else
	      {
		  /* write dark red line from last viral load reading to current one */
		  glColor4fv(Orange);
		  glVertex2d(timeCoordp,vitCoord);
		  val = (theState->points->vit[i] > 0)?log10((double)(theState->points->vit[i])):0;
		  vitCoord = (val/log10((double)theState->max_vl)) *MAX_Y_COORD;
		  glVertex2d(timeCoord,vitCoord);
	      }
	  }

	  if (theState->plotInf && traces != INF_GRAPH)
	  {
	      traces=INF_GRAPH;
	      /* write red line from last infected cell reading to current one */
	      glColor4fv(DarkGreen);
	      glVertex2d(timeCoordp,infCoord);
	      timeCoord = ((double)(theState->points->time[i]-start_time)/time_spread)
			    *MAX_X_COORD;
	      val = (theState->points->inf[i] > 0)?log10((double)(theState->points->inf[i])):0;
	      infCoord = (val/log10((double)theState->max_inf)) *MAX_Y_COORD;
	      glVertex2d(timeCoord,infCoord);
	  }

	  if (theState->plotCd8 && traces != CD8_GRAPH)
	  {
	      traces=CD8_GRAPH;
	      /* write black line from last cd8 reading to current one */
	      glColor4fv(Black);
	      glVertex2d(timeCoordp,cd8Coord);
	      timeCoord = ((double)(theState->points->time[i]-start_time)/time_spread)
			    *MAX_X_COORD;
	      val = theState->points->cd8cells[i];
	      cd8Coord = (val/(double)theState->max_cd8s) *MAX_Y_COORD;
	      glVertex2d(timeCoord,cd8Coord);
	  }

	  if (theState->plotACV && traces != ACV_GRAPH)
	  {
	      traces=ACV_GRAPH;
	      /* write black line from last cd8 reading to current one */
	      glColor4fv(DarkBlue);
	      glVertex2d(timeCoordp,ACVCoord);
	      timeCoord = ((double)(theState->points->time[i]-start_time)/time_spread)
			    *MAX_X_COORD;
	      val = theState->points->ACV[i];
	      ACVCoord = (val/(double)theState->max_ACV) *MAX_Y_COORD;
	      glVertex2d(timeCoord,ACVCoord);
	  }
	  timeCoordp = timeCoord;
	}
	//}
  }
  glEnd();
  if (curr_sample <= end_sample && curr_sample >= start_sample)
  {
      timeCoord = ((double)(theState->points->time[curr_sample]-start_time)/time_spread)
			    *MAX_X_COORD;
      glBegin(GL_TRIANGLES);
      glColor4fv(Black);
      glVertex2d(timeCoord,2);
      glVertex2d(timeCoord-1,0);
      glVertex2d(timeCoord+1,0);
      glEnd();
  }
  glFlush();
}

void draw_cells(int option)
{
  double xfactor;
  double yfactor;
  double factor;

  double total;

  xfactor=50. / hexcell::max_x_vertex;
  yfactor=50. / hexcell::max_y_vertex;

  factor = MIN(xfactor,yfactor);

  if (option == 1 || option == 6)
  {
    glColor4fv(Black);
    if (theState->plotColor)
    {
      if (option == 1)
	  renderText(-30,50, string("log 10 cell free HHV8 DNA copies/mL"));
      else
	  renderText(-30,50, string("log 10 cell assoc HHV8 DNA copies/mL"));

      if (theState->plotRegions)
      {
	  /* print a legend using all pColors */
	  /* a heatmap legend (10 shades of each pColor) */
	  for (int i=0; i < 10; i++) 
	  {
	    double val = 0.0+i*1.0;
	    char val_str[100];
	    sprintf(val_str,"%g",val);

	    double scaleFactor = val/MAX_VL_LOG;

	    glColor4fv(Black);
	    renderText(-54,-58+10*i, val_str);

	    for (int j=1; j < 7;j++)
	    {
		GLfloat thisColor[4]={
			1-(1-(pColors[j][0]))*scaleFactor,
			1-(1-(pColors[j][1]))*scaleFactor,
			1-(1-(pColors[j][2]))*scaleFactor,1.0};
		draw_filled_box(thisColor,-58.0+2.0*(j-1), -54.+10.*i, 2.0, 5.0);
	    }
	  }
      } 
      else
      {
	  /* print a heatmap legend (10 shades of red or blue) */
	  for (int i=0; i < 10; i++) 
	  {
	    double val = i;
	    char val_str[100];
	    sprintf(val_str,"%g",val);

	    double scaleFactor = val/MAX_VL_LOG;

	    glColor4fv(Black);
	    renderText(-54,-58+10*i, val_str);

	    if (option == 1)
	    {
		GLfloat thisColor[4]={1-(1-(DarkRed[0]))*scaleFactor,1-(1-(DarkRed[1]))*scaleFactor,1-(1-(DarkRed[2]))*scaleFactor,1};
		draw_filled_box(thisColor,-58.0, -54.+10.*i, 10.0, 5.0);
	    }
	    else
	    {
		GLfloat thisColor[4]={1-(1-(DarkBlue[0]))*scaleFactor,1-(1-(DarkBlue[1]))*scaleFactor,1-(1-(DarkBlue[2]))*scaleFactor,1};
		draw_filled_box(thisColor,-58.0, -54.+10.*i, 10.0, 5.0);
	    }
	  }
      }
    }
    else if (theState->plotLogs)
    {
      if (option == 1)
	  renderText(-30,50, string("log 10 cell free HHV8 DNA copies/mL"));
      else
	  renderText(-30,50, string("log 10 cell assoc HHV8 DNA copies/mL"));
    }
    else
    {
      if (option == 1)
	  renderText(-30,50, string("cell free HHV8 DNA copies/mL"));
      else
	  renderText(-30,50, string("cell assoc HHV8 DNA copies/mL"));
    }
  }
  else if (option == 2)
  {
    glColor4fv(Black);
    if (theState->plotColor)
    {
      renderText(-20,50, string("log 10 Infected Cells"));

      /* print a heatmap legend (6 shades of green) */
      for (int i=0; i < 7; i++) 
      {
	double val = 0.0+i;
	char val_str[100];
	sprintf(val_str,"%g",val);

	double scaleFactor = val/MAX_INF_LOG;

	glColor4fv(Black);
	renderText(-54,-38+10*i, val_str);

	GLfloat thisColor[4]={1-(1-(DarkGreen[0]))*scaleFactor,1-(1-(DarkGreen[1]))*scaleFactor,1-(1-(DarkGreen[2]))*scaleFactor,1};
	draw_filled_box(thisColor,-58.0, -34.+10.*i, 10.0, 5.0);
      }
    }
    else if (theState->plotLogs)
      renderText(-20,50, string("log 10 Infected Cells"));
    else
      renderText(-10,50, string("Infected Cells"));
  }
  else if (option == 3)
  {
    glColor4fv(DarkBlue);
    if (theState->plotColor)
    {
      renderText(-20,50, string("CD8+ T-cells/region" ));

      /* print a heatmap legend (6 shades of black) */
      for (int i=0; i < 7; i++) 
      {
	double val = MIN_CD8+i*1000.;
	char val_str[100];
	sprintf(val_str,"%g",val);

	double scaleFactor = (val-MIN_CD8)/(MAX_CD8-MIN_CD8);

	glColor4fv(DarkBlue);
	renderText(-54,-38+10*i, val_str);

	GLfloat thisColor[4]={1-(1-(DarkBlue[0]))*scaleFactor,1-(1-(DarkBlue[1]))*scaleFactor,1-(1-(DarkBlue[2]))*scaleFactor,1};
	draw_filled_box(thisColor,-58.0, -34.+10.*i, 10.0, 5.0);
      }
    }
    else if (theState->plotLogs)
      renderText(-30,50, string("log 10 CD8+ T-cells/region" ));
    else
      renderText(-20,50, string("CD8+ T-cells/region" ));
  }
  else if (option == 4)
  {
    glColor4fv(Black);
    if (theState->plotColor)
    {
      renderText(-20,50, string("log 10 Reproductive number"  ));

      /* print a heatmap legend (8 shades of purple) */
      for (int i=0; i < 12; i++) 
      {
	double val = -1.0+i*0.25;
	char val_str[100];
	sprintf(val_str,"%g",val);

	double scaleFactor;
	GLfloat thisColor[4];

	glColor4fv(Black);
	renderText(-54,-58+10*i, val_str);

	if (val < 0)
	{
	    scaleFactor = -val;
	    thisColor[0]=(1-(1-(Black[0]))*scaleFactor);
	    thisColor[1]=(1-(1-(Black[1]))*scaleFactor);
	    thisColor[2]=(1-(1-(Black[2]))*scaleFactor);
	    thisColor[3]=1;
	}
	else
	{
	    scaleFactor = val/2.0;
	    thisColor[0]=(1-(1-(DarkGreen[0]))*scaleFactor);
	    thisColor[1]=(1-(1-(DarkGreen[1]))*scaleFactor);
	    thisColor[2]=(1-(1-(DarkGreen[2]))*scaleFactor);
	    thisColor[3]=1;
	}
	draw_filled_box(thisColor,-58.0, -54.+10.*i, 10.0, 5.0);
      }
    }
    else
      renderText(-20,50, string("log 10 Reproductive number" ));
  }
  else if (option == 5)
  {
  }
  else if (option == 7)
  {
    glColor4fv(Black);
    renderText(-20,50, string("Region numbering"));
  }

  if (!drawnTime)
  {
	  glColor4fv(Black);
	  char time_str[100];
	  sprintf(time_str,"Time = %g",theState->time+theState->hex_time_bias);
	  renderText(-10,-50, time_str);
	  drawnTime = 1;
  }
  glPushMatrix();
  glTranslated(4.0, 0.0,0.0);
  glScaled(factor, factor,0.0);

  total=0;
  for (int i=0; i < hexcell::num_hex_cells; i++)
    total += theState->cells[i]->graph_cell(theState,option,pColors);

  glPopMatrix();

  glFlush();
}

static bool
draw_routine ( GLfloat width, GLfloat height)
{
  /*** OpenGL BEGIN ***/
  glClear (GL_COLOR_BUFFER_BIT);

  /* Set the foreground colour. */
  glColor3f(0.0,0.0,0.0);

  /* count num active plots (1-6) */

  int plotCnt = 0;
  if (theState->plotStyle1)
      plotCnt++;
  if (theState->cells != NULL)
  {
    if (theState->plotOpt1)
	plotCnt++;
    if (theState->plotOpt2)
	plotCnt++;
    if (theState->plotOpt3)
	plotCnt++;
    if (theState->plotOpt4)
	plotCnt++;
    if (theState->plotOpt5)
	plotCnt++;
    if (theState->plotOpt6)
	plotCnt++;
  }

  int rows, cols;
  switch(plotCnt)
  {
      case 1:
	  cols=1;
	  rows=1;
	  break;
      case 2:
	  cols=2;
	  rows=1;
	  break;
      case 3:
	  cols=3;
	  rows=1;
	  break;
      case 4:
	  cols=2;
	  rows=2;
	  break;
      default:
	  cols=3;
	  rows=2;
	  break;
  }
  double xfactor = 1./(double)cols;;
  double yfactor = 1./(double)rows;;

  double factor = MIN(xfactor,yfactor);

  double xshift_size = 120.;
  double yshift_size = 120.;

  double xshift = 0.;
  double yshift = 0.;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  glPushMatrix();
  glScaled(xfactor, yfactor,0.0);

  if (plotCnt > 3)
    yshift = 60.;

  if (cols > 2)
  {
    xshift = -120.;
  }
  else if (cols > 1)
    xshift = -60.;

  glTranslated(xshift, yshift,0.0);
  plotCnt = 0;
  if (theState->plotStyle1)
  {
	plotCnt++;
	draw_outline();
	glPushMatrix();
	glTranslated(-40.0, -40.0,0.0);
	glPushMatrix();
	draw_axes();
	draw_graph();
	glPopMatrix();
	glPopMatrix();
	glTranslated(xshift_size, 0.0,0.0);
  }
  if (theState->cells != NULL )
  {
    drawnTime = 0;
    if (theState->plotOpt1)
    {
	draw_outline();
	draw_cells(1);
	plotCnt++;
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
    }
    if (theState->plotOpt2)
    {
	draw_outline();
	draw_cells(2);
	plotCnt++;
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
    }
    if (theState->plotOpt3)
    {
	draw_outline();
	draw_cells(3);
	plotCnt++;
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
    }
    if (theState->plotOpt4)
    {
	draw_outline();
	draw_cells(4);
	plotCnt++;
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
    }
    if (theState->plotOpt5)
    {
	draw_outline();
	draw_cells(5);
	plotCnt++;
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
    }
    if (theState->plotOpt6)
    {
	draw_outline();
	draw_cells(6);
	plotCnt++;
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
    }
  }
  glPopMatrix();

  
  glFlush ();

  return true;
}

void usage(char *prog_name)
{

    fprintf(stderr,
	"Usage: %s [-h][-d][-b][-f <input_file>][-c <crit file>][-n <cells>[-r][-s <seed>][-v][-w <write_mask>][-W <crit><weight>]\n",prog_name);
    fprintf(stderr,"\t-h = this help\n");
    fprintf(stderr,"\t-f = optional input file\n");
    fprintf(stderr,"\t-c = optional criteria file\n");
    fprintf(stderr,"\t\tFormat: target lower and upper CIs and mean values for...\n");
    fprintf(stderr,"\t\t\t cumul percent pos swabs (9 bins)\n");
    fprintf(stderr,"\t\t\t peak log VL histograms (8 bins)\n");
    fprintf(stderr,"\t-w <write_mask> = which output (csv) files to generate (Bit-mask 1-13)\n");
    fprintf(stderr,"\t\t1= cumul & regional ve vs time\n");
    fprintf(stderr,"\t\t2= ve current episode episode\n");
    fprintf(stderr,"\t\t3= cumul & regional vi vs time\n");
    fprintf(stderr,"\t\t4= vi current episode episode\n");
    fprintf(stderr,"\t\t5= cumul & regional inf cells vs time\n");
    fprintf(stderr,"\t\t6= inf cells current episode episode\n");
    fprintf(stderr,"\t\t8= regional R0s vs time\n");
    fprintf(stderr,"\t\t9= run stats (per run)\n");
    fprintf(stderr,"\t\t10= transmission info\n");
    fprintf(stderr,"\t\t11= episode stats (per run)\n");
    fprintf(stderr,"\t\t12= more detailed regional info\n");
    fprintf(stderr,"\t\t13= ACV dosing info\n");
    fprintf(stderr,"\t-b = batch mode (no GUI)\n");
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
    char dat_file13[] = "hhv8_sim.dat13";

    char criteria_file[] = "hhv8_sim.crit";
    char *crit_file;

    FILE *out;
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

    batchMode = 0;

    stoch = 1;
    int writeData = 0;
    int writeMask = 0;

    int crit_num = 0;
    double crit_weight = 0;

    double best;
    gsl_vector *params;

    hexcell *the_cells[MAX_HEXCELLS];
    int next_available = 1;

    int num_cells;
    hexcell::num_hex_cells = 0;

    pColors[0] = White;
    pColors[1] = DarkRed;
    pColors[2] = DarkGreen;
    pColors[3] = DarkBlue;
    pColors[4] = Orange;
    pColors[5] = Gold;
    pColors[6] = DarkPurple;

#ifndef NO_GUI
    scapture.pipefd = -1;
#endif

    for(int i = 0; i < argc; i++) {
	cout << "argv[" << i << "] = " << argv[i] << endl;
	if (!strcmp(argv[i],"-b") || !strcmp(argv[i],"-B")) {
	    batchMode = 1;
	}
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
    settings.dataF13 = NULL;

    if (writeData) {

	if (batchMode != 0)
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
	if(((writeMask & (1<<12)) && ((settings.dataF13 = fopen(dat_file13,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file13<<"\n";
	    exit(1);
	}
    }

#ifndef NO_GUI
    if (batchMode == 0)
    {
	GtkWidget *window;

	GdkGLConfig *glconfig;

	/* Initialize GTK. */
	gtk_init (&argc, &argv);

	/* Initialize GtkGLExt. */
	gtk_gl_init (&argc, &argv);

	/* Configure OpenGL framebuffer. */
	//glconfig = configure_gl ();

	/* Create and show the application window. */
	//window = create_window (glconfig);
	// gtk_widget_show (window);

	settings.N_runs = 1;

	make_gui();

	gtk_main ();
    } 
    else 
    {
	params=gsl_vector_alloc(NUM_PARAMS);

	best =fitmodel(params, &settings);
    }
#else 
    {
        bgr = initialize_hidden(DEFAULT_IMAGE_WIDTH, DEFAULT_IMAGE_HEIGHT);

	params=gsl_vector_alloc(NUM_PARAMS);

	our_font.init("Test.ttf", 16);      

	best =fitmodel(params, &settings);

	if (settings.Verbose)
	{
	    fprintf(stderr,"Final score = %lf for the following values...\n\n",best);
	    fprintf(stderr,"beta = %g\n",gsl_vector_get(params,0));
	    fprintf(stderr,"latent_inf = %g\n",gsl_vector_get(params,2));
	    fprintf(stderr,"p = %g\n",gsl_vector_get(params,2));
	    fprintf(stderr,"c = %g\n",gsl_vector_get(params,3));
	    fprintf(stderr,"theta = %g\n",gsl_vector_get(params,4));
	    fprintf(stderr,"inf = %g\n",gsl_vector_get(params,5));
	    fprintf(stderr,"r = %g\n",gsl_vector_get(params,6));
	    fprintf(stderr,"rinf = %g\n",gsl_vector_get(params,7));
	    fprintf(stderr,"delta = %g\n",gsl_vector_get(params,8));
	    fprintf(stderr,"rho = %g\n",gsl_vector_get(params,10));
	    fprintf(stderr,"eclipse = %g\n",gsl_vector_get(params,11));
	}
	free(bgr);
    }
#endif

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
	if(settings.dataF13 != NULL) fclose(settings.dataF13);
    }


}//end of main
