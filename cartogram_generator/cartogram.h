/******************************** Inclusions. ********************************/

#include <fftw3.h>

/******************************** Definitions. *******************************/

/* Areas on cartogram differ at most by an absolute relative error of        */
/* MAX_PERMITTED_AREA_ERROR. That is,                                        */
/* |area_on_cartogram / target_area - 1| <= MAX_PERMITTED_AREA_ERROR.        */

#define MAX_PERMITTED_AREA_ERROR (0.01)
#define L (512)            /* Maximum dimension of the FFT lattice is L x L. */

/* The program checks whether the area for each region is smaller than a     */
/* particular threshold. If it is, the region's area will be set to that     */
/* threshold. The threshold is scaled to each region's perimeter and is      */
/* calculated by                                                             */
/* MAX(region_perimeter[i]/total_perimeter * MIN_PERIMETER_FAC, 0.00025)     */
/* This is only applicable when the -n flag is not set.                      */
#define MIN_PERIMETER_FAC (0.025)/* Minimum area threshold determined by this*/

/* If the -n flag is used, areas will be calculated accurately. However, if  */
/* region contains exactly zero population, it will be replaced by           */
/* MIN_POP_FAC times the smallest positive population in any region.         */
/* This is only applicable when the -n flag is set.                          */
#define MIN_POP_FAC (0.2)       /* Replace area 0 by the minimum times this. */
#define PADDING (1.5)          /* Determines space between map and boundary. */
#define BLUR_WIDTH (5e0)  /* Width of Gaussian blur to smoothen the density. */
#define MAX_STRING_LENGTH (1000)
#define PI (3.14159265358979323846264338327950288419716939937510)

/********************************** Macros. **********************************/

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)>(b)) ? (b) : (a))

/*********************************** Types. **********************************/

typedef struct {      /* Useful structure for coordinates in two dimensions. */
  double x;
  double y;
} POINT;
typedef enum {          /* Declares an enumeration data type called BOOLEAN. */
  FALSE,                /* FALSE = 0, TRUE = 1 */
  TRUE
} BOOLEAN; 
typedef struct {
  char name[100];
  int id;
} name_id_pair;

/***************************** Global variables. *****************************/

/* Variables for map. */

extern double *area_err, *cart_area, map_maxx, map_maxy, map_minx, map_miny,
  *target_area, *region_perimeter;
extern int max_id, n_poly, *n_polycorn, *n_polyinreg, n_reg, *polygon_id, *poly_is_hole,
  **polyinreg, *region_id, *region_id_inv, *region_na;
extern POINT **cartcorn, **origcorn, **polycorn, *proj, *proj2;
extern BOOLEAN use_perimeter_threshold;

/* Variables for digitizing the density. */

extern double *rho_ft, *rho_init;
extern fftw_plan plan_fwd;
extern int lx, ly;

/**************************** Function prototypes. ***************************/

void set_inside_values_for_polygon (int region, int n_polycorn,
				    POINT *polycorn, int **inside);
double polygon_area (int ncrns, POINT *polygon);
double polygon_perimeter (int ncrns, POINT *polygon);
BOOLEAN fill_with_density1 (char *map_file_name, char *area_file_name,
			    BOOLEAN inv, BOOLEAN eps);
void fill_with_density2 (void);
void read_map (char *map_file);
void ps_figure (char *ps_name, POINT **corn, POINT *proj, BOOLEAN grat);
double interpol (double x, double y, double *grid, char zero);
void ffb_integrate (void);
void diff_integrate (void);
void project (BOOLEAN proj_graticule);
double max_area_err (double *area_err, double *cart_area, POINT **corn,
		     double *sum_cart_area);
void output_to_gen (BOOLEAN usestd, POINT **corn);
void output_to_geojson (BOOLEAN usestd, POINT **corn, char *map_file_name);
void output_error (void);
void inv_project (void);
void print_usage (char *program_name);
void process_json (char *json_file_name, char *program_name);
