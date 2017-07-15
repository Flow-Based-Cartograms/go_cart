/****************************** Inclusions. **********************************/

#include <stdio.h>
#include <stdlib.h>
#include "cartogram.h"

/**************************** Function prototypes. ***************************/

void count_poly (FILE *in_file);
void count_corn (FILE *in_file);
void read_corn (FILE *in_file);
void make_region (void);

/*****************************************************************************/
/******************* Function to count the number of polygons. ***************/

void count_poly (FILE *in_file)
{
  char line[MAX_STRING_LENGTH];

  n_poly = 0;
  while (fgets(line, MAX_STRING_LENGTH, in_file) != NULL)
    if (line[0] == 'E') n_poly++;
  n_poly--;               /* The .gen file ends with two consecutive "END"s. */

  return;
}

/*****************************************************************************/
/**** Function to count polygon corners. Also determines minimum/maximum  ****/
/**** x-/y-coordinate.                                                    ****/

void count_corn (FILE *in_file)
{
  char line[MAX_STRING_LENGTH];
  double x, y;
  int polyctr = 0;

  n_polycorn = (int*) calloc(n_poly, sizeof(int));
  polycorn = (POINT**) malloc(n_poly * sizeof(POINT*));
  if (fgets(line, MAX_STRING_LENGTH, in_file) == NULL) { /* Skip first line. */
    fprintf(stderr, "ERROR: in_file empty.\n");
    exit(1);
  }
  if (fgets(line, MAX_STRING_LENGTH, in_file) != NULL) {
    sscanf(line, "%lf %lf", &map_minx, &map_miny);
    map_maxx = map_minx;
    map_maxy = map_miny;
    n_polycorn[0] = 1;
  }
  else {
    fprintf(stderr, "ERROR: in_file has only one line?\n");
    exit(1);
  }
  while (fgets(line, MAX_STRING_LENGTH, in_file) != NULL) {
    if (line[0] != 'E') {
      sscanf(line, "%lf %lf", &x, &y);
      map_minx = MIN(map_minx, x);
      map_maxx = MAX(map_maxx, x);
      map_miny = MIN(map_miny, y);
      map_maxy = MAX(map_maxy, y);
      n_polycorn[polyctr]++;
    }
    else {
      if (fgets(line, MAX_STRING_LENGTH, in_file) == NULL) {
      	fprintf(stderr, "ERROR: in_file not in proper format.\n");
      	exit(1);
      }
      polycorn[polyctr] = (POINT*) malloc(n_polycorn[polyctr] * sizeof(POINT));
      polyctr++;
    }
  }

  return;
}

/*****************************************************************************/
/**** Function to read polygon corners. The first and last vertex of each ****/
/**** polygon must be identical.                                          ****/

void read_corn (FILE *in_file)
{
  char line[MAX_STRING_LENGTH];
  int i, id, polyctr=0;

  polygon_id = (int*) malloc(n_poly * sizeof(int));
  if (fgets(line, MAX_STRING_LENGTH, in_file) == NULL) { /* Skip first line. */
    fprintf(stderr, "ERROR: in_file empty.\n");
    exit(1);
  }
  sscanf(line, "%i", &id);
  polygon_id[polyctr] = id;
  i = 0;
  while (fgets(line, MAX_STRING_LENGTH, in_file) != NULL) {
    if (line[0] != 'E') {
      sscanf(line, "%lf %lf",
	     &polycorn[polyctr][i].x, &polycorn[polyctr][i].y);
      i++;
    }
    else {
      /* Are first and last corner identical? */
      if (polycorn[polyctr][0].x !=
    	           polycorn[polyctr][n_polycorn[polyctr]-1].x ||
    	           polycorn[polyctr][0].y !=
    	           polycorn[polyctr][n_polycorn[polyctr]-1].y) {
        fprintf(stderr, "WARNING: %i-th polygon does not close upon itself.\n",
		        polyctr+1);
        fprintf(stderr, "Identifier %i, first point (%f, %f).\n",
		        polygon_id[polyctr], polycorn[polyctr][0].x,
		        polycorn[polyctr][0].y);
        fprintf(stderr, "%i points.\n", n_polycorn[polyctr]);
      }
      if (fgets(line, MAX_STRING_LENGTH, in_file) == NULL) {
	      fprintf(stderr, "ERROR: in_file not in proper format.\n");
	      exit(1);
      }
      sscanf(line, "%i", &id);
      i = 0;
      polyctr++;
      if (polyctr < n_poly)
	      polygon_id[polyctr] = id;
    }
  }

  return;
}

/*****************************************************************************/
/* Function to make regions from polygons. Region IDs in the .gen file must  */
/* be nonnegative.                                                           */

void make_region (void)
{
  BOOLEAN repeat;
  int i, j, last_id, min_id;

  n_reg = 0;                                 /* Count the number of regions. */
  max_id = min_id = polygon_id[0];
  for (j=0; j<n_poly; j++) {       /* -99999 is a special ID. Such polygons  */
    if (polygon_id[j] == -99999)   /* will be assigned to the same region as */
      continue;                    /* the polygon immediately before it.     */
    if (polygon_id[j]>max_id || polygon_id[j]<min_id)
      n_reg++;
    else {
      repeat = FALSE;
      for (i=0; i<j; i++)
	      if (polygon_id[j] == polygon_id[i]) {
      	  repeat = TRUE;
      	  break;
        }
      if (!repeat)
        n_reg++;
    }
    max_id = MAX(max_id, polygon_id[j]);
    min_id = MIN(min_id, polygon_id[j]);
  }
  if (min_id < 0) {
    fprintf(stderr, "ERROR: Negative region identifier %i.\n", min_id);
    exit(1);
  }
  region_id = (int*) malloc(n_reg * sizeof(int));       /* Match region IDs. */
  n_reg = 0;
  max_id = min_id = polygon_id[0];
  for (j=0; j<n_poly; j++) {
    if (polygon_id[j] == -99999)
      continue;
    if (polygon_id[j]>max_id || polygon_id[j]<min_id)
      region_id[n_reg++] = polygon_id[j];
    else {
      repeat = FALSE;
      for (i=0; i<j; i++)
      	if (polygon_id[j] == polygon_id[i]) {
      	  repeat = TRUE;
      	  break;
      	}
      if (!repeat)
	      region_id[n_reg++] = polygon_id[j];
    }
    max_id = MAX(max_id, polygon_id[j]);
    min_id = MIN(min_id, polygon_id[j]);
  }

  /* region_id[i] takes as input C's internal identifier for the region and  */
  /* assumes the value of the ID in the .gen file.                           */
  /* region_id_inv[i] does the opposite. Its input is an ID from the .gen    */
  /* file. Its value is the internal identifier used by C.                   */

  region_id_inv = (int*) malloc((max_id+1) * sizeof(int));
  for (i=0; i<=max_id; i++)                 /* -1 for unused IDs that do not */
    region_id_inv[i] = -1;                  /* appear in the .gen file.      */
  for (i=0; i<n_reg; i++)
    region_id_inv[region_id[i]] = i;
  n_polyinreg = (int*) calloc(n_reg, sizeof(int));  /* Which polygons con-   */
  polyinreg = (int**) malloc(n_reg * sizeof(int*)); /* tribute to which      */
  last_id = polygon_id[0];                          /* region?               */
  for (j=0; j<n_poly; j++) {
    if (polygon_id[j] != -99999) {
      n_polyinreg[region_id_inv[polygon_id[j]]]++;
      last_id = polygon_id[j];
    }
    else
      n_polyinreg[region_id_inv[last_id]]++;
  }
  for (j=0; j<n_reg; j++)
    polyinreg[j] = (int*) malloc(n_polyinreg[j] * sizeof(int));
  for (j=0; j<n_reg; j++)
    n_polyinreg[j] = 0;
  last_id = polygon_id[0];
  for (j=0; j<n_poly; j++) {
    if (polygon_id[j] != -99999) {
      polyinreg[region_id_inv[polygon_id[j]]]
	            [n_polyinreg[region_id_inv[polygon_id[j]]]++] = j;
      last_id = polygon_id[j];
    }
    else
      polyinreg[region_id_inv[last_id]]
	            [n_polyinreg[region_id_inv[last_id]]++] = j;
  }

  return;
}

/*****************************************************************************/
/******************* Function to process map information. ********************/

void read_gen (char *gen_file_name)
{
  FILE *gen_file;

  if ((gen_file = fopen(gen_file_name,"r")) == NULL) {
    fprintf(stderr,"ERROR: Cannot find gen-file.\n");
    exit(1);
  }
  count_poly(gen_file);
  fclose(gen_file);
  gen_file = fopen(gen_file_name,"r");
  count_corn(gen_file);
  fclose(gen_file);
  gen_file = fopen(gen_file_name,"r");
  read_corn(gen_file);
  fclose(gen_file);
  make_region();
  printf("%i polygon(s), %i region(s)\n", n_poly, n_reg);

  return;
}
