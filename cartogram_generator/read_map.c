/****************************** Inclusions. **********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cjson/cJSON.h>
#include "cartogram.h"

/******************************** Definitions. *******************************/

/* We remove areas less than AREA_THRESHOLD * (area of bounding box). */
#define AREA_THRESHOLD (1e-12)

/**************************** Function prototypes. ***************************/

void read_geojson (FILE *in_file);
void count_poly (FILE *in_file);
void count_corn (FILE *in_file);
void read_corn (FILE *in_file);
void make_region (void);
void remove_tiny_polygons (void);

/*****************************************************************************/
/******************* Function to read geojson. ***************/
void read_geojson (FILE *in_file)
{
  char *strJson = NULL;
  int line_len;
  int MAX_STRING_LENGTH_JSON = 1000000; /* Large in order to read files more quickly */
  char* str = (char *) malloc(MAX_STRING_LENGTH_JSON);
  int strJson_size = 0;
  
  while (fgets(str, MAX_STRING_LENGTH_JSON, in_file) != NULL){
    line_len = strlen(str);
    if(line_len > 0 && str[line_len - 1] == '\n'){
      str[line_len - 1] = '\0';
    }

    if(strJson == NULL){
      strJson_size = line_len + 1;
      strJson = (char *) calloc(line_len + 1, sizeof(char));
      strncpy(strJson, str, strJson_size - strlen(strJson) - 1);
    }else{
      strJson_size = (strlen(strJson)) + line_len + 1;
      strJson = (char *) realloc(strJson, (strlen(strJson)) + line_len + 1);
      strncat(strJson, str, strJson_size - strlen(strJson) - 1);
    }
  }
  free(str);
  cJSON *root = cJSON_Parse(strJson);
  free(strJson);
  cJSON *features = cJSON_GetObjectItemCaseSensitive(root, "features");
  cJSON *bbox = cJSON_GetObjectItemCaseSensitive(root, "bbox"); /* Accesses bbox attribute from geojson */
  n_reg = cJSON_GetArraySize(features); /* Calculates number of regions */
  cJSON *feature0 = cJSON_GetArrayItem(features, 0);
  cJSON *feature_properties_0 = cJSON_GetObjectItemCaseSensitive(feature0, "properties");
  cJSON *feature_cartogramid = cJSON_GetObjectItemCaseSensitive(feature_properties_0, "cartogram_id");
  if(!feature_cartogramid){
    fprintf(stderr, "ERROR: GeoJSON file needs to be processed first.\nPlease process the file using \"cartogram -p [json_file_name]\"\n");
    exit(1);
  }
  if (!bbox){
    fprintf(stderr, "ERROR: GeoJSON file needs to contain bbox values.\nPlease process the file using \"cartogram -p [json_file_name]\"\n");
    exit(1);
  }else{
    map_minx = cJSON_GetArrayItem(bbox, 0)->valuedouble; /* Stores bounding box values */
    map_miny = cJSON_GetArrayItem(bbox, 1)->valuedouble;
    map_maxx = cJSON_GetArrayItem(bbox, 2)->valuedouble;
    map_maxy = cJSON_GetArrayItem(bbox, 3)->valuedouble;
  }
  n_poly = 0;
  // This for loop counts the number of polygons in the geojson file
  cJSON *feature_iterator = NULL;
  cJSON_ArrayForEach(feature_iterator, features){
    cJSON *feature_geometry = cJSON_GetObjectItemCaseSensitive(feature_iterator, "geometry");
    char * type = cJSON_GetObjectItemCaseSensitive(feature_geometry, "type")->valuestring;
    cJSON * feature_geometry_coordinates = cJSON_GetObjectItemCaseSensitive(feature_geometry, "coordinates");    
    if(strcmp(type, "Polygon") == 0){
      cJSON * polygon_array_of_linear_rings = feature_geometry_coordinates;
      n_poly += cJSON_GetArraySize(polygon_array_of_linear_rings);
    }else if(strcmp(type, "MultiPolygon") == 0){
      cJSON * multipolygon_array_of_polygons = feature_geometry_coordinates;
      cJSON * polygon_array_of_linear_rings = NULL;
      cJSON_ArrayForEach(polygon_array_of_linear_rings, multipolygon_array_of_polygons){
        n_poly += cJSON_GetArraySize(polygon_array_of_linear_rings);
      }
    }else{
      fprintf(stderr, "Error: Region contains geometry other than polygons and multipolygons.\n");
      exit(1);
    }
  }

  //This for loop counts the number of polygon corners in the geojson file
  n_polycorn = (int*) calloc(n_poly, sizeof(int));
  polycorn = (POINT**) malloc(n_poly * sizeof(POINT*));
  int polyctr = 0;
  feature_iterator = NULL;
  cJSON_ArrayForEach(feature_iterator, features){
    cJSON *feature_geometry = cJSON_GetObjectItemCaseSensitive(feature_iterator, "geometry");
    char * type = cJSON_GetObjectItemCaseSensitive(feature_geometry, "type")->valuestring;
    cJSON * feature_geometry_coordinates = cJSON_GetObjectItemCaseSensitive(feature_geometry, "coordinates");
    if(strcmp(type, "Polygon") == 0){
      cJSON * polygon_array_of_linear_rings = feature_geometry_coordinates;
      cJSON * linear_ring_array_of_positions = NULL;
      cJSON_ArrayForEach(linear_ring_array_of_positions, polygon_array_of_linear_rings){
        n_polycorn[polyctr] = cJSON_GetArraySize(linear_ring_array_of_positions);
        polycorn[polyctr] = (POINT*) malloc(n_polycorn[polyctr] * sizeof(POINT));
        polyctr++;
      }
    }else if(strcmp(type, "MultiPolygon") == 0){
      cJSON * multipolygon_array_of_polygons = feature_geometry_coordinates;
      cJSON * polygon_array_of_linear_rings = NULL;
      cJSON_ArrayForEach(polygon_array_of_linear_rings, multipolygon_array_of_polygons){
        cJSON * linear_ring_array_of_positions = NULL;
        cJSON_ArrayForEach(linear_ring_array_of_positions, polygon_array_of_linear_rings){
          n_polycorn[polyctr] = cJSON_GetArraySize(linear_ring_array_of_positions);
          polycorn[polyctr] = (POINT*) malloc(n_polycorn[polyctr] * sizeof(POINT));
          polyctr++;
        }
      }
    }else{
      fprintf(stderr, "Error: Region contains geometry other than polygons and multipolygons.\n");
      exit(1);
    }
  }
  polyctr = 0;
  polygon_id = (int*) malloc(n_poly * sizeof(int));
  //This for loop reads in and stores the polygon corners in the geojson file
  feature_iterator = NULL;
  cJSON_ArrayForEach(feature_iterator, features){
    cJSON * feature_properties = cJSON_GetObjectItemCaseSensitive(feature_iterator, "properties");
    char * feature_id = cJSON_GetObjectItemCaseSensitive(feature_properties, "cartogram_id")->valuestring;
    int id;
    sscanf(feature_id, "%d", &id);
    cJSON *feature_geometry = cJSON_GetObjectItemCaseSensitive(feature_iterator, "geometry");
    char * type = cJSON_GetObjectItemCaseSensitive(feature_geometry, "type")->valuestring;
    cJSON * feature_geometry_coordinates = cJSON_GetObjectItemCaseSensitive(feature_geometry, "coordinates");
    if(strcmp(type, "Polygon") == 0){
      cJSON * polygon_array_of_linear_rings = feature_geometry_coordinates;
      cJSON * linear_ring_array_of_positions = NULL;
      cJSON_ArrayForEach(linear_ring_array_of_positions, polygon_array_of_linear_rings){
        cJSON * position = NULL;
        int k = 0;
        cJSON_ArrayForEach(position, linear_ring_array_of_positions){
          polycorn[polyctr][k].x = cJSON_GetArrayItem(position, 0)->valuedouble;
          polycorn[polyctr][k].y = cJSON_GetArrayItem(position, 1)->valuedouble;
          k++;
        }
        /* Are first and last corner identical? */
        if (polycorn[polyctr][0].x !=
          polycorn[polyctr][n_polycorn[polyctr]-1].x ||
          polycorn[polyctr][0].y !=
          polycorn[polyctr][n_polycorn[polyctr]-1].y) {
          fprintf(stderr, "WARNING: %i-th polygon does not close upon itself.\n", polyctr+1);
          fprintf(stderr, "Identifier %i, first point (%f, %f).\n",
          polygon_id[polyctr], polycorn[polyctr][0].x,
          polycorn[polyctr][0].y);
          fprintf(stderr, "%i points.\n", n_polycorn[polyctr]);
        }
        polygon_id[polyctr] = id;
        polyctr++;
      }
    }else if(strcmp(type, "MultiPolygon") == 0){
      cJSON * multipolygon_array_of_polygons = feature_geometry_coordinates;
      cJSON * polygon_array_of_linear_rings = NULL;
      cJSON_ArrayForEach(polygon_array_of_linear_rings, multipolygon_array_of_polygons){
        cJSON * linear_ring_array_of_positions = NULL;
        cJSON_ArrayForEach(linear_ring_array_of_positions, polygon_array_of_linear_rings){
          cJSON * position = NULL;
          int k = 0;
          cJSON_ArrayForEach(position, linear_ring_array_of_positions){
            polycorn[polyctr][k].x = cJSON_GetArrayItem(position, 0)->valuedouble;
            polycorn[polyctr][k].y = cJSON_GetArrayItem(position, 1)->valuedouble;
            k++;
          }
          /* Are first and last corner identical? */	  
          if (polycorn[polyctr][0].x !=
            polycorn[polyctr][n_polycorn[polyctr]-1].x ||
            polycorn[polyctr][0].y !=
            polycorn[polyctr][n_polycorn[polyctr]-1].y) {
            fprintf(stderr, "WARNING: %i-th polygon does not close upon itself.\n", polyctr+1);
            fprintf(stderr, "Identifier %i, first point (%f, %f).\n",
            polygon_id[polyctr], polycorn[polyctr][0].x,
            polycorn[polyctr][0].y);
            fprintf(stderr, "%i points.\n", n_polycorn[polyctr]);
          }
          polygon_id[polyctr] = id;
          polyctr++;
        }
      }
    }else{
      fprintf(stderr, "Error: Region contains geometry other than polygons and multipolygons.\n");
      exit(1);
    }
  }
  cJSON_Delete(root);
  return;
}

/*****************************************************************************/
/**** Function to count number of polygons (for gen file only)            ****/

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
/**** x-/y-coordinate. (for gen file only)                                ****/

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
/**** polygon must be identical.   (for gen file only)                    ****/

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
/* Function to determine polygon area. This is needed to remove polygons     */
/* zero area and determine the average population.                            */
/* The problem in short is to find the area of a polygon whose vertices are  */
/* given. Recall Stokes' theorem in 3d for a vector field v:                 */
/* integral[around closed curve dA]v(x,y,z).ds =                             */
/*                                          integral[over area A]curl(v).dA. */
/* Now let v(x,y,z) = (0,Q(x,y),0) and dA = (0,0,dx*dy). Then                */
/* integral[around closed curve dA]Q(x,y)dy =                                */
/*                                         integral[over area A]dQ/dx*dx*dy. */
/* If Q = x:                                                                 */
/* A = integral[over area A]dx*dy = integral[around closed curve dA]x dy.    */
/* For every edge from (x[i],y[i]) to (x[i+1],y[i+1]) there is a             */
/* parametrization                                                           */
/* (x(t),y(t)) = ((1-t)x[i]+t*x[i+1],(1-t)y[i]+t*y[i+1]), 0<t<1              */
/* so that the path integral along this edge is                              */
/* int[from 0 to 1]{(1-t)x[i]+t*x[i+1]}(y[i+1]-y[i])dt =                     */
/*                                          0.5*(y[i+1]-y[i])*(x[i]+x[i+1]). */
/* Summing over all edges yields:                                            */
/* Area = 0.5*[(x[0]+x[1])(y[1]-y[0]) + (x[1]+x[2])(y[2]-y[1]) + ...         */
/*               ... + (x[n-1]+x[n])(y[n]-y[n-1]) + (x[n]+x[0])(y[0]-y[n])]. */
/* ArcGIS treats a clockwise direction as positive, so that there is an      */
/* additional minus sign.                                                    */

double polygon_area (int ncrns, POINT *polygon)
{
  double area = 0.0;
  int i;

  for (i=0; i<ncrns-1; i++)
    area -=
      0.5 * (polygon[i].x+polygon[i+1].x) * (polygon[i+1].y-polygon[i].y);
  return area -= 0.5 * (polygon[ncrns-1].x+polygon[0].x) *
          (polygon[0].y-polygon[ncrns-1].y);
}

double polygon_perimeter (int ncrns, POINT *polygon) {
  double perimeter = 0.0;
  for (int i=0; i<ncrns - 1; i++){
    perimeter += sqrt((polygon[i+1].x - polygon[i].x) * (polygon[i+1].x - polygon[i].x) +
      (polygon[i+1].y - polygon[i].y) * (polygon[i+1].y - polygon[i].y));
  }
  return perimeter + sqrt((polygon[0].x - polygon[ncrns-1].x) * (polygon[0].x - polygon[ncrns-1].x) +
    (polygon[0].y - polygon[ncrns-1].y) * (polygon[0].y - polygon[ncrns-1].y));
}

/*****************************************************************************/
/********************* Function to remove tiny polygons. *********************/

void remove_tiny_polygons (void) {
  
  /* Find out whether there are any tiny polygons. */
  
  BOOLEAN *poly_has_tiny_area = (BOOLEAN*) calloc(n_poly, sizeof(BOOLEAN));
  for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
    poly_has_tiny_area[poly_indx] =
      (fabs(polygon_area(n_polycorn[poly_indx], polycorn[poly_indx])) <
       AREA_THRESHOLD * (map_maxx - map_minx) * (map_maxy - map_miny));
  }
  int n_non_tiny_poly = 0;
  for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
    n_non_tiny_poly += !(poly_has_tiny_area[poly_indx]);
  }
  if (n_non_tiny_poly < n_poly) {
    fprintf(stderr, "Removing tiny polygons.\n");
    
    /* If there are tiny polygons, we replace the original polygons by the   */
    /* subset of non-tiny polygons.                                          */
    
    int *n_non_tiny_polycorn = (int*) malloc(n_non_tiny_poly * sizeof(int));
    int *non_tiny_polygon_id = (int*) malloc(n_non_tiny_poly * sizeof(int));
    n_non_tiny_poly = 0;
    for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
      if (!poly_has_tiny_area[poly_indx]) {
	n_non_tiny_polycorn[n_non_tiny_poly] = n_polycorn[poly_indx];
	non_tiny_polygon_id[n_non_tiny_poly] = polygon_id[poly_indx];
	n_non_tiny_poly++;
      }
    }
    POINT **non_tiny_polycorn =
      (POINT**) malloc(n_non_tiny_poly * sizeof(POINT*));
    for (int poly_indx = 0; poly_indx < n_non_tiny_poly; poly_indx++) {
      non_tiny_polycorn[poly_indx] =
	(POINT*) malloc(n_non_tiny_polycorn[poly_indx] * sizeof(POINT));
    }
    n_non_tiny_poly = 0;
    for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
      if (!poly_has_tiny_area[poly_indx]) {
	for (int corn_indx = 0;
	     corn_indx < n_polycorn[poly_indx];
	     corn_indx++) {
	  non_tiny_polycorn[n_non_tiny_poly][corn_indx].x =
	    polycorn[poly_indx][corn_indx].x;
	  non_tiny_polycorn[n_non_tiny_poly][corn_indx].y =
	    polycorn[poly_indx][corn_indx].y;
	}
	n_non_tiny_poly++;
      }
    }
    
    /* Free the memory used by the original polygons. */
    
    free(polygon_id);
    free(n_polycorn);
    for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
      free(polycorn[poly_indx]);
    }
    free(polycorn);
    
    /* Copy the non-tiny polygons to the variables used by the original      */
    /* polygons.                                                             */
    
    n_poly = n_non_tiny_poly;
    polygon_id = (int*) malloc(n_poly * sizeof(int));
    n_polycorn = (int*) malloc(n_poly * sizeof(int));
    for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
      polygon_id[poly_indx] = non_tiny_polygon_id[poly_indx];
      n_polycorn[poly_indx] = n_non_tiny_polycorn[poly_indx];
    }
    polycorn = (POINT**) malloc(n_poly * sizeof(POINT*));
    for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
      polycorn[poly_indx] =
	(POINT*) malloc(n_polycorn[poly_indx] * sizeof(POINT));
    }    
    for (int poly_indx = 0; poly_indx < n_poly; poly_indx++) {
      for (int corn_indx = 0; corn_indx < n_polycorn[poly_indx]; corn_indx++) {
	polycorn[poly_indx][corn_indx].x =
	  non_tiny_polycorn[poly_indx][corn_indx].x;
	polycorn[poly_indx][corn_indx].y =
	  non_tiny_polycorn[poly_indx][corn_indx].y;	
      }
    }
    
    /* Free the memory used by the temporary storage for the non-tiny        */
    /* polygons.                                                             */
    
    free(non_tiny_polygon_id);
    for (int poly_indx = 0; poly_indx < n_non_tiny_poly; poly_indx++) {
      free(non_tiny_polycorn[poly_indx]);
    }
    free(non_tiny_polycorn);
    free(n_non_tiny_polycorn);
  }
  free(poly_has_tiny_area);
  
  return;
}

/*****************************************************************************/
/* Function to make regions from polygons. Region IDs in the .gen file must  */
/* be nonnegative.                                                           */

void make_region (void)
{
  BOOLEAN repeat;
  int i, j, last_id, min_id;
  
   /* Compute which polygons are holes */
  poly_is_hole = calloc(n_poly, sizeof(int));
  for(j=0; j<n_poly; j++) {
    if(polygon_area(n_polycorn[j], polycorn[j]) < 0 ) {
      poly_is_hole[j] = TRUE;
    } else {
      poly_is_hole[j] = FALSE;
    }
  }
  
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
  region_na = (int*) calloc(n_reg, sizeof(int));
  region_perimeter = (double*) calloc(n_reg, sizeof(double));
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

void read_map (char *map_file_name)
{
  FILE *map_file;
  if ((map_file = fopen(map_file_name,"r")) == NULL) {
    fprintf(stderr,"ERROR: Cannot find map file.\n");
    exit(1);
  }
  size_t namelen = strlen(map_file_name), extlenjson = strlen(".json"), extlengeojson = strlen(".geojson"), extlengen = strlen(".gen");
  if(namelen >= extlengen && !strcmp(map_file_name + namelen - extlengen, ".gen")){
    count_poly(map_file);
    fclose(map_file);
    map_file = fopen(map_file_name,"r");
    count_corn(map_file);
    fclose(map_file);
    map_file = fopen(map_file_name,"r");
    read_corn(map_file);
    fclose(map_file);
  }else if(namelen >= extlenjson && !strcmp(map_file_name + namelen - extlenjson, ".json")){
    read_geojson(map_file);
    fclose(map_file);
  }else if(namelen >= extlengeojson && !strcmp(map_file_name + namelen - extlengeojson, ".geojson")){
    read_geojson(map_file);
    fclose(map_file);
  }else{
    fprintf(stderr,"ERROR: Map file not in proper file format. Map file needs to be a .gen file, .json file or .geojson file.\n");
    exit(1);
  }
  remove_tiny_polygons();
  make_region();
  fprintf(stderr, "%i polygon(s), %i region(s)\n", n_poly, n_reg);
  
  return;
}
