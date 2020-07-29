/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <string.h>
#include "cartogram.h"

/***************************** Global variables. *****************************/

int **xyhalfshift2reg;

/**************************** Function prototypes. ***************************/

void rescale_map (BOOLEAN inv);
void set_inside_value_at_y (int region, POINT pk, POINT pn, int l,
          double poly_minx, int **inside);
void set_inside_values_between_points (int region, POINT pk, POINT pn,
               double poly_minx, int **inside);
void interior (void);
void gaussian_blur (double tot_init_area, double avg_dens);

/*****************************************************************************/
/* Function to change coordinates from (minx, miny, maxx, maxy) to           */
/* (0, 0, LX, LY).                                                           */

void rescale_map (BOOLEAN inv)
{
  double latt_const, new_maxx, new_maxy, new_minx, new_miny;
  int i, j;
  
  /* Minimum dimensions that leave enough space between map and rectangular  */
  /* boundaries.                                                             */

  new_maxx = 0.5 * ((1.0+PADDING)*map_maxx + (1.0-PADDING)*map_minx);
  new_minx = 0.5 * ((1.0-PADDING)*map_maxx + (1.0+PADDING)*map_minx);
  new_maxy = 0.5 * ((1.0+PADDING)*map_maxy + (1.0-PADDING)*map_miny);
  new_miny = 0.5 * ((1.0-PADDING)*map_maxy + (1.0+PADDING)*map_miny);  
  if (map_maxx-map_minx > map_maxy-map_miny) {    
    lx = L;
    latt_const = (new_maxx-new_minx) / L;
    ly = 1 << ((int)ceil(log2((new_maxy-new_miny)/latt_const)));
    new_maxy = 0.5*(map_maxy+map_miny) + 0.5*ly*latt_const;
    new_miny = 0.5*(map_maxy+map_miny) - 0.5*ly*latt_const;
  }
  else {
    ly = L;
    latt_const = (new_maxy-new_miny) / L;
    lx = 1 << ((int) ceil(log2((new_maxx-new_minx) / latt_const)));
    new_maxx = 0.5*(map_maxx+map_minx) + 0.5*lx*latt_const;
    new_minx = 0.5*(map_maxx+map_minx) - 0.5*lx*latt_const;
  }
  fprintf(stderr,
    "Using a %d x %d lattice with bounding box\n\t(%f %f %f %f).\n",
    lx, ly, new_minx, new_miny, new_maxx, new_maxy);
  
  /********************* Rescale all polygon coordinates. ********************/

  for (i=0; i<n_poly; i++)
    for (j=0; j<n_polycorn[i]; j++) {
      polycorn[i][j].x = (polycorn[i][j].x - new_minx) / latt_const;
      polycorn[i][j].y = (polycorn[i][j].y - new_miny) / latt_const;
    }

  /*** If we wish to plot the inverse transform, we must save the original ***/
  /*** polygon coordinates.                                                ***/
  if (inv) {
    origcorn = (POINT**) malloc(n_poly * sizeof(POINT*));
    for (i=0; i<n_poly; i++)
      origcorn[i] = (POINT*) malloc(n_polycorn[i] * sizeof(POINT));
    for (i=0; i<n_poly; i++)
      for (j=0; j<n_polycorn[i]; j++) {
  origcorn[i][j].x = polycorn[i][j].x;
  origcorn[i][j].y = polycorn[i][j].y;
      }
  }
  
  return;
}

/*****************************************************************************/
/* Function to set values of inside[][], used in                             */
/* set_inside_values_for_polygon() below. It sets the value in inside[][]    */
/* for all x-values between poly_minx and the x-value (of the point on the   */
/* line connecting the given two coordinates) that corresponds to the        */
/* current y-value l.                                                        */

void set_inside_value_at_y (int region, POINT pk, POINT pn, int l,
          double poly_minx, int **inside)
{
  double intersection;
  int m;

  /* x-value of the intersection between y = l and the line formed by the    */
  /* coordinates (pkx, pky) and (pnx, pny).                                  */

  intersection = (pn.x-0.5 - (pk.x-0.5)) * (l - (pk.y-0.5)) /
          (pn.y-0.5 - (pk.y-0.5)) + (pk.x-0.5);  
  for (m = (int) poly_minx; m < intersection; m++)
    inside[m][l] = region - inside[m][l] - 1;

  return;
}

/*****************************************************************************/
/* Function that takes two polygon coordinates and loops over the y-values   */
/* between the two input y-coordinates. It updates the value of inside[][]   */
/* for all points between polyminx (for this polygon) and the x-value at all */
/* coordinates on the horizontal line to the left of the line segment        */
/* connecting the input coordinates.                                         */

void set_inside_values_between_points (int region, POINT pk, POINT pn,
               double poly_minx, int **inside)
{
  int l;
  
  /* Loop over all integer y-values between the two input y-coordinates.     */

  for (l = ceil(MIN(pn.y, pk.y) - 0.5); l < MAX(pn.y - 0.5, pk.y - 0.5); l++)
    set_inside_value_at_y(region, pk, pn, l, poly_minx, inside);

  return;
}

/*****************************************************************************/
/**** Function to set values in inside[][] for a particular polygon in a  ****/
/**** region.                                                             ****/

void set_inside_values_for_polygon (int region, int n_polycorn,
            POINT *polycorn, int **inside)
{
  double poly_minx = polycorn[0].x;
  int k, n;

  /************ Determine the minimum x-coordinate of the polygon. ***********/
  
  for (k=0; k<n_polycorn; k++)
    poly_minx = MIN(poly_minx, polycorn[k].x);

  /* Loop over all pairs of consecutive coordinates of polygon.              */

  for (k=0, n=n_polycorn-1; k<n_polycorn; n=k++)
    set_inside_values_between_points(region, polycorn[k], polycorn[n],
             poly_minx, inside);
  return;
}

/*****************************************************************************/
/* Function to determine if a grid point at (x+0.5, y+0.5) is inside one of  */
/* the polygons and, if yes, in which region. The result is stored in the    */
/* array xyhalfshift2reg[x][y]. If (x+0.5, y+0.5) is outside all polygons,   */
/* then xyhalfshift2reg[x][y] = -1. If it is inside region i, then           */
/* xyhalfshift2reg[x][y] = i.                                                */

void interior (void)
{
  int i, j, poly;

  /********************* Initialize xyhalfshift2reg[][]. *********************/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++)
      xyhalfshift2reg[i][j] = -1;

  /**** Use the concept of the "crossing number" to determine the correct ****/
  /**** values of xyhalfshift2reg[][].                                    ****/

  for (i=0; i<n_reg; i++)
    for (j=0; j<n_polyinreg[i]; j++) {
      poly = polyinreg[i][j];
      set_inside_values_for_polygon(i, n_polycorn[poly], polycorn[poly],
            xyhalfshift2reg);
    }
  return;
}

/*****************************************************************************/
/************ Function to smoothen the density near polygon edges. ***********/

void gaussian_blur (double tot_init_area, double avg_dens)
{
  double prefactor, scale_i, scale_j;
  fftw_plan plan_bwd;
  int i, j;

  /* Prepare the backward transform (i.e. from Fourier to real space). */

  plan_bwd = fftw_plan_r2r_2d(lx, ly, rho_ft, rho_init,
          FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);

  /* Upon forward and backward transform, rho_init is multiplied by the      */
  /* factor 4*lx*ly.                                                         */

  for (i=0; i<lx*ly; i++)
    rho_init[i] /= 4*lx*ly;

  /* At the start of this function, rho_ft is the unsmoothed density. We     */
  /* first have to transform it to Fourier space.                            */

  fftw_execute(plan_fwd);

  /* Now perform the Gaussian blur.                                          */

  prefactor = -0.5 * BLUR_WIDTH * BLUR_WIDTH * PI * PI;
  for (i=0; i<lx; i++) {
    scale_i = (double) i / lx;
    for (j=0; j<ly; j++) {
      scale_j = (double) j /ly;
      rho_ft[i*ly+j] *=
  exp(prefactor * (scale_i*scale_i + scale_j*scale_j));
    }
  }
  fftw_execute(plan_bwd);      /* Transform back from Fourier to real space. */

  /* Free memory. */

  fftw_destroy_plan(plan_bwd);

  return;
}

/*****************************************************************************/
/* Function to fill the lx-times-ly-grid with the initial density. It reads  */
/* the input target areas and produces a .eps image of the input polygons.   */
/* It also performs a Gaussian blur on the input density. This function      */
/* should only be called once, namely before the first round of integration. */
/* Afterwards use fill_with_density2() below.                                */
/* Input: - map_file_name: name of the .gen-file/.json-file containing the polygon      */
/*                         coordinates                                       */
/*        - area_file_name: name of the .dat-file containing the target      */
/*                          areas                                            */
/*        - inv: do we need the inverse transform? If yes, we must store the */
/*               original polygon coordinates.                               */
/* Return value: is there only one region? If yes, we neither need to        */
/*               rasterize the map nor read the target areas because the     */
/*               output cartogram is simply a linear transformation of the   */
/*               input map.                                                  */

BOOLEAN fill_with_density1 (char *map_file_name, char *area_file_name,
          BOOLEAN inv, BOOLEAN eps)
{
  char line[MAX_STRING_LENGTH];
  double area, avg_dens, *dens, *init_area, tot_init_area,
          tot_target_area, r, g, b;;

  FILE *area_file;
  int i, id, j;

  /************************** Read the coordinates. **************************/

  read_map(map_file_name);

  /******* Allocate memory for the area and area error and the color *********/
  
  cart_area = (double*) malloc(n_reg * sizeof(double));
  area_err = (double*) malloc(n_reg * sizeof(double));
  target_area = (double*) malloc(n_reg * sizeof(double));
  color = (rgb_color*) malloc(n_reg * sizeof(rgb_color));

  /***************** Fit the map on an (lx)*(ly)-square grid. ****************/

  rescale_map(inv);

  /* Print a map of the input coordinates in eps-format. The last argument   */
  /* is FALSE because we do not need to show the graticule.                  */
  
  if (eps)
    ps_figure("map.eps", polycorn, proj, FALSE);
  if (n_reg == 1) {

    /* Placeholder for target area so that max_area_error() produces a       */
    /* finite numeric result.                                                */

    target_area[0] = 1.0;
    return TRUE;  /* It is true that there is only one region. */
  }
  
  /***************************** Allocate memory. ****************************/

  xyhalfshift2reg = (int**) malloc(lx * sizeof(int*));
  for (i=0; i<lx; i++)
    xyhalfshift2reg[i] = (int*) malloc(ly * sizeof(int));
  dens = (double*) malloc(n_reg * sizeof(double));
  init_area = (double*) calloc(n_reg, sizeof(double));
  
  /******* Determine inside which regions the grid points are located. *******/

  interior();

  /************ Read target areas and color information from file. ***********/

  for (i=0; i<n_reg; i++) target_area[i] = -1.0;
  if (area_file_name == NULL){
    area_file = stdin;
  }else if ((area_file = fopen(area_file_name, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open area-file.\n");
    exit(1);
  }
  while (fgets(line, MAX_STRING_LENGTH, area_file) != NULL) {
    id = -1;
    area = -1.0;
    r = 0.96;  
    g = 0.92;
    b = 0.70;
    sscanf(line, "%d%*[, ]%lf%*[, a-zA-Z]%lf%*[, ]%lf%*[, ]%lf", &id, &area, &r, &g, &b);
    if(id != -1 && area != -1.0){
      if (id>max_id || region_id_inv[id]<0) {
        fprintf(stderr, "ERROR: Identifier %d in area-file does not match\n",
         id);
        fprintf(stderr, "       any identifier in map-file.\n");
        exit(1);
      }
      target_area[region_id_inv[id]] = area;
      color[region_id_inv[id]].r = r/255.0;
      color[region_id_inv[id]].g = g/255.0;
      color[region_id_inv[id]].b = b/255.0;
    }else if(id != -1 && area == -1.0){
      char tmpchar1 = '\0', tmpchar2 = '\0';
      sscanf(line, "%*d%*[, ]%c%c", &tmpchar1, &tmpchar2);
      if(tmpchar1 == 'N' && tmpchar2 == 'A'){
        target_area[region_id_inv[id]] = -2.0;
        region_na[region_id_inv[id]] = 1;
      }
    }
  }

  if (area_file != stdin){
    fclose(area_file);
  }

  for (i=0; i<n_reg; i++){
    if (target_area[i] < 0.0 && target_area[i] != -2.0) {
      fprintf(stderr, "ERROR: No target area for region %d.\n", region_id[i]);
      exit(1);
    }
  }

  if (eps){
    ps_figure("map.eps", polycorn, proj, FALSE);
  }

  /****** Replace target areas equal to zero by a small positive value. ******/
  int na_ctr = 0;
  double tmp_tot_target_area = 0.0;
  tot_init_area = 0.0;
  for (i=0; i<n_reg; i++){
    if (region_na[i] == 1){
      na_ctr++;
    }else{
      tmp_tot_target_area += target_area[i];
    }
    for (j=0; j<n_polyinreg[i]; j++){
      init_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
           polycorn[polyinreg[i][j]]);
    }
    tot_init_area += init_area[i];
  }
  
  /****** Calculate region perimeter *********/
  for(i=0; i<n_reg; i++){
    for (j=0; j<n_polyinreg[i]; j++){
      region_perimeter[i] += polygon_perimeter(n_polycorn[polyinreg[i][j]],
       polycorn[polyinreg[i][j]]);
    }
  }
  int first_region = 1;
  double total_NA_ratio = 0;

  for(i=0; i<n_reg; i++){
    if(region_na[i] == 1){
      total_NA_ratio += init_area[i] / tot_init_area;
    }
  }
  
  double total_NA_area = (total_NA_ratio * tmp_tot_target_area) / (1 - total_NA_ratio);
  tmp_tot_target_area += total_NA_area;

  for (i=0; i<n_reg; i++){
    /****** Set target area for regions with NA values *******/
    if(region_na[i] == 1){
      if(first_region == 1){
        fprintf(stderr, "\nSetting area for NA regions:\n");
        first_region = 0;
      }
      target_area[i] = (init_area[i] / tot_init_area) / total_NA_ratio * total_NA_area;
      fprintf(stderr, "%d: %lf\n", region_id[i], target_area[i]);
    }
  }

  fprintf(stderr, "\n");

  /****** Increase target area for regions which will be too small in order to speed ******/
  /****** up cartogram generation process. This happens when -n flag is not set      ******/
  if(use_perimeter_threshold == TRUE){
    fprintf(stderr, "Note: Enlarging extremely small regions using scaled");
    fprintf(stderr, " perimeter threshold. Areas for these regions will be");
    fprintf(stderr, " scaled up. To disable this, please add the -n flag.\n\n");
    int *region_small = (int*) calloc(n_reg, sizeof(int));
    int region_small_ctr = 0;
    double *region_threshold, *region_threshold_area, tot_region_small_area = 0,
      total_perimeter = 0, total_threshold = 0;
    region_threshold = (double*) calloc(n_reg, sizeof(double));
    region_threshold_area = (double*) calloc(n_reg, sizeof(double));
    for(i=0; i<n_reg; i++){
      total_perimeter += region_perimeter[i];
    }
    for(i=0; i<n_reg; i++){
      region_threshold[i] = MAX((region_perimeter[i]/total_perimeter) * MIN_PERIMETER_FAC, 0.00025);
      if((target_area[i]/tmp_tot_target_area < region_threshold[i])){
        region_small[i] = 1;
        region_small_ctr++;
        tot_region_small_area += target_area[i];
      }
    }
    for (i=0; i<n_reg; i++){
      if(region_small[i] == 1){
        total_threshold += region_threshold[i];
      }
    }
    double total_threshold_area = (total_threshold * (tmp_tot_target_area - tot_region_small_area)) / (1 - total_threshold);

    if(region_small_ctr > 0){
      fprintf(stderr, "Enlarging small regions:\n");
    }

    for(i=0; i<n_reg; i++){
      if(region_small[i] == 1){
        region_threshold_area[i] = (region_threshold[i]/total_threshold) * total_threshold_area;
        double old_target_area = target_area[i];
        target_area[i] = region_threshold_area[i];
        tmp_tot_target_area += target_area[i];
        tmp_tot_target_area -= old_target_area;
        fprintf(stderr, "%d: %lf\n", region_id[i], target_area[i]);
      }
    }
    if(region_small_ctr > 0){
      fprintf(stderr, "\n");
    }else{
      fprintf(stderr, "No regions below minimum threshold.\n\n");
    }
    free(region_small);
  }else{
    /* If -n flag is set, regions with zero area will be replaced by MIN_POP_FAC * min_area */
    fprintf(stderr, "Note: Not using scaled perimeter threshold.\n\n");
    double min_area = target_area[0];
    for (i=1; i<n_reg; i++){
      if (target_area[i] > 0.0){
        min_area = MIN(min_area, target_area[i]);
      }
    }
    for (i=0; i<n_reg; i++){
      if (target_area[i] == 0.0){
        target_area[i] = MIN_POP_FAC * min_area;
      }
    }
  }
  
  /**************** Calculate all densities ***************/

  for (i=0; i<n_reg; i++){
    dens[i] = target_area[i] / init_area[i];
  }

  for (i=0, tot_target_area=0.0; i<n_reg; i++){
    tot_target_area += target_area[i];
  }
  avg_dens = tot_target_area / tot_init_area;


  /***** Allocate memory for the Fourier transform of the input density. *****/

  rho_ft = (double*) fftw_malloc(lx * ly * sizeof(double));
  rho_init = (double*) fftw_malloc(lx * ly * sizeof(double));

  /************************** Digitize the density. **************************/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      if (xyhalfshift2reg[i][j]==-1)
  rho_init[i*ly+j] = avg_dens;
      else
  rho_init[i*ly+j] = dens[xyhalfshift2reg[i][j]];
    }

  /*** Plan the Fourier transform already now so that it can be shared with **/
  /*** gaussian_blur().                                                     **/

  plan_fwd = fftw_plan_r2r_2d(lx, ly, rho_init, rho_ft,
          FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);

  /** Smoothen the density profile to avoid uncontrolled distortions around **/
  /** the edges of the polygons.                                            **/

  gaussian_blur(tot_init_area, avg_dens);

  /* Compute rho_ft[], the two dimensional cosine forward Fourier transform  */
  /* of rho_init[]. After the transform we have                              */
  /* rho_ft[i*ly+j] = 4 \sum_{m=0}^{lx-1} \sum_{n=0}^{ly-1}                  */
  /*           rho_init[m][n] cos(pi i (m+1/2) / lx) cos(pi j (n+1/2) / ly). */
  /* We must bear in mind that rho_0[m][n] is the density associated with    */
  /* the coordinates (m+1/2, n+1/2).                                         */

  fftw_execute(plan_fwd);

  /******************************* Free memory. ******************************/

  for (i=0; i<lx; i++)
    free(xyhalfshift2reg[i]);
  free(xyhalfshift2reg);
  free(dens);
  free(init_area);

  return FALSE;        /* The statement "there is only one region" is false. */
}

/*****************************************************************************/
/* Function to fill the lx-times-ly grid with density *after* the first      */
/* round of integration. The main differences compared to                    */
/* fill_with_density1() are that fill_with_density2()                        */
/* - does not assign the target areas,                                       */
/* - does not produce a map,                                                 */
/* - does not perform a Gaussian blur.                                       */

void fill_with_density2 (void)
{
  double avg_dens, *dens, *tmp_area, tot_target_area, tot_tmp_area;
  int i, j;

  /* Copy cartcorn[][] to polycorn[][]. */

  for (i=0; i<n_poly; i++)
    for (j=0; j<n_polycorn[i]; j++)
      polycorn[i][j] = cartcorn[i][j];

  /***************************** Allocate memory. ****************************/

  xyhalfshift2reg = (int**) malloc(lx * sizeof(int*));
  for (i=0; i<lx; i++)
    xyhalfshift2reg[i] = (int*) malloc(ly * sizeof(int));
  dens = (double*) malloc(n_reg * sizeof(double));
  tmp_area = (double*) calloc(n_reg, sizeof(double));

  /******* Determine inside which regions the grid points are located. *******/

  interior();

  /** Calculate all region areas and densities up to this point in the      **/
  /** algorithm.                                                            **/

  for (i=0; i<n_reg; i++)
    for (j=0; j<n_polyinreg[i]; j++)
      tmp_area[i] +=
  polygon_area(n_polycorn[polyinreg[i][j]], polycorn[polyinreg[i][j]]);
  for (i=0; i<n_reg; i++) dens[i] = target_area[i] / tmp_area[i];

  /******************** Calculate the "average density" = ********************/
  /**************** (total target area)/(total initial area). ****************/

  for (i=0, tot_tmp_area=0.0; i<n_reg; i++)
    tot_tmp_area += tmp_area[i];
  for (i=0, tot_target_area=0.0; i<n_reg; i++)
    tot_target_area += target_area[i];
  avg_dens = tot_target_area / tot_tmp_area;

  /************************** Digitize the density. **************************/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      if (xyhalfshift2reg[i][j]==-1)
  rho_init[i*ly+j] = avg_dens;
      else
  rho_init[i*ly+j] = dens[xyhalfshift2reg[i][j]];
    }
  fftw_execute(plan_fwd);

  /******************************* Free memory. ******************************/

  for (i=0; i<lx; i++)
    free(xyhalfshift2reg[i]);
  free(xyhalfshift2reg);
  free(dens);
  free(tmp_area);

  return;
}
