/* Program for the construction of cartograms. The program needs polygon     */
/* coordinates and a value for the (relative) target area for each region as */
/* input, calculates the new coordinates, and writes these to a file.        */
/* Postscript images of the original map and the cartogram are created.      */

/* If you use output created by this program please acknowledge the use of   */
/* this code and its first publication in:                                   */
/* Michael T. Gastner, Vivien Seguy and Pratyush More, "A fast flow-based    */
/* algorithm for creating density-equalizing map  projections" (submitted    */
/* for publication).                                                         */

/* The code must be called as "./cartogram.exe input.gen target_area.dat".   */
/* Additional arguments are                                                  */
/* - diff: replaces the fast flow-based algorithm by the Gastner-Newman      */
/*         diffusion method,                                                 */
/* - inv: calculate the inverse transformation and prepare a map             */
/*        (invtransf.eps) that shows the original polygon coordinates        */
/*        together with a graticule where each deformed quadrilateral        */
/*        contains equal population.                                         */
/* The input coordinates in the first command-line argument must be in       */
/* ArcInfo "generate" format of the type:                                    */

/* 2302 Maine02              ID for region followed by optional description. */
/* 0.302204 -0.188090        Pairs of x-, y-coordinates. Orientation along   */
/* 0.302716 -0.187835        outer boundaries must be clockwise. If a        */
/* ...                       polygon has a hole, the inner boundary must be  */
/* 0.303897 -0.193159        anticlockwise.                                  */
/* 0.302204 -0.188090        Regions are not permitted to overlap. This is   */
/* END                       *not* checked by this code!                     */
/* 2301 Maine01                                                              */
/* 0.333358 -0.200693                                                        */
/* ...                                                                       */
/* 0.333358 -0.200693                                                        */
/* END                                                                       */
/* 2301 Maine01              IDs can repeat if a region consists of          */
/* 0.334699 -0.204771        multiple polygons.                              */
/* ...                                                                       */
/* 0.334699 -0.204771                                                        */
/* END                       Each polygon terminates with END.               */
/* END                       One more END signals the end of the file.       */

/* The target area in the second command-line argument must be given in the  */
/* space-delimited format                                                    */
/* region_ID target_area (optional comment)                                  */
/* For example,                                                              */

/* 1 9.0 Alabama                                                             */
/* 4 11.0 Arizona                                                            */
/* 5 6.0 Arkansas                                                            */
/* ...                                                                       */

/* Output:                                                                   */
/* (1) The output coordinates are written to cartogram.gen in ArcInfo        */
/*     "generate" format.                                                    */
/* (2) A postscript image of the original input is prepared as map.eps, an   */
/*     image of the cartogram as cartogram.eps.                              */
/* (3) The relative area errors are printed to area_error.dat.               */

/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>
#include "cartogram.h"

/************* Error message for unknown command-line arguments. *************/

void addl_comm_arg_err (void)
{
  fprintf(stderr, "Additional optional arguments that can be passed:\n");
  fprintf(stderr, "diff: use Gastner-Newman (i.e. diffusion) ");
  fprintf(stderr, "instead of fast flow-based method\n");
  fprintf(stderr, "inv:  print inverse transform to file\n");
  exit(1);
}

/*********************************** Main. ***********************************/

int main (int argc, char* argv[])
{
  BOOLEAN diff, inv;
  double mae;
  int i, integration, j;
  
  /**************************** Parse the command. ***************************/

  if (argc < 3) {
    fprintf(stderr, "ERROR: Use as ./cartogram input.gen target_area.dat\n");
    addl_comm_arg_err();
  }
  diff = FALSE;
  inv = FALSE;
  for (i=3; i<argc; i++) {
    if (!diff && strcmp(argv[i], "diff") == 0) {
      diff = TRUE;
      printf("Using Gastner-Newman (i.e. diffusion) method\n");
    }
    else if (!inv && strcmp(argv[i], "inv") == 0) {
      inv = TRUE;
      printf("Will print inverse transform\n");
    }
    else if (strcmp(argv[i], "diff") != 0 && strcmp(argv[i], "inv") != 0) {
      fprintf(stderr, "\nERROR: Unknown command-line argument %s\n", argv[i]);
      addl_comm_arg_err();
    }
  }
    
  /* The equations of motion are integrated on an lx-times-ly square grid.   */
  /* The parameter L will become max(lx, ly). Optimal values for lx and ly   */
  /* are computed in fill_with_density.c. According to the FFTW              */
  /* documentation, "transforms whose sizes are powers of 2 are especially   */
  /* fast". For reasons of efficiency, we only allow powers of 2.            */
  
  if ((L <= 0) || ((L & (~L + 1)) != L)) {
    fprintf(stderr,"ERROR: L must be an integer power of 2.\n");
    exit(1);
  }
  
  /* If a region contains exactly zero population, it will be replaced by    */
  /* MIN_POP_FAC times the smallest positive population in any region.       */
  /* Obviously, MIN_POP_FAC should be <1. A value around 0.2 is sensible.    */
  /* Regions of zero area cause mathematical problems and are almost never   */
  /* aesthetically desirable.                                                */
  
  if (MIN_POP_FAC >= 1.0 || MIN_POP_FAC <=0.0) {
    fprintf(stderr,"ERROR: MIN_POP_FAC must be < 1.0.\n");
    exit(1);
  }
  
  /* We leave some space between the edges of the lx-times-ly grid and the   */
  /* mapped regions so that the cartogram is unaffected by the particular    */
  /* choice of boundary conditions. The parameter PADDING determines how     */
  /* much space we leave. A value around 1.5 is sensible.                    */
  
  if (PADDING < 1.0) {
    fprintf(stderr,"ERROR: PADDING must be >= 1.0.\n");
    exit(1);
  }
  
  /***************************** Read input data. ****************************/
  
  /* Read the original polygon coordinates, fill the lx-times-ly grid with   */
  /* density and print a map. rho_ft[] will be filled with the Fourier       */
  /* transform of the initial density.                                       */
  
  fill_with_density1(argv[1], argv[2], inv);
  
  /*************** Allocate memory for the projected positions. **************/
  
  proj = (POINT*) malloc(lx * ly * sizeof(POINT));
  cartcorn = (POINT**) malloc(n_poly * sizeof(POINT*));
  for (i=0; i<n_poly; i++)
    cartcorn[i] = (POINT*) malloc(n_polycorn[i] * sizeof(POINT));
  area_err = (double*) malloc(n_reg * sizeof(double));
  cart_area = (double*) malloc(n_reg * sizeof(double));
  
  /* proj[i*ly+j] will store the current position of the point that started  */
  /* at (i+0.5, j+0.5).                                                      */
  
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      proj[i*ly + j].x = i + 0.5;
      proj[i*ly + j].y = j + 0.5;
    }
  
  /* Print a map of the input coordinates in eps-format. The last argument   */
  /* is FALSE because we do not need to show the graticule.                  */
  ps_figure("map.eps", polycorn, proj, FALSE); 
  
  /************** First integration of the equations of motion. **************/
  
  printf("Starting integration 1\n");
  if (!diff)
    ffb_integrate();
  else
    diff_integrate();
  project(FALSE);  /* FALSE because we do not need to project the graticule. */
  mae = max_area_err(area_err, cart_area);
  printf("max. abs. area error: %f\n", mae);
  
  /********* Additional integrations to come closer to target areas. *********/
  
  proj2 = (POINT*) malloc(lx * ly * sizeof(POINT));
  integration = 1;  
  while (mae > MAX_PERMITTED_AREA_ERROR) {
    fill_with_density2();
    
    /* Copy the current graticule before resetting. We will construct the    */
    /* final graticule by interpolating proj2 on the basis of proj.          */
    
    for (i=0; i<lx; i++)
      for (j=0; j<ly; j++) {
	proj2[i*ly + j].x = proj[i*ly + j].x;
	proj2[i*ly + j].y = proj[i*ly + j].y;
      }
    for (i=0; i<lx; i++)
      for (j=0; j<ly; j++) {
      	proj[i*ly + j].x = i + 0.5;
      	proj[i*ly + j].y = j + 0.5;
      }
    integration++;
    printf("Starting integration %d\n", integration);
    if (!diff)
      ffb_integrate();
    else
      diff_integrate();
    project(TRUE);     /* TRUE because we need to project the graticule too. */

    /* Overwrite proj with proj2. */
    
    for (i=0; i<lx; i++)
      for (j=0; j<ly; j++) {
      	proj[i*ly + j].x = proj2[i*ly + j].x;
      	proj[i*ly + j].y = proj2[i*ly + j].y;
      }
    mae = max_area_err(area_err, cart_area);
    printf("max. abs. area error: %f\n", mae);
  }

  /* Print the cartogram in eps-format. Set the final argument to TRUE if    */
  /* you want to add the graticule.                                          */
  
  ps_figure("cartogram.eps", cartcorn, proj, FALSE);

  /* Print additional output files. */
  output_to_ascii();    /* Print coordinates in .gen format and area errors. */
  if (inv)
    inv_project();        /* Show the graticules from the inverse transform. */
  
  /******************************* Free memory. ******************************/
  
  fftw_destroy_plan(plan_fwd);
  fftw_free(rho_ft);
  fftw_free(rho_init);
  for (i=0; i<n_poly; i++)
    free(polycorn[i]);
  free(polycorn);
  for (i=0; i<n_poly; i++)
    free(cartcorn[i]);
  free(cartcorn);
  free(n_polycorn);
  free(polygon_id);
  free(region_id);
  free(region_id_inv);
  for (i=0; i<n_reg; i++)
    free(polyinreg[i]);
  free(polyinreg);
  free(n_polyinreg);
  free(proj);
  free(proj2);
  free(target_area);
  free(area_err);
  free(cart_area);
  if (inv) {
    for (i=0; i<n_poly; i++)
      free(origcorn[i]);
    free(origcorn);
  }
  return 0;
}
