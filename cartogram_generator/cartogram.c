/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include "cartogram.h"

/**************************** Function prototypes. ***************************/

void project (BOOLEAN proj_graticule);
void output_to_ascii (void);

/*****************************************************************************/
/********** Function to project the polygons in the input .gen file. *********/

void project (BOOLEAN proj_graticule)
{
  double *xdisp, x2, *ydisp, y2;
  int i, j;

  /* The displacement vector (xdisp[i*ly+j], ydisp[i*ly+j]) is the point     */
  /* that was initially at (i+0.5, j+0.5). We work with (xdisp, ydisp)       */
  /* instead of (xproj, yproj) so that we can use the function interpol()    */
  /* defined in integrate.c.                                                 */

  xdisp = (double*) malloc(lx * ly * sizeof(double));
  ydisp = (double*) malloc(lx * ly * sizeof(double));
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      xdisp[i*ly + j] = xproj[i*ly + j] - i - 0.5;
      ydisp[i*ly + j] = yproj[i*ly + j] - j - 0.5;
    }

  /********************* Project the polygon coordinates. ********************/
  for (i=0; i<n_poly; i++)
    for (j=0; j<n_polycorn[i]; j++) {
      cartcorn[i][j].x =
            	interpol(polycorn[i][j].x, polycorn[i][j].y, xdisp, "x")
            	+ polycorn[i][j].x;
      cartcorn[i][j].y =
            	interpol(polycorn[i][j].x, polycorn[i][j].y, ydisp, "y")
            	+ polycorn[i][j].y;
    }
  if (proj_graticule) {

    /******** Project (xproj2, yproj2) on the basis of (xproj, yproj). *******/

    for (i=0; i<lx*ly; i++) {
      x2 = xproj2[i];
      y2 = yproj2[i];
      xproj2[i] = interpol(x2, y2, xdisp, "x") + x2;
      yproj2[i] = interpol(x2, y2, ydisp, "y") + y2;
    }
  }

  /******************************* Free memory. ******************************/

  free(xdisp);
  free(ydisp);

  return;
}

/*****************************************************************************/
/** Function to write cartogram polygons and relative area errors to files. **/

void output_to_ascii (void)
{
  double *cart_area, obj_area, sum_cart_area, sum_target_area;
  FILE *err_file = fopen("error.dat", "w"),
          *gen_file = fopen("cartogram.gen", "w");
  int i, j;

  /***************** Output of coordinates to cartogram.gen. *****************/

  for (i=0; i<n_poly; i++) {
    fprintf(gen_file, "%i\n", polygon_id[i]);
    for (j=0; j<n_polycorn[i]; j++)
      fprintf(gen_file, "%f %f\n", cartcorn[i][j].x, cartcorn[i][j].y);
    fprintf(gen_file, "END\n");
  }
  fprintf(gen_file, "END\n");
  fclose(gen_file);

  /*************** Output of relative area errors to error.dat. **************/

  err_file = fopen("error.dat", "w");
  cart_area = (double*) calloc(n_reg, sizeof(double));
  for (i=0; i<n_reg; i++)
    for (j=0; j<n_polyinreg[i]; j++)
      cart_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
              cartcorn[polyinreg[i][j]]);
  for (i=0, sum_target_area=0.0; i<n_reg; i++)
    sum_target_area += target_area[i];
  for (i=0, sum_cart_area=0.0; i<n_reg; i++)
    sum_cart_area += cart_area[i];
  for (i=0; i<n_reg; i++) {
    fprintf(err_file, "region %i: ", region_id[i]);
    obj_area = target_area[i] * sum_cart_area / sum_target_area;
    fprintf(err_file,
	          "objective area = %f, actual area = %f, relative error = %f\n",
	          obj_area, cart_area[i], (cart_area[i]-obj_area)/obj_area);
  }
  fclose(err_file);

  /******************************* Free memory. ******************************/

  free(cart_area);

  return;
}
