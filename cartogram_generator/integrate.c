/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cartogram.h"

/******************************** Definitions. *******************************/

#define ABS_TOL (MIN(lx, ly) * 1e-6)
#define INC_AFTER_ACC (1.1)
#define DEC_AFTER_NOT_ACC (0.75)

/***************************** Global variables. *****************************/

double *gridvx, *gridvy, *grid_fluxx_init, *grid_fluxy_init;
fftw_plan plan_grid_fluxx_init, plan_grid_fluxy_init;

/**************************** Function prototypes. ***************************/

void init_gridv (void);
void calcv (double t);

/*****************************************************************************/
/* Function to initialize the Fourier transforms of gridvx[] and gridvy[] at */
/* every point on the lx-times-ly grid at t = 0. After this function has     */
/* finished, we do not need to do any further Fourier transforms for this    */
/* round of integration                                                      */

void init_gridv (void)
{
  double di, dlx, dly;
  int i, j;

  dlx = (double) lx;            /* We must typecast. Otherwise the ratios in */
  dly = (double) ly;            /* the denominator will evaluate as zero.    */

  /* There is a bit less typing later on if we divide now by 4*lx*ly because */
  /* then the REDFT01 transform of rho_ft[0] will directly be the mean of    */
  /* rho_init[].                                                             */

  for (i=0; i<lx*ly; i++)
    rho_ft[i] /= 4*lx*ly;

  /* We temporarily insert the Fourier coefficients for the x- and           */
  /* y-components of the flux vector in the arrays grid_fluxx_init[] and     */
  /* grid_fluxy_init[].                                                      */

  for (i=0; i<lx-1; i++) {
    di = (double) i;
    for (j=0; j<ly; j++)
      grid_fluxx_init[i*ly + j] =
      	- rho_ft[(i+1)*ly + j] /
      	(PI * ((di+1)/dlx + (j/(di+1)) * (j/dly) * (dlx/dly)));
  }
  for (j=0; j<ly; j++)
    grid_fluxx_init[(lx-1)*ly + j] = 0.0;
  for (i=0; i<lx; i++) {
    di = (double) i;
    for (j=0; j<ly-1; j++)
      grid_fluxy_init[i*ly + j] =
	      -rho_ft[i*ly + j + 1] /
	      (PI * ((di/(j+1)) * (di/dlx) * (dly/dlx) + (j+1)/dly));
  }
  for (i=0; i<lx; i++)
    grid_fluxy_init[i*ly + ly - 1] = 0.0;

  /* Compute the flux vector and store the result in grid_fluxx_init[] and   */
  /* grid_fluxy_init[].                                                      */

  fftw_execute(plan_grid_fluxx_init);
  fftw_execute(plan_grid_fluxy_init);

  return;
}

/*****************************************************************************/
/* Function to calculate the velocity at the grid points (x, y) with x =     */
/* 0.5, 1.5, ..., lx-0.5 and y = 0.5, 1.5, ..., ly-0.5 at time t.            */

void calcv (double t)
{
  double rho;
  int k;

  #pragma omp parallel for private(rho)
  for (k=0; k<lx*ly; k++) {
    rho = rho_ft[0] + (1.0-t) * (rho_init[k] - rho_ft[0]);
    gridvx[k] = -grid_fluxx_init[k] / rho;
    gridvy[k] = -grid_fluxy_init[k] / rho;
  }

  return;
}

/*****************************************************************************/
/* Function to bilinearly interpolate a numerical array grid[0..lx*ly-1]     */
/* whose entries are numbers for the positions:                              */
/* x = (0.5, 1.5, ..., lx-0.5), y = (0.5, 1.5, ..., ly-0.5).                 */
/* The final argument "zero" can take two possible values: "x" or "y". If    */
/* zero==x, the interpolated function is forced to return 0 if x=0 or x=lx.  */
/* This option is suitable fo interpolating from gridvx because there can be */
/* no flow through the boundary. If zero==y, the interpolation returns 0 if  */
/* y=0 or y=ly, suitable for gridvy. The unconstrained boundary will be      */
/* determined by continuing the function value at 0.5 (or lx-0.5 or ly-0.5)  */
/* all the way to the edge (i.e. the slope is 0 consistent with a cosine     */
/* transform).                                                               */

double interpol (double x, double y, double *grid, char *zero)
{
  double delta_x, delta_y, fx0y0, fx0y1, fx1y0, fx1y1, x0, x1, y0, y1;

  if (x<0 || x>lx || y<0 || y>ly) {
    fprintf(stderr, "ERROR: coordinate outside bounding box in interpol().\n");
    fprintf(stderr, "x=%f, y=%f\n", x, y);
    exit(1);
  }
  if (strcmp(zero, "x")!=0 && strcmp(zero, "y")!=0) {
    fprintf(stderr, "ERROR: unknown argument zero in interpol().\n");
    exit(1);
  }

  x0 =                             /* Nearest grid point smaller than x.     */
    MAX(0.0, floor(x+0.5) - 0.5);  /* Exception: if x<0.5, x0 becomes 0.0.   */
  x1 =                             /* Nearest grid point larger than x.      */
    MIN(lx, floor(x+0.5) + 0.5);   /* Exception: if x>lx-0.5, x1 becomes lx. */
  y0 = MAX(0.0, floor(y+0.5) - 0.5);                     /* Similarly for y. */
  y1 = MIN(ly, floor(y+0.5) + 0.5);
  delta_x = (x-x0) / (x1-x0);   /* On a scale from 0 to 1, how far is x (or */
  delta_y = (y-y0) / (y1-y0);   /* y) away from x0 (or y0)? 1 means x=x1.   */

  /* Function value at (x0, y0). */

  if ((x<0.5 && y<0.5) || (x<0.5 && strcmp(zero, "x")==0) ||
          (y<0.5 && strcmp(zero, "y")==0))
    fx0y0 = 0.0;
  else
    fx0y0 = grid[(int)x0*ly + (int)y0];

  /* Function value at (x0, y1). */

  if ((x<0.5 && y>=ly-0.5) || (x<0.5 && strcmp(zero, "x")==0) ||
          (y>=ly-0.5 && strcmp(zero, "y")==0))
    fx0y1 = 0.0;
  else if (x>=0.5 && y>=ly-0.5 && strcmp(zero, "x")==0)
    fx0y1 = grid[(int)x0*ly + ly -1];
  else
    fx0y1 = grid[(int)x0*ly + (int)y1];

  /* Function value at (x1, y0). */

  if ((x>=lx-0.5 && y<0.5) || (x>=lx-0.5 && strcmp(zero, "x")==0) ||
          (y<0.5 && strcmp(zero, "y")==0))
    fx1y0 = 0.0;
  else if (x>=lx-0.5 && y>=0.5 && strcmp(zero, "y")==0)
    fx1y0 = grid[(lx-1)*ly + (int)y0];
  else
    fx1y0 = grid[(int)x1*ly + (int)y0];

  /* Function value at (x1, y1). */

  if ((x>=lx-0.5 && y>=ly-0.5) || (x>=lx-0.5 && strcmp(zero, "x")==0) ||
          (y>=ly-0.5 && strcmp(zero, "y")==0))
    fx1y1 = 0.0;
  else if (x>=lx-0.5 && y<ly-0.5 && strcmp(zero, "y")==0)
    fx1y1 = grid[(lx-1)*ly + (int)y1];
  else if (x<lx-0.5 && y>=ly-0.5 && strcmp(zero, "x")==0)
    fx1y1 = grid[(int)x1*ly + ly - 1];
  else
    fx1y1 = grid[(int)x1*ly + (int)y1];

  return (1-delta_x)*(1-delta_y)*fx0y0 + (1-delta_x)*delta_y*fx0y1
          + delta_x*(1-delta_y)*fx1y0 + delta_x*delta_y*fx1y1;
}

/*****************************************************************************/
/*************** Function to integrate the equations of motion. **************/

void integrate (void)
{
  BOOLEAN accept;
  double delta_t, t, *vx_intp, *vx_intp_half, *vy_intp, *vy_intp_half,
          *xeul, *xmid, *yeul, *ymid;
  int iter, k;

  /****************** Allocate memory for the velocity grid. *****************/

  gridvx = (double*) malloc(lx * ly * sizeof(double));
  gridvy = (double*) malloc(lx * ly * sizeof(double));

  /*************** Allocate memory for the Fourier transforms. ***************/

  grid_fluxx_init = (double*) fftw_malloc(lx * ly * sizeof(double));
  grid_fluxy_init = (double*) fftw_malloc(lx * ly * sizeof(double));

  /************ Prepare the fftw plans for the Fourier transforms. ***********/
  plan_grid_fluxx_init =
          fftw_plan_r2r_2d(lx, ly, grid_fluxx_init, grid_fluxx_init,
		      FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  plan_grid_fluxy_init =
          fftw_plan_r2r_2d(lx, ly, grid_fluxy_init, grid_fluxy_init,
		      FFTW_REDFT01, FFTW_RODFT01, FFTW_ESTIMATE);

  /* xeul[i*ly+j] and yeul[i*ly+j] will be the new displacement proposed by  */
  /* a simple Euler step (i.e. moving a full time interval delta_t with the  */
  /* velocity at time t and position (i + 0.5 + xdisp[i*ly+j],               */
  /* j + 0.5 + xdisp[i*ly+j]).                                               */

  xeul = (double*) malloc(lx * ly * sizeof(double));
  yeul = (double*) malloc(lx * ly * sizeof(double));

  /* xmid[i*ly+j] and yeul[i*ly+j] will be the new displacement proposed by  */
  /* the midpoint method (see comment below for the formula).                */

  xmid = (double*) malloc(lx * ly * sizeof(double));
  ymid = (double*) malloc(lx * ly * sizeof(double));

  /* (vx_intp, vy_intp) will be the velocity at position (xproj, yproj) at   */
  /* time t.                                                                 */

  vx_intp = (double*) malloc(lx * ly * sizeof(double));
  vy_intp = (double*) malloc(lx * ly * sizeof(double));

  /* (vx_intp_half, vy_intp_half) will be the velocity at the midpoint       */
  /* (xproj + 0.5*delta_t*vx_intp, yproj + 0.5*delta_t*vy_intp) at time      */
  /* t + 0.5*delta_t.                                                        */

  vx_intp_half = (double*) malloc(lx * ly * sizeof(double));
  vy_intp_half = (double*) malloc(lx * ly * sizeof(double));

  /********* Initialize grids for vx and vy using Fourier transforms. ********/

  init_gridv();
  t = 0.0;
  delta_t = 1e-2;
  iter = 0;

  /******************************** Integrate. *******************************/

  do {
    calcv(t);
    #pragma omp parallel for
    for (k=0; k<lx*ly; k++) {

      /* We know, either because of the initialization or because of the     */
      /* check at the end of the last iteration, that (xproj[k], yproj[k])   */
      /* is inside the rectangle [0, lx] x [0, ly]. This fact guarantees     */
      /* that interpol() is given a point that cannot cause it to fail.      */

      vx_intp[k] = interpol(xproj[k], yproj[k], gridvx, "x");
      vy_intp[k] = interpol(xproj[k], yproj[k], gridvy, "y");
    }
    accept = FALSE;
    while (!accept) {

      /* Simple Euler step. */

      #pragma omp parallel for
      for (k=0; k<lx*ly; k++) {
      	xeul[k] = xproj[k] + vx_intp[k] * delta_t;
      	yeul[k] = yproj[k] + vy_intp[k] * delta_t;
      }

      /* Use "explicit midpoint method".                                     */
      /* x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),                  */
      /*                        y + 0.5*delta_t*v_y(x,y,t),                  */
      /*                        t + 0.5*delta_t)                             */
      /* and similarly for y.                                                */

      calcv(t + 0.5*delta_t);

      /* Make sure we do not pass a point outside [0, lx] x [0, ly] to       */
      /* interpol(). Otherwise we decrease the time step further down and    */
      /* try again.                                                          */

      accept = TRUE;
      for (k=0; k<lx*ly; k++)
      	if (xproj[k] + 0.5*delta_t*vx_intp[k] < 0.0 ||
      	    xproj[k] + 0.5*delta_t*vx_intp[k] > lx ||
      	    yproj[k] + 0.5*delta_t*vy_intp[k] < 0.0 ||
      	    yproj[k] + 0.5*delta_t*vy_intp[k] > ly)
      	  accept = FALSE;
      if (accept) {

      	/* OK, we can run interpol(). */

        #pragma omp parallel for
      	for (k=0; k<lx*ly; k++) {
      	  vx_intp_half[k] = interpol(xproj[k] + 0.5*delta_t*vx_intp[k],
      				    yproj[k] + 0.5*delta_t*vy_intp[k],
      				    gridvx, "x");
      	  vy_intp_half[k] = interpol(xproj[k] + 0.5*delta_t*vx_intp[k],
      				    yproj[k] + 0.5*delta_t*vy_intp[k],
      				    gridvy, "y");
      	  xmid[k] = xproj[k] + vx_intp_half[k] * delta_t;
      	  ymid[k] = yproj[k] + vy_intp_half[k] * delta_t;

      	  /* Do not accept the integration step if the maximum distance      */
      	  /* between the Euler and midpoint proposals exceeds ABS_TOL.       */
      	  /* Neither should we accept the integration step if one of these   */
      	  /* positions wandered out of the boundaries. If it happened,       */
      	  /* decrease the time step.                                         */

      	  if ((xmid[k]-xeul[k]) * (xmid[k]-xeul[k]) +
      	          (ymid[k]-yeul[k]) * (ymid[k]-yeul[k]) > ABS_TOL ||
      	          xmid[k] < 0.0 || xmid[k] > lx || ymid[k] < 0.0 || ymid[k] > ly)
      	    accept = FALSE;
        }
      }
      if (!accept)
	      delta_t *= DEC_AFTER_NOT_ACC;
    }

    /* Control output. */

    if (iter % 10 == 0)
      printf("iter = %i, t = %e, delta_t = %e\n", iter, t, delta_t);

    /* When we get here, the integration step was accepted. */

    t += delta_t;
    iter++;
    for (k=0; k<lx*ly; k++) {
	    xproj[k] = xmid[k];
	    yproj[k] = ymid[k];
    }
    delta_t *= INC_AFTER_ACC;           /* Try a larger step size next time. */

  } while (t < 1.0);

  /* Free memory. */

  fftw_destroy_plan(plan_grid_fluxx_init);
  fftw_destroy_plan(plan_grid_fluxy_init);
  free(gridvx);
  free(gridvy);
  fftw_free(grid_fluxx_init);
  fftw_free(grid_fluxy_init);
  free(xeul);
  free(yeul);
  free(xmid);
  free(ymid);
  free(vx_intp);
  free(vy_intp);
  free(vx_intp_half);
  free(vy_intp_half);

  return;
}
