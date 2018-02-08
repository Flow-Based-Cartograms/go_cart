/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cartogram.h"

/******************************** Definitions. *******************************/

#define ABS_TOL (MIN(lx, ly) * 1e-6)
#define INC_AFTER_ACC (1.1)
#define DEC_AFTER_NOT_ACC (0.75)
#define CONV_MAX_CHANGE (MIN(lx, ly) * 1e-9)       /* Convergence criterion. */
#define MAX_ITER (10000)
#define MIN_T (1e3)
#define MAX_T (1e12)

/***************************** Global variables. *****************************/

double *gridvx, *gridvy, *rho;
fftw_plan plan_gridvx, plan_gridvy, plan_rho;

/**************************** Function prototypes. ***************************/
void diff_calcv (double t);

/*****************************************************************************/
/* Function to calculate the velocity at the grid points (x, y) with x =     */
/* 0.5, 1.5, ..., lx-0.5 and y = 0.5, 1.5, ..., ly-0.5 at time t.            */

void diff_calcv (double t)
{
  double dlx, dly;
  int i, j;

  dlx = (double) lx;                        /* We must typecast to prevent   */
  dly = (double) ly;                        /* integer division by lx or ly. */
  
  /* Fill rho with Fourier coefficients.                                     */
  
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++)
      rho[i*ly+j] =
	exp((- (i/dlx)*(i/dlx) - (j/dly)*(j/dly)) * t) * rho_ft[i*ly+j];

  /* Replace rho by cosine Fourier backtransform in both variables.          */
  /* rho[i*ly+j] is the density at position (i+1/2, j+1/2).                  */
  
  fftw_execute(plan_rho);

  /* We temporarily insert the Fourier coefficients for the x- and           */
  /* y-components of the flux vector in the arrays gridvx and gridvy.        */
  
  for (i=0; i<lx-1; i++)
    for (j=0; j<ly; j++)
      gridvx[i*ly + j] =
	rho_ft[(i+1)*ly+j] * (i+1) *
	exp((- ((i+1)/dlx)*((i+1)/dlx) - (j/dly)*(j/dly)) * t) / (PI * dlx);
  for (j=0; j<ly; j++)
    gridvx[(lx-1)*ly + j] = 0.0;
  for (i=0; i<lx; i++) {
    for (j=0; j<ly-1; j++)
      gridvy[i*ly + j] =
	rho_ft[i*ly+j+1] * (j+1) *
	exp((- (i/dlx)*(i/dlx) - ((j+1)/dly)*((j+1)/dly)) * t) / (PI * dly);
  }
  for (i=0; i<lx; i++)
    gridvy[i*ly + ly - 1] = 0.0;

  /* Compute the flux vector and temporarily store the result in gridvx and  */
  /* gridvy.                                                                 */
  
  fftw_execute(plan_gridvx);
  fftw_execute(plan_gridvy);

  /* The velocity is the flux divided by the density.                        */

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      if (rho[i*ly + j] <= 0.0) {
	fprintf(stderr, "ERROR: division by zero in diff_calcv()\n");
	fprintf(stderr, "rho[%d, %d] = %e\n", i, j, rho[i*ly + j]);
	//exit(1);
      }
      gridvx[i*ly + j] /= rho[i*ly + j];      
      gridvy[i*ly + j] /= rho[i*ly + j];
    }
  
  return;
}

/*****************************************************************************/
/* Function to integrate the equations of motion with the diffusion method.  */

void diff_integrate (void)
{
  BOOLEAN accept;
  double delta_t, max_change, t, *vx_intp, *vx_intp_half, *vy_intp,
    *vy_intp_half;
  int iter, k;
  POINT *eul, *mid;
  
  /*************** Allocate memory for the Fourier transforms. ***************/
  rho = fftw_malloc(lx * ly * sizeof(double));
  gridvx = fftw_malloc(lx * ly * sizeof(double));
  gridvy = fftw_malloc(lx * ly * sizeof(double));

  /************ Prepare the fftw plans for the Fourier transforms. ***********/
  
  plan_rho = fftw_plan_r2r_2d(lx, ly, rho, rho,
			      FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);
  plan_gridvx = fftw_plan_r2r_2d(lx, ly, gridvx, gridvx,
				 FFTW_RODFT01, FFTW_REDFT01, FFTW_MEASURE);
  plan_gridvy = fftw_plan_r2r_2d(lx, ly, gridvy, gridvy,
				 FFTW_REDFT01, FFTW_RODFT01, FFTW_MEASURE);

  /* eul[i*ly+j] will be the new position of proj[i*ly+j] proposed by a      */
  /* simple Euler step: move a full time interval delta_t with the velocity  */
  /* at time t and position (proj[i*ly+j].x, proj[i*ly+j].y).                */
  
  eul = (POINT*) malloc(lx * ly * sizeof(POINT));

  /* mid[i*ly+j] will be the new displacement proposed by the midpoint       */
  /* method (see comment below for the formula).                             */
  
  mid = (POINT*) malloc(lx * ly * sizeof(POINT));

  /* (vx_intp, vy_intp) will be the velocity at position (proj.x, proj.y) at */
  /* time t.                                                                 */
  
  vx_intp = (double*) malloc(lx * ly * sizeof(double));
  vy_intp = (double*) malloc(lx * ly * sizeof(double));

  /* (vx_intp_half, vy_intp_half) will be the velocity at the midpoint       */
  /* (proj.x + 0.5*delta_t*vx_intp, proj.y + 0.5*delta_t*vy_intp) at time    */
  /* t + 0.5*delta_t.                                                        */

  vx_intp_half = (double*) malloc(lx * ly * sizeof(double));
  vy_intp_half = (double*) malloc(lx * ly * sizeof(double));
  
  t = 0.0;
  delta_t = 1e-2;                                      /* Initial time step. */
  iter = 0;  
  do {    
    diff_calcv(t);
    for (k=0; k<lx*ly; k++) {
      
      /* We know, either because of the initialization or because of the     */
      /* check at the end of the last iteration, that (proj.x[k], proj.y[k]) */
      /* is inside the rectangle [0, lx] x [0, ly]. This fact guarantees     */
      /* that interpol() is given a point that cannot cause it to fail.      */
      
      vx_intp[k] = interpol(proj[k].x, proj[k].y, gridvx, "x");
      vy_intp[k] = interpol(proj[k].x, proj[k].y, gridvy, "y");
    }
    
    accept = FALSE;
    while (!accept) {
      
      /* Simple Euler step. */
      
      for (k=0; k<lx*ly; k++) {
	eul[k].x = proj[k].x + vx_intp[k] * delta_t;
	eul[k].y = proj[k].y + vy_intp[k] * delta_t;
      }
      
      /* Use "explicit midpoint method".                                     */
      /* x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),                  */
      /*                        y + 0.5*delta_t*v_y(x,y,t),                  */
      /*                        t + 0.5*delta_t)                             */
      /* and similarly for y.                                                */
      
      diff_calcv(t + 0.5*delta_t);
      
      /* Make sure we do not pass a point outside [0, lx] x [0, ly] to       */
      /* interpol(). Otherwise decrease the time step below and try again.   */
      
      accept = TRUE;
      for (k=0; k<lx*ly; k++)
	if (proj[k].x + 0.5*delta_t*vx_intp[k] < 0.0 ||
	    proj[k].x + 0.5*delta_t*vx_intp[k] > lx ||
	    proj[k].y + 0.5*delta_t*vy_intp[k] < 0.0 ||
	    proj[k].y + 0.5*delta_t*vy_intp[k] > ly) {
	  accept = FALSE;	  
	  delta_t *= DEC_AFTER_NOT_ACC;
	  break;
	}      
      if (accept) {                            /* OK, we can run interpol(). */
	for (k=0; k<lx*ly; k++) {
	  vx_intp_half[k] = interpol(proj[k].x + 0.5*delta_t*vx_intp[k],
				     proj[k].y + 0.5*delta_t*vy_intp[k],
				     gridvx, "x");
	  vy_intp_half[k] = interpol(proj[k].x + 0.5*delta_t*vx_intp[k],
				     proj[k].y + 0.5*delta_t*vy_intp[k],
				     gridvy, "y");
	  mid[k].x = proj[k].x + vx_intp_half[k] * delta_t;
	  mid[k].y = proj[k].y + vy_intp_half[k] * delta_t;
	  
	  /* Do not accept the integration step if the maximum squared       */
	  /* difference between the Euler and midpoint proposals exceeds     */
	  /* ABS_TOL. Neither should we accept the integration step if one   */
	  /* of the positions wandered out of the boundaries. If it          */
	  /* happened, decrease the time step.                               */
	  
	  if ((mid[k].x-eul[k].x) * (mid[k].x-eul[k].x) +
	      (mid[k].y-eul[k].y) * (mid[k].y-eul[k].y) > ABS_TOL ||
	      mid[k].x < 0.0 || mid[k].x > lx ||
	      mid[k].y < 0.0 || mid[k].y > ly) {
	    accept = FALSE;
	    delta_t *= DEC_AFTER_NOT_ACC;
	    break;
	  }
	}
      }
    }
    
    /* What is the maximum change in squared displacements between this and  */
    /* the previous integration step?                                        */
    
    for (k=0, max_change=0.0; k<lx*ly; k++)
      max_change = MAX((mid[k].x-proj[k].x)*(mid[k].x-proj[k].x) +
		       (mid[k].y-proj[k].y)*(mid[k].y-proj[k].y),
		       max_change);
    if (iter % 10 == 0)
      printf("iter = %d, t = %e, delta_t = %e, max_change = %e\n",
	     iter, t, delta_t, max_change);
    
    /* When we get here, the integration step was accepted. */
    
    t += delta_t;
    iter++;
    for (k=0; k<lx*ly; k++) {
	proj[k].x = mid[k].x;
	proj[k].y = mid[k].y;
      }
    delta_t *= INC_AFTER_ACC;           /* Try a larger step size next time. */
  } while ((max_change > CONV_MAX_CHANGE && t < MAX_T && iter < MAX_ITER)
	   || t < MIN_T);

  /* Free memory. */
  
  fftw_destroy_plan(plan_rho);
  fftw_destroy_plan(plan_gridvx);
  fftw_destroy_plan(plan_gridvy);
  fftw_free(rho);
  fftw_free(gridvx);
  fftw_free(gridvy);
  free(eul);
  free(mid);
  free(vx_intp);
  free(vy_intp);
  free(vx_intp_half);
  free(vy_intp_half);
  
  return;
}
