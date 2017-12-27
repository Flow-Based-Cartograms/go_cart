/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cartogram.h"

/**************************** Function prototypes. ***************************/

double min4 (double a, double b, double c, double d);
double max4 (double a, double b, double c, double d);
POINT affine_transf(int triid, POINT *tri, double x, double y);

/*****************************************************************************/
/********** Function to project the polygons in the input .gen file. *********/

void project (BOOLEAN proj_graticule)
{
  double *xdisp, x2, *ydisp, y2;
  int i, j;

  /* The displacement vector (xdisp[i*ly+j], ydisp[i*ly+j]) is the point     */
  /* that was initially at (i+0.5, j+0.5). We work with (xdisp, ydisp)       */
  /* instead of (proj.x, proj.y) so that we can use the function interpol()  */
  /* defined in integrate.c.                                                 */

  xdisp = (double*) malloc(lx * ly * sizeof(double));
  ydisp = (double*) malloc(lx * ly * sizeof(double));
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      xdisp[i*ly + j] = proj[i*ly + j].x - i - 0.5;
      ydisp[i*ly + j] = proj[i*ly + j].y - j - 0.5;
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
  if (proj_graticule)

    /****************** Project proj2 on the basis of proj. ******************/

    for (i=0; i<lx*ly; i++) {
      x2 = proj2[i].x;
      y2 = proj2[i].y;
      proj2[i].x = interpol(x2, y2, xdisp, "x") + x2;
      proj2[i].y = interpol(x2, y2, ydisp, "y") + y2;
    }
  
  /******************************* Free memory. ******************************/

  free(xdisp);
  free(ydisp);

  return;
}

/*****************************************************************************/
/* Function to return the maximum absolute relative area error. The relative */
/* area error is defined by:                                                 */
/* area_on_cartogram / target_area - 1.                                      */
/* The function also updates the arrays cart_area[] and area_err[] that are  */
/* passed by reference.                                                      */

double max_area_err (double *area_err, double *cart_area)
{
  double max, obj_area, sum_cart_area, sum_target_area;
  int i, j;
  
  for (i=0; i<n_reg; i++) {
    cart_area[i] = 0.0;
    for (j=0; j<n_polyinreg[i]; j++)
      cart_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
				   cartcorn[polyinreg[i][j]]);
  }
  for (i=0, sum_target_area=0.0; i<n_reg; i++)
    sum_target_area += target_area[i];
  for (i=0, sum_cart_area=0.0; i<n_reg; i++)
    sum_cart_area += cart_area[i];
  for (i=0; i<n_reg; i++) {
    obj_area =                         /* Objective area in cartogram units. */
      target_area[i] * sum_cart_area / sum_target_area;
    area_err[i] = cart_area[i] / obj_area - 1.0;
  }
  max = 0.0;                   /* Determine the maximum absolute area error. */
  for (i=0; i<n_reg; i++)
    max = MAX(max, fabs(area_err[i]));
  
  return max;
}

/*****************************************************************************/
/** Function to write cartogram polygons and relative area errors to files. **/

void output_to_ascii (void)
{
  FILE *err_file = fopen("area_error.dat", "w"),
    *gen_file = fopen("cartogram.gen", "w");
  int i, j;

  /***************** Output of coordinates to cartogram.gen. *****************/

  for (i=0; i<n_poly; i++) {
    fprintf(gen_file, "%d\n", polygon_id[i]);
    for (j=0; j<n_polycorn[i]; j++)
      fprintf(gen_file, "%f %f\n", cartcorn[i][j].x, cartcorn[i][j].y);
    fprintf(gen_file, "END\n");
  }
  fprintf(gen_file, "END\n");
  fclose(gen_file);

  /*************** Output of relative area errors to err_file. ***************/
  for (i=0; i<n_reg; i++) {
    fprintf(err_file, "region %d: ", region_id[i]);
    fprintf(err_file,
    	    "cartogram area = %f, relative error = %f\n",
	    cart_area[i], area_err[i]);
  }
  fclose(err_file);
  
  return;
}

/*****************************************************************************/
/* How do we compute the inverse projection? We partition the untransformed  */
/* lattice into right triangles. The following figure shows where the        */
/* vertices are and how we label the triangles.                              */
/*                                                                           */
/* ly   -------------------------------------------  ----------------------  */
/*  |  |\       4ly-1      /|\       8ly-1      /|    |\     4lx*ly-1     /| */
/*  |  |   \            /   |   \            /   |    |   \            /   | */
/*  |  |      \      /      |      \      /      |    |      \      /      | */
/* ly- | 4ly-3   \/   4ly-2 | 8ly-3   \/   8ly-2 |... | 4lx*ly  \/  4lx*ly | */
/* 0.5 |         /\         |         /\         |    | -3      /\  -2     | */
/*  |  |      /      \      |      /      \      |    |      /      \      | */
/*  |  |   /            \   |   /            \   |    |   /            \   | */
/*  |  |/       4ly-4      \|/       8ly-4      \|    |/     4lx*ly-4     \| */
/* ly  |--------------------|--------------------|-  -|--------------------| */
/* -1  |          .         |          .         |    |          .         | */
/*  .             .                    .                         .           */
/*  .  |          .         |          .         |    |          .         | */
/*  2  |--------------------|--------------------|-  -|--------------------| */
/*  |  |\         7        /|\       4ly+7      /|    |\    4(lx-1)ly+7   /| */
/*  |  |   \            /   |   \            /   |    |   \            /   | */
/*  |  |      \      /      |      \      /      |    |      \      /      | */
/* 1.5 |  5      \/      6  | 4ly+5   \/   4ly+6 |... | 4(lx-1) \/ 4(lx-1) | */
/*  |  |         /\         |         /\         |    | *ly + 5 /\ *ly + 6 | */
/*  |  |      /      \      |      /      \      |    |      /      \      | */
/*  |  |   /            \   |   /            \   |    |   /            \   | */
/*  |  |/         4        \|/       4ly+4      \|    |/    4(lx-1)ly+4   \| */
/*  1  |--------------------|--------------------|-  -|--------------------| */
/*  |  |\         3        /|\       4ly+3      /|    |\    4(lx-1)ly+3   /| */
/*  |  |   \            /   |   \            /   |    |   \            /   | */
/*  |  |      \      /      |      \      /      |    |      \      /      | */
/* 0.5 |  1      \/      2  | 4ly+1   \/   4ly+2 |... | 4(lx-1) \/ 4(lx-1) | */
/*  |  |         /\         |         /\         |    | *ly + 1 /\ *ly + 2 | */
/*  |  |      /      \      |      /      \      |    |      /      \      | */
/*  |  |   /            \   |   /            \   |    |   /            \   | */
/*  |  |/         0        \|/        4ly       \|    |/     4(lx-1)ly    \| */
/*  0   -------------------------------------------  ----------------------  */
/*     0 ------- 0.5 ------ 1 ------- 1.5 ------ 2 . lx-1 --- lx-0.5 ---- lx */
/*                                                                           */
/* We project all the vertices on this lattice, bearing in mind that we      */
/* already have saved almost half of these coordinates in proj[][].          */
/* Suppose we find that, after the cartogram transformation, a point (x, y)  */
/* is in the projected triangle (a, b, c). We want to find its position in   */
/* the original triangle (p, q, r). We locally approximate the cartogram     */
/* transformation by an affine transformation T such that T(a) = p,          */
/* T(b) = q and T(c) = r. We can think of T as a 3x3 matrix                  */
/*  /t11 t12 t13\                                                            */
/* | t21 t22 t23 |  such that                                                */
/*  \ 0   0   1 /                                                            */
/*  /t11 t12 t13\   /a1 b1 c1\     /p1 q1 r1\                                */
/* | t21 t22 t23 | | a2 b2 c2 | = | p2 q2 r2 | or TA = P. Hence T = PA^{-1}. */
/*  \ 0   0   1 /   \ 1  1  1/     \ 1  1  1/                                */
/*                              /b2-c2 c1-b1 b1*c2-b2*c1\                    */
/* We have A^{-1} = (1/det(A)) | c2-a2 a1-c1 a2*c1-a1*c2 |. By multiplying   */
/*                              \a2-b2 b1-a1 a1*b2-a2*b1/                    */
/* PA^{-1} we obtain t11, t12, t13, t21, t22, t23. The preimage of (x, y) in */
/* the unprojected map is then "pre" with coordinates                        */
/* pre.x = t11*x + t12*y + t13, pre.y = t21*x + t22*y + t23.                 */

/***************** Function to perform the affine transform. *****************/
/* Input: triid - ID of the triangle, see figure above.                      */
/*        tri   - the vertices of the cartogram-transformed triangle.        */
/*        x, y  - coordinates on the cartogram.                              */

POINT affine_transf(int triid, POINT *tri, double x, double y)
{
  double ainv11, ainv12, ainv13, ainv21, ainv22, ainv23, ainv31, ainv32,
    ainv33, t11, t12, t13, t21, t22, t23, det;
  POINT p, pre, q, r;

  /* Determine the vertices p, q, r of the unprojected triangle from the ID  */
  /* of the triangle. Note that the order of the three points must match the */
  /* order of the vertices in tri[].                                         */
  
  switch (triid % 4) {
  case 0:
    p.x = triid / (4 * ly);
    p.y = (triid / 4) % ly;
    q.x = p.x + 0.5;
    q.y = p.y + 0.5;
    r.x = p.x + 1;
    r.y = p.y;
    break;
  case 1:
    p.x = triid / (4 * ly);
    p.y = (triid / 4) % ly;
    q.x = p.x;
    q.y = p.y + 1;
    r.x = p.x + 0.5;
    r.y = p.y + 0.5;
    break;
  case 2:
    p.x = triid / (4 * ly) + 0.5;
    p.y = (triid / 4) % ly + 0.5;
    q.x = p.x + 0.5;
    q.y = p.y + 0.5;
    r.x = q.x;
    r.y = q.y - 1;
    break;
  default:
    p.x = triid / (4 * ly);
    p.y = (triid / 4) % ly + 1;
    q.x = p.x + 1;
    q.y = p.y;
    r.x = p.x + 0.5;
    r.y = p.y - 0.5;
  }

  /**************************** Determinant of A. ****************************/
  
  det = tri[0].x * tri[1].y + tri[1].x * tri[2].y + tri[2].x * tri[0].y
    - tri[1].x * tri[0].y - tri[2].x * tri[1].y - tri[0].x * tri[2].y;
  
  /*********** Compute det(A) * A^{-1}. We divide by det(A) later. ***********/
  
  ainv11 = tri[1].y - tri[2].y;
  ainv12 = tri[2].x - tri[1].x;
  ainv13 = tri[1].x * tri[2].y - tri[1].y * tri[2].x;
  ainv21 = tri[2].y - tri[0].y;
  ainv22 = tri[0].x - tri[2].x;
  ainv23 = tri[0].y * tri[2].x - tri[0].x * tri[2].y;
  ainv31 = tri[0].y - tri[1].y;
  ainv32 = tri[1].x - tri[0].x;
  ainv33 = tri[0].x * tri[1].y - tri[0].y * tri[1].x;

  /******************************** Compute T. *******************************/
  
  t11 = p.x * ainv11 + q.x * ainv21 + r.x * ainv31;
  t12 = p.x * ainv12 + q.x * ainv22 + r.x * ainv32;
  t13 = p.x * ainv13 + q.x * ainv23 + r.x * ainv33;
  t21 = p.y * ainv11 + q.y * ainv21 + r.y * ainv31;
  t22 = p.y * ainv12 + q.y * ainv22 + r.y * ainv32;
  t23 = p.y * ainv13 + q.y * ainv23 + r.y * ainv33;

  /********************* Transform the input coordinates. ********************/
  
  pre.x = (t11*x + t12*y + t13) / det;
  pre.y = (t21*x + t22*y + t23) / det;
  
  return(pre);
}

/*************** Helper function: return min/max of 4 numbers. ***************/
double min4 (double a, double b, double c, double d)
{
  if (a <= b && a <= c && a <= d)
    return(a);
  if (b <= a && b <= c && b <= d)
    return(b);
  if (c <= a && c <= b && c <= d)
    return(c);
  return(d);
}
double max4 (double a, double b, double c, double d)
{
  if (a >= b && a >= c && a >= d)
    return(a);
  if (b >= a && b >= c && b >= d)
    return(b);
  if (c >= a && c >= b && c >= d)
    return(c);
  return(d);
}

/********** Function to calculate the inverse projection invproj[]. **********/

void inv_project (void)
{
  double *xdisp, *ydisp;
  int i, j, k, **xyhalfshift2tri;
  POINT *invproj, *invproj2, **projgrid, **tri;
  
  /**************************** Memory allocation. ***************************/

  xdisp = (double*) malloc(lx * ly * sizeof(double));
  ydisp = (double*) malloc(lx * ly * sizeof(double));
  invproj = (POINT*) malloc(lx * ly * sizeof(POINT));
  invproj2 = (POINT*) malloc(lx * ly * sizeof(POINT));
  projgrid = (POINT**) malloc((lx+1) * sizeof(POINT*));
  for (i=0; i<=lx; i++)
    projgrid[i] = (POINT*) malloc((ly+1) * sizeof(POINT));
  tri = (POINT**) malloc(4 * lx * ly * sizeof(POINT*));
  for (i=0; i<4*lx*ly; i++)
    tri[i] = (POINT*) malloc(3 * sizeof(POINT));
  xyhalfshift2tri = (int**) malloc(lx * sizeof(int*));
  for (i=0; i<lx; i++)
    xyhalfshift2tri[i] = (int*) malloc(ly * sizeof(int));
  
  /* The displacement vector (xdisp[i*ly+j], ydisp[i*ly+j]) is the point     */
  /* that was initially at (i+0.5, j+0.5). We work with (xdisp, ydisp)       */
  /* instead of (proj.x, proj.y) so that we can use the function interpol()  */
  /* defined in integrate.c.                                                 */
  
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      xdisp[i*ly + j] = proj[i*ly + j].x - i - 0.5;
      ydisp[i*ly + j] = proj[i*ly + j].y - j - 0.5;
    }
  
  /* projgrid[i][j] is the projected position of (i, j) without half-shift.  */
  
  for (i=0; i<=lx; i++)
    for (j=0; j<=ly; j++) {      
      projgrid[i][j].x = interpol(i, j, xdisp, "x") + i;      
      projgrid[i][j].y = interpol(i, j, ydisp, "y") + j;
    }
  
  /************ Project the triangles shown in the lattice above. ************/
  
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      tri[4*(i*ly + j)][0].x =                      /* Lower left of square. */
	tri[4*(i*ly + j) + 1][0].x = projgrid[i][j].x;
      tri[4*(i*ly + j)][0].y =
	tri[4*(i*ly + j) + 1][0].y = projgrid[i][j].y;
      tri[4*(i*ly + j) + 1][1].x =                            /* Upper left. */
	tri[4*(i*ly + j) + 3][0].x = projgrid[i][j+1].x;
      tri[4*(i*ly + j) + 1][1].y =
	tri[4*(i*ly + j) + 3][0].y = projgrid[i][j+1].y;
      tri[4*(i*ly + j)][2].x =                               /* Lower right. */
	tri[4*(i*ly + j) + 2][2].x = projgrid[i+1][j].x;
      tri[4*(i*ly + j)][2].y =
	tri[4*(i*ly + j) + 2][2].y = projgrid[i+1][j].y;
      tri[4*(i*ly + j) + 2][1].x =                           /* Upper right. */
	tri[4*(i*ly + j) + 3][1].x = projgrid[i+1][j+1].x;
      tri[4*(i*ly + j) + 2][1].y =
	tri[4*(i*ly + j) + 3][1].y = projgrid[i+1][j+1].y;
      tri[4*(i*ly + j)][1].x =                                  /* Midpoint. */
	tri[4*(i*ly + j) + 1][2].x =
	tri[4*(i*ly + j) + 2][0].x =
	tri[4*(i*ly + j) + 3][2].x = proj[i*ly + j].x;
      tri[4*(i*ly + j)][1].y =
	tri[4*(i*ly + j) + 1][2].y =
	tri[4*(i*ly + j) + 2][0].y =
	tri[4*(i*ly + j) + 3][2].y = proj[i*ly + j].y;
    }

  /***** xyhalfshift2tri[i][j]=k means that (i+0.5, j+0.5) is in tri[k]. *****/
  
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++)
      xyhalfshift2tri[i][j] = -1;
  for (i=0; i<4*lx*ly; i++)
    set_inside_values_for_polygon(i, 3, tri[i], xyhalfshift2tri);

  /**** Inverse projection for a point at (i+0.5, j+0.5) on the cartogram. ***/
  
  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      k = xyhalfshift2tri[i][j];
      invproj[i*ly + j] = affine_transf(k, tri[k], i+0.5, j+0.5);
    }

  /* Near sharp density gradients, there can be numerical artifacts. We      */
  /* polish them here. If the preimage of (i+0.5, j+0.5) has an x-coordinate */
  /* less than                                                               */
  /* min(invproj[i*ly + j - 1].x, invproj[i*ly + j + 1].x,                   */
  /*     invproj[(i-1)*ly + j].x, invproj[(i+1)*ly + j].x) - 1,              */
  /* greater than                                                            */
  /* max(invproj[i*ly + j - 1].x, invproj[i*ly + j + 1].x,                   */
  /*     invproj[(i-1)*ly + j].x, invproj[(i+1)*ly + j].x) + 1,              */
  /* a y-coordinate less than                                                */
  /* min(invproj[i*ly + j - 1].y, invproj[i*ly + j + 1].y,                   */
  /*     invproj[(i-1)*ly + j].y, invproj[(i+1)*ly + j].y) - 1               */
  /* or greater than                                                         */
  /* max(invproj[i*ly + j - 1].y, invproj[i*ly + j + 1].y,                   */
  /*     invproj[(i-1)*ly + j].y, invproj[(i+1)*ly + j].y) + 1,              */
  /* we replace invproj[i*ly + j] by the centroid of the four neighbouring   */
  /* preimages.                                                              */
  
  for (j=0; j<ly-1; j++) {
    invproj2[j].x = invproj[j].x;
    invproj2[j].y = invproj[j].y;
  }
  for (i=0; i<lx-1; i++) {
    invproj2[i*ly + ly - 1].x = invproj[i*ly + ly - 1].x;
    invproj2[i*ly + ly - 1].y = invproj[i*ly + ly - 1].y;
  }
  for (j=1; j<ly; j++) {
    invproj2[(lx-1)*ly + j].x = invproj[(lx-1)*ly + j].x;
    invproj2[(lx-1)*ly + j].y = invproj[(lx-1)*ly + j].y;    
  }
  for (i=1; i<lx; i++) {
    invproj2[i*ly].x = invproj[i*ly].x;
    invproj2[i*ly].y = invproj[i*ly].y;
  }  
  for (i=1; i<lx-1; i++)
    for (j=1; j<ly-1; j++) {      
      if (invproj[i*ly + j].x < min4(invproj[i*ly + j - 1].x,
				     invproj[i*ly + j + 1].x,
				     invproj[(i-1)*ly + j].x,
				     invproj[(i+1)*ly + j].x) - 1 ||
	  invproj[i*ly + j].x > max4(invproj[i*ly + j - 1].x,
				     invproj[i*ly + j + 1].x,
				     invproj[(i-1)*ly + j].x,
				     invproj[(i+1)*ly + j].x) + 1 ||
	  invproj[i*ly + j].y < min4(invproj[i*ly + j - 1].y,
				     invproj[i*ly + j + 1].y,
				     invproj[(i-1)*ly + j].y,
				     invproj[(i+1)*ly + j].y) - 1 ||
	  invproj[i*ly + j].y > max4(invproj[i*ly + j - 1].y,
				     invproj[i*ly + j + 1].y,
				     invproj[(i-1)*ly + j].y,
				     invproj[(i+1)*ly + j].y) + 1) {
	invproj2[i*ly + j].x =
	  0.25 * (invproj[i*ly + j - 1].x + invproj[i*ly + j + 1].x +
		  invproj[(i-1)*ly + j].x + invproj[(i+1)*ly + j].x);
	invproj2[i*ly + j].y =
	  0.25 * (invproj[i*ly + j - 1].y + invproj[i*ly + j + 1].y +
		  invproj[(i-1)*ly + j].y + invproj[(i+1)*ly + j].y);
      }
      else {
	invproj2[i*ly + j].x = invproj[i*ly + j].x;
	invproj2[i*ly + j].y = invproj[i*ly + j].y;
      }
    }
  
  /****** Represent the inverse projection by an image of the graticule. *****/
  
  ps_figure("invproj.eps", origcorn, invproj2, TRUE);
  
  /******************************* Free memory. ******************************/
  
  free(xdisp);
  free(ydisp);
  free(invproj);
  free(invproj2);
  for (i=0; i<=lx; i++)
    free(projgrid[i]);
  free(projgrid);
  for (i=0; i<4*lx*ly; i++)
    free(tri[i]);
  free(tri);
  for (i=0; i<lx; i++)
    free(xyhalfshift2tri[i]);
  free(xyhalfshift2tri);

  return;
}
