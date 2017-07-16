/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include "cartogram.h"

/******************** Function to prepare postscript map. ********************/

void ps_figure (char *ps_name, POINT **corn)
{
  FILE *ps_file;
  int i, j, k;
  
  if ((ps_file = fopen(ps_name, "w")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open ps-file.\n");
    exit(1);
  }
  
  /**************************** Postscript header. ***************************/
  
  fprintf(ps_file, "%%!PS-Adobe-2.0 EPSF-2.0\n%%%%Title: %s\n", ps_name);
  fprintf(ps_file, "%%%%Creator: Michael T. Gastner et al.\n");
  fprintf(ps_file,
          "%%%%BoundingBox: 0 0 %d %d\n%%%%Magnification: 1.0000\n", lx, ly);
  fprintf(ps_file, "%%%%EndComments\n");
  fprintf(ps_file, "/m {moveto} def\n/l {lineto} def\n/s {stroke} def\n");
  fprintf(ps_file, "/n {newpath} def\n/c {closepath} def\n/f {fill} def\n");
  fprintf(ps_file, "/SLW {setlinewidth} def\n/SGRY {setgray} def\n");
  fprintf(ps_file, "/SRGB {setrgbcolor} def\n");
  
  /**************************** Plot the regions. ****************************/
  
  fprintf(ps_file, "0.7 SLW\n");
  for (i=0; i<n_reg; i++) {
    fprintf(ps_file, "n\n");
    for (j=0; j<n_polyinreg[i]; j++) {
      fprintf(ps_file, "%.3f %.3f m\n",
	      corn[polyinreg[i][j]][0].x, corn[polyinreg[i][j]][0].y);
      for (k=0; k<n_polycorn[polyinreg[i][j]]; k++)
      	fprintf(ps_file, "%.3f %.3f l\n",
		corn[polyinreg[i][j]][k].x, corn[polyinreg[i][j]][k].y);
      fprintf(ps_file, "c\n");
    }
    fprintf(ps_file, "gsave\n0.5 SGRY f\ngrestore\n0 SGRY s\n");
  }
  fprintf(ps_file, "showpage\n");
  
  fclose(ps_file);
  
  return;
}
