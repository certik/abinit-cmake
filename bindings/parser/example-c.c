#include "ab_dtset_c.h"
#include <stdlib.h>
#include <stdio.h>

static void onError(AbDtsets *dt, AbError error);

int main(int argc, const char *argv[])
{
  AbError error;
  AbDtsets *dt;
  int natom, ndtset, idtset;
  double *coord, rprimd[3][3];
  int dims[7], ndims;
  int i, j, n, brvltt;

  dt = ab_dtset_new(argv[1]);

  error = ab_dtset_get_ndtset(dt, &ndtset);
  if (error != AB_NO_ERROR) onError(dt, error);

  for (idtset = 0; idtset <= ndtset; idtset++)
    {
      error = ab_dtset_get_integer(dt, AB_DTSET_NATOM, idtset, &natom);
      if (error != AB_NO_ERROR) onError(dt, error);
      error = ab_dtset_get_shape(dt, &ndims, dims, AB_DTSET_XRED_ORIG, idtset);
      if (error != AB_NO_ERROR) onError(dt, error);
      n     = dims[0] * dims[1];
      coord = malloc(sizeof(double) * n);
      error = ab_dtset_get_real_array(dt, coord, n, AB_DTSET_XRED_ORIG, idtset);
      if (error != AB_NO_ERROR) onError(dt, error);
      error = ab_dtset_get_real_array(dt, (double*)rprimd, 9, AB_DTSET_RPRIMD_ORIG, idtset);
      if (error != AB_NO_ERROR) onError(dt, error);
      error = ab_dtset_get_integer(dt, AB_DTSET_BRVLTT, idtset, &brvltt);
      if (error != AB_NO_ERROR) onError(dt, error);

      printf("### DATASET %d/%d ###\n", idtset, ndtset);
      printf("Number of atoms in dataset %d: %d\n",
	     idtset, natom);
      printf("box definition: ( %f %f %f )\n",
	     rprimd[0][0], rprimd[0][1], rprimd[0][2]);
      printf("                ( %f %f %f )\n",
	     rprimd[1][0], rprimd[1][1], rprimd[1][2]);
      printf("                ( %f %f %f )\n",
	     rprimd[2][0], rprimd[2][1], rprimd[2][2]);
      printf("Size of coordiantes array in dataset %d: %d\n", idtset, n);
      printf("Coordinates in dataset %d:\n", idtset);
      for (j = 0; j < dims[1]; j++)
	{
	  for (i = 0; i < dims[0]; i++)
	    printf("  %g", coord[j * dims[0] + i]);
	  printf("\n");
	}
      free(coord);
      printf("Bravais lattice in dataset %d: %d\n", idtset, brvltt);
      printf("\n");
    }

  ab_dtset_free(dt);

  return 0;
}

static void onError(AbDtsets *dt, AbError error)
{
  fprintf(stderr, "Error %d\n", (int)error);
  ab_dtset_free(dt);
  exit(error);
}
