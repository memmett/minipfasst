/* Dump routines for VisIt BOV format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define BUFLEN 256

void dump_velocity_c(char *dname, char *fname, int nx, int ny, int nz, double *array)
{
  int    i, j, k;
  FILE*  fp;
  char   errmsg[BUFLEN];
  char   buf[BUFLEN];

  double *wbuf;

  /* bov file */
  snprintf(buf, BUFLEN, "%s/%s.bov", dname, fname);
  fp = fopen(buf, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create BOV file (%s)", buf);
    perror(errmsg);
    return;
  }

  /* dat file name */
  snprintf(buf, BUFLEN, "%s.dat", fname);

  fprintf(fp, "VARIABLE: velocity\n");
  fprintf(fp, "DATA_FILE: %s\n", buf);
  fprintf(fp, "DATA_SIZE: %d %d %d\n", nx, ny, nz);
  fprintf(fp, "DATA_COMPONENTS: 3\n");
  fprintf(fp, "DATA_FORMAT: DOUBLE\n");
  fprintf(fp, "DATA_ENDIAN: LITTLE\n");
  fclose(fp);

  /* dat file */
  snprintf(buf, BUFLEN, "%s/%s.dat", dname, fname);
  fp = fopen(buf, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create DAT file (%s)", buf);
    perror(errmsg);
    return;
  }

  wbuf = malloc(3*nx*ny*nz*sizeof(double));
  if (wbuf == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to allocate write buffer");
    fprintf(stderr, errmsg);
    perror(errmsg);
    return;
  }

  for (i=0; i<3*nx*ny*nz; i++)
    wbuf[i] = array[2*i];

  fwrite(wbuf, sizeof(double), 3*nx*ny*nz, fp);

  free(wbuf);
  fclose(fp);
}

void dump_vorticity_c(char *dname, char *fname, int nx, int ny, int nz, double *array)
{
  int    i, j, k;
  FILE*  fp;
  char   errmsg[BUFLEN];
  char   buf[BUFLEN];

  /* bov file */
  snprintf(buf, BUFLEN, "%s/%s.bov", dname, fname);
  fp = fopen(buf, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create BOV file (%s)", buf);
    perror(errmsg);
    return;
  }

  /* dat file name */
  snprintf(buf, BUFLEN, "%s.dat", fname);

  fprintf(fp, "VARIABLE: vorticity\n");
  fprintf(fp, "DATA_FILE: %s\n", buf);
  fprintf(fp, "DATA_SIZE: %d %d %d\n", nx, ny, nz);
  fprintf(fp, "DATA_COMPONENTS: 1\n");
  fprintf(fp, "DATA_FORMAT: DOUBLE\n");
  fprintf(fp, "DATA_ENDIAN: LITTLE\n");
  fclose(fp);

  /* dat file */
  snprintf(buf, BUFLEN, "%s/%s.dat", dname, fname);
  fp = fopen(buf, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create DAT file (%s)", buf);
    perror(errmsg);
    return;
  }

  fwrite(array, sizeof(double), nx*ny*nz, fp);

  fclose(fp);
}
