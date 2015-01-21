/* Dump routines for VisIt BOV format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#define BUFLEN 256

void dump_mkdir(char *dname)
{
  if (mkdir(dname, 0755)) {
    if (errno != EEXIST) {
      char errmsg[BUFLEN];
      snprintf(errmsg, BUFLEN, "WARNING: Unable to create directory: %s", dname);
      perror(errmsg);
    }
  }
}

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

  wbuf = (double*) malloc(3*nx*ny*nz*sizeof(double));
  if (wbuf == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to allocate write buffer");
    perror(errmsg);
    return;
  }

  for (i=0; i<3*nx*ny*nz; i++)
    wbuf[i] = array[2*i];

  fwrite(wbuf, sizeof(double), 3*nx*ny*nz, fp);

  free(wbuf);
  fclose(fp);
}

void dump_velocity_component_c(char *fname, char *comp, int nx, int ny, int nz, double *u)
{
  int    i, j, k;
  FILE*  fp;
  char   errmsg[BUFLEN];
  char   buf[BUFLEN];

  double *wbuf;

  /* bov file */
  snprintf(buf, BUFLEN, "%s_%s.bov", fname, comp);
  fp = fopen(buf, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create BOV file (%s)", buf);
    perror(errmsg);
    return;
  }

  /* dat file name */
  snprintf(buf, BUFLEN, "%s_%s.dat", fname, comp);

  fprintf(fp, "VARIABLE: %s\n", comp);
  fprintf(fp, "DATA_FILE: %s\n", buf);
  fprintf(fp, "DATA_SIZE: %d\n", nx);
  fprintf(fp, "DATA_COMPONENTS: 1\n");
  fprintf(fp, "DATA_FORMAT: DOUBLE\n");
  fprintf(fp, "DATA_ENDIAN: LITTLE\n");
  fclose(fp);

  /* dat file */
  snprintf(buf, BUFLEN, "%s.dat", fname);
  fp = fopen(buf, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create DAT file (%s)", buf);
    perror(errmsg);
    return;
  }

  fwrite(u, sizeof(double), nx*ny*nz, fp);

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

void read_velocity_c(char *fname, int nx, int ny, int nz, double *array)
{
  int    i, j, k;
  FILE*  fp;
  char   errmsg[BUFLEN];

  double *wbuf;

  fp = fopen(fname, "rb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to open DAT file (%s)", fname);
    perror(errmsg);
    return;
  }

  wbuf = (double*) malloc(3*nx*ny*nz*sizeof(double));
  if (wbuf == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to allocate read buffer");
    perror(errmsg);
    return;
  }

  fread(wbuf, sizeof(double), 3*nx*ny*nz, fp);

  for (i=0; i<3*nx*ny*nz; i++) {
    array[2*i]   = wbuf[i];
    array[2*i+1] = 0.0;
  }

  free(wbuf);
  fclose(fp);
}
