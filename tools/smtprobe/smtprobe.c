/* SMT feasibility probe: fixed work per process, prints elapsed seconds.
   usage: smtprobe triad|gather <elements> <reps>
   triad : a[i] = b[i] + s*c[i]                (streaming, bandwidth-bound)
   gather: s += x[idx[i]], idx local-random within +-4096 of i (latency-bound,
           like a renumbered unstructured mesh) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
int main(int c, char **v){
  if (c < 4) return 2;
  long n = atol(v[2]); int reps = atoi(v[3]);
  double *a = malloc(n*sizeof(double)), *b = malloc(n*sizeof(double)), *x = malloc(n*sizeof(double));
  int *idx = malloc(n*sizeof(int));
  unsigned s = 12345u;
  for (long i = 0; i < n; i++){ a[i]=0; b[i]=i*1e-9; x[i]=i*1e-9; s = s*1103515245u+12345u; long j = i + (long)(s%8193) - 4096; if (j<0) j=0; if (j>=n) j=n-1; idx[i]=(int)j; }
  double sum = 0, t0 = now();
  for (int r = 0; r < reps; r++){
    if (!strcmp(v[1],"triad")) { for (long i = 0; i < n; i++) a[i] = b[i] + 1.000001*x[i]; sum += a[r%n]; }
    else { for (long i = 0; i < n; i++) sum += x[idx[i]]; }
  }
  printf("%.4f %g\n", now()-t0, sum);
  return 0;
}
