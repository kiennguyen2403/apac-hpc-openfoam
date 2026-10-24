/* Diagnostic only: LD_PRELOAD constructor that sets FTZ and DAZ in MXCSR
   (denormal results flushed to zero, denormal inputs read as zero). Not
   IEEE-compliant; used to test whether denormal arithmetic slows some ranks.
   FTZ_DAZ_OFF=1: leave MXCSR alone, only report.
   FTZ_DAZ_REPORT=1: rank 0 prints MXCSR at load and at exit. */
#include <xmmintrin.h>
#include <stdio.h>
#include <stdlib.h>

static int report(void)
{
    return getenv("FTZ_DAZ_REPORT") && getenv("PMI_RANK") && atoi(getenv("PMI_RANK")) == 0;
}

static void atExit(void)
{
    if (report()) fprintf(stderr, "ftzDaz: MXCSR at exit = 0x%x\n", _mm_getcsr());
}

__attribute__((constructor)) static void ftzDazOn(void)
{
    if (!getenv("FTZ_DAZ_OFF")) _mm_setcsr(_mm_getcsr() | 0x8040);
    if (report()) fprintf(stderr, "ftzDaz: MXCSR at load = 0x%x\n", _mm_getcsr());
    atexit(atExit);
}
