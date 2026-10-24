/* STREAM-style triad a[i] = b[i] + s*c[i] on every MPI rank at once, to measure the memory bandwidth a
   fully populated node can sustain (reference for the solver's memory traffic). Each rank uses three
   N-element double arrays; bandwidth is counted as 24 bytes per element (STREAM convention: two reads,
   one write; build with streaming stores so the write needs no read-for-ownership).
   usage: mpirun -np 104 ./stream_triad [N (default 40M)] [reps (default 10)] */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const long n = argc > 1 ? atol(argv[1]) : 40000000L;
    const int reps = argc > 2 ? atoi(argv[2]) : 10;
    double* a = aligned_alloc(64, n * sizeof(double));
    double* b = aligned_alloc(64, n * sizeof(double));
    double* c = aligned_alloc(64, n * sizeof(double));
    for (long i = 0; i < n; ++i) { a[i] = 0.0; b[i] = 1.0; c[i] = 2.0; }
    const double s = 3.0;
    double best = 1e30;
    for (int r = 0; r < reps; ++r)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        const double t0 = MPI_Wtime();
        for (long i = 0; i < n; ++i) a[i] = b[i] + s * c[i];
        MPI_Barrier(MPI_COMM_WORLD);
        const double t = MPI_Wtime() - t0;
        if (t < best) best = t;
    }
    double check = 0.0;
    for (long i = 0; i < n; i += 4096) check += a[i];
    const double bytes = 24.0 * (double)n * size;   /* all ranks, one triad pass */
    double total;
    MPI_Reduce(&check, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
        printf("STREAM_TRIAD ranks %d n %ld best_s %.6f node_GBps %.1f per_rank_GBps %.2f check %.1f\n",
               size, n, best, bytes / best / 1e9, bytes / best / 1e9 / size, total);
    free(a); free(b); free(c);
    MPI_Finalize();
    return 0;
}
