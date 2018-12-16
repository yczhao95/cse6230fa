
void compute_thread_bounds (const int *restrict bl, const int *restrict t, int *restrict bounds);
void jacobi_sweep (const int *restrict bl, const int *restrict bounds,
                   const double *restrict *restrict *restrict F,
                   double *restrict *restrict *restrict U,
                   double *restrict *restrict *restrict Ucopy);
