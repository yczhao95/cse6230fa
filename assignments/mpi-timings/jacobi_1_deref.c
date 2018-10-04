
void
jacobi_sweep (const int bl[], const double *restrict *restrict *restrict F, double *restrict *restrict *restrict U)
{
  #pragma omp parallel for collapse(1) schedule(static)
  for (int i = 0; i < bl[0]; i++) {
    for (int j = 0; j < bl[1]; j++) {
      const double *restrict Fxy = F[i][j];
      double *restrict Uxy[3][3];

      for (int d = 0; d < 3; d++) {
        for (int e = 0; e < 3; e++) {
          Uxy[d][e] = U[i + d - 1][j + e - 1];
        }
      }
      for (int k = 0; k < bl[2]; k++) {
        Uxy[1][1][k] = (1./6.) * (Fxy[k] + Uxy[0][1][k] +
                                           Uxy[2][1][k] +
                                           Uxy[1][0][k] +
                                           Uxy[1][2][k] +
                                           Uxy[1][1][k-1] +
                                           Uxy[1][1][k+1]);
      }
    }
  }
}
