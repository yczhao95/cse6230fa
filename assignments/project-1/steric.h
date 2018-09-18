#if !defined(STERIC_H)
#define      STERIC_H


/* This kernel should be called if the distance between two particles is less
 * than twice the particle radius */
static inline void
force_in_range (double k, /* The interaction strength (may be scaled by the time step already) */
                double r, /* The radius of a particle.  Two particles interact if they intersect */
                double R, /* The distance between these two particles */
                double dx, double dy, double dz, /* The displacement from particle 2 to particle 1 */
                double f[3]) /* The output force exerted on particle 1 by particle 2 */
{
  /* The interaction strength starts at 0 when they are just touching,
   * becoming infinite as the distance becomes zero */
  double strength = (2. * r - R) / (2. * R);

  f[0] = k * strength * dx;
  f[1] = k * strength * dy;
  f[2] = k * strength * dz;
}

static inline void
force (double k, /* The interaction strength (may be scaled by the time step already) */
       double r, /* The radius of a particle.  Two particles interact if they intersect */
       double L, /* The width of the periodic domain */
       double x1, double y1, double z1, /* The center of the first particle */
       double x2, double y2, double z2, /* The center of the first particle */
       double f[3]) /* The output force exterted on particle 1 by particle 2 */
{
  double dx, dy, dz;
  double R2, ir3;
  double interactR2 = 4. * r * r;

  dx = remainder(x1 - x2, L);
  dy = remainder(y1 - y2, L);
  dz = remainder(z1 - z2, L);

  R2 = dx*dx + dy*dy + dz*dz;

  /* If the distance between the centers is less than twice the radius, they
   * interact */
  if (R2 < interactR2) {
    double R = sqrt(R2);

    force_in_range (k, r, R, dx, dy, dz, f);
  }
}


#endif
