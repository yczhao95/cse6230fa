#if !defined(CLOUD_H)
#define      CLOUD_H
#include <math.h>
#include <cse6230rand.h>

/* a quadratic modification of periodic distance,
 * with a periodicity of 1 */
#define qdist(x) ((x) * (1. - fabs(x)))

/* The derivative of qdist wrt x */
#define dqdist(x) (1. - 2. * fabs(x))

/* computes the potential between two particles,
 * assuming they are moving in a periodic domain.
 * This model problem asssumes both masses are 1.
 * To make the potential smooth, we smooth out the
 * elbow in the periodic distance function with
 * a quadratic approximation */
static inline double
potential (double k,
           double x1, double y1, double z1,
           double x2, double y2, double z2)
{
  double dx, dy, dz;
  double r2, ir;

#if 1
  dx = remainder(x1 - x2, 1.);
  dy = remainder(y1 - y2, 1.);
  dz = remainder(z1 - z2, 1.);

  dx = qdist(dx);
  dy = qdist(dy);
  dz = qdist(dz);
#else
  dx = x1 - x2;
  dy = y1 - y2;
  dz = z1 - z2;
#endif

  r2 = dx*dx + dy*dy + dz*dz;
  ir = 1. / sqrt(r2);

  return k * ir;
}

/* computes the force particle 2 exerts on particle 1,
 * according to the potential above */
static inline void
force (double k,
       double x1, double y1, double z1,
       double x2, double y2, double z2,
       double f[])
{
  double dx, dy, dz;
  double qdx, qdy, qdz;
  double r2, ir3;

#if 1
  dx = remainder(x1 - x2, 1.);
  dy = remainder(y1 - y2, 1.);
  dz = remainder(z1 - z2, 1.);

  qdx = qdist(dx);
  qdy = qdist(dy);
  qdz = qdist(dz);
#else
  qdx = x1 - x2;
  qdy = y1 - y2;
  qdz = z1 - z2;
#endif

  r2 = qdx*qdx + qdy*qdy + qdz*qdz;
  ir3 = 1. / (sqrt(r2) * r2);

#if 1
  f[0] = k * ir3 * qdx * dqdist(dx);
  f[1] = k * ir3 * qdy * dqdist(dy);
  f[2] = k * ir3 * qdz * dqdist(dz);
#else
  f[0] = k * ir3 * qdx;
  f[1] = k * ir3 * qdy;
  f[2] = k * ir3 * qdz;
#endif
}

#endif
