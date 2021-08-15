#include "adi.h"

int main ()
{
  adi::DateTime epoch = adi::DateTime (2020, 1, 1, 0, 0, 0);
  adi::DateTime start = adi::DateTime (2020, 1, 1, 0, 0, 0);
  adi::DateTime stop = adi::DateTime (2020, 1, 7, 0, 0, 0);
  adi::TimeSpan minute = adi::TimeSpan (0, 1, 0);
  adi::TimeSpan second = adi::TimeSpan (0, 0, 1);
  Interval interval;
  interval.epoch = epoch;
  interval.start = start;
  interval.stop = stop;
  interval.step = minute;
  double sma = K_RE + 500.0;
  double ecc = 0.0;
  double inc = adi::DegreesToRadians (97.0);
  double raan = adi::DegreesToRadians (10.0);
  double aop = 0.0;
  double ma = 0.0;
  size_t N_SAT = 10;
  adi::SatelliteHelper satHelper;
  adi::Satellites plane = satHelper.DistributeEvenly (adi::Satellite::Element {sma, ecc, inc, raan, aop, ma}, N_SAT);

  satHelper.CalcAllTrajectories (interval);
  return 0;
}
