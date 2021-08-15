/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2021 Innovation Academy for Microsatellites of CAS
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Wang Junyong (wangjunyong@microsate.com)
 */

#include "adi-satellite.h"
#include "adi-constant.h"

using namespace adi;

int main ()
{
  DateTime epoch = DateTime (2020, 1, 1, 0, 0, 0);
  DateTime start = DateTime (2020, 1, 1, 0, 0, 0);
  DateTime stop = DateTime (2020, 1, 2, 0, 0, 0);
  TimeSpan minute = TimeSpan (0, 1, 0);
  TimeSpan second = TimeSpan (0, 0, 1);
  Intervals intervals;
  Interval interval = Interval (start, stop);
  interval.InsertToList (intervals);
  interval.Shift (Day).Shift (Day);
  interval.InsertToList (intervals);
  double sma = K_RE + 500.0;
  double ecc = 0.0;
  double inc = 97.0 * M_PI / 180.0;
  double raan = 100.0 * M_PI / 180.0;
  double aop = 0.0;
  double ma = 0.0;
  Satellite testSatellite;
  testSatellite.SetEpoch (epoch);
  testSatellite.SetStep (Second);
  testSatellite.SetElement (sma, ecc, inc, raan, aop, ma, false);
  testSatellite.SetIntervalList (intervals);
  clock_t cpuStart, cpuStop;
  cpuStart = clock ();
  testSatellite.CalcTrajectory ();
  cpuStop = clock ();
  std::cout << (double)(cpuStop - cpuStart) / CLOCKS_PER_SEC << "sec";
  return 0;
}