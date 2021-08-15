/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "adi-sun-kernel.h"
#include "adi-util.h"

namespace adi {
namespace cuda {
namespace sun {

__global__ void CalcEciPosition (
  Vector*   position,
  int64_t*  ticks,
  size_t    N)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= N)
  {
    return;
  }
  const double mjd = ToJ2000(ticks[tid]);
  double delta_ET;
  {
    const double year = 1900 + mjd / 365.25;
    delta_ET = 26.465 + 0.747622 * (year - 1950) + 1.886913
          * sin (kTWOPI * (year - 1975) / 33);
  }
  const double T = (mjd + delta_ET / kSECONDS_PER_DAY) / 36525.0;
  const double M = DegreesToRadians (Wrap360 (358.47583
              + Wrap360 (35999.04975 * T)
              - (0.000150 + 0.0000033 * T) * T * T));
  const double L = DegreesToRadians (Wrap360 (279.69668
              + Wrap360 (36000.76892 * T)
              + 0.0003025 * T * T));
  const double e = 0.01675104 - (0.0000418 + 0.000000126 * T) * T;
  const double C = DegreesToRadians ((1.919460
              - (0.004789 + 0.000014 * T) * T) * sin (M)
              + (0.020094 - 0.000100 * T) * sin (2 * M)
              + 0.000293 * sin (3 * M));
  const double O = DegreesToRadians (Wrap360 (259.18 - 1934.142 * T));
  const double Lsa = WrapTwoPI (L + C - DegreesToRadians (0.00569 - 0.00479 * sin (O)));
  const double nu = WrapTwoPI (M + C);
  double R = 1.0000002 * (1 - e * e) / (1 + e * cos (nu));
  const double eps = DegreesToRadians (23.452294 - (0.0130125 + (0.00000164 - 0.000000503 * T) * T) * T + 0.00256 * cos(O));
  R = R * kAU;
  position[tid].x = R * cos (Lsa);
  position[tid].y = R * sin (Lsa) * cos (eps);
  position[tid].z = R * sin (Lsa) * sin (eps);
}

}
}
}