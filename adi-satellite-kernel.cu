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

#include "adi-satellite-kernel.h"
#include "adi-util.h"
#include "adi-constant.h"

namespace adi {
namespace cuda {
namespace satellite {

__global__ void CalcParam(
  adi::Satellite::Par*  param,
  size_t                rows)
{
  uint32_t sid = blockIdx.x;
  if (sid >= rows)
  {
    return;
  }
  adi::Satellite::Ele& elem = param[sid].elem;
  double tmp1 = (1 - elem.ecc * elem.ecc);
  double tmp2 = tmp1 * tmp1;
  double tmp3 = K_RE * K_RE / (elem.sma * elem.sma);
  double tmp4 = cos (elem.inc);
  tmp4 = tmp4 * tmp4;
  double tmp5 = tmp1 * sqrt (tmp1);
  // Calculate the angular velocity
  param[sid].ar = sqrt (K_MU / (elem.sma * elem.sma * elem.sma));
  // Calculate the semi-latus rectum
  param[sid].slr = elem.sma * tmp1;
  // Calculate the rate of right ascension of ascending node
  param[sid].draan = -1.5 * param[sid].ar * K_J2 * cos (elem.inc) * tmp3 / tmp2;
  // Calculate the rate of argument of perigee
  param[sid].daop = -0.75 * param[sid].ar * K_J2 * (1 - 5 * tmp4) * tmp3 / tmp2;
  // Calculate the rate of mean anomaly
  param[sid].dma = param[sid].ar + 0.75 * param[sid].ar * K_J2 * (3 * tmp4 - 1) * tmp3 / tmp5;
  return;
}

__global__ void CalcState (
  adi::Satellite::Sta*  state,
  adi::Satellite::Par*  param,
  double*               times,
  size_t                rows,
  size_t                cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  double t = times[tid];
  adi::Satellite::Ele& elem = param[sid].elem;
  uint32_t i = tid + sid * cols;
  // calculate current raan
  state[i].raan = WrapTwoPI (elem.raan + param[sid].draan * t);
  // calculate current argument of perigee
  state[i].aop = WrapTwoPI (elem.aop + param[sid].daop * t);
  // calculate current mean anomaly
  state[i].ma = WrapTwoPI (elem.ma + param[sid].dma * t);
}

__global__ void CalcTureAnomaly (
  adi::Satellite::Sta*  state,
  adi::Satellite::Par*  param,
  size_t                rows,
  size_t                cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  adi::Satellite::Ele& elem = param[sid].elem;
  uint32_t i = tid + sid * cols;
  double M = state[i].ma;
  double e = elem.ecc;
  if (e == 0)
  {
    state[i].ta = M;
    state[i].dta = param[sid].dma;
    return;
  }
  double M1 = M / (1 - e);
  double M2 = M + e;
  double M3 = M2 * (M_PI - M) / (1 + e);
  double E = M1 > M2 ? (M2) : M1;
  E = E > M3 ? M3 : E;
  double Eold = E + 0.1;
  for (int count = 0; abs (E - Eold) > 1e-10 || count < 10;count++)
  {
    Eold = E;
    double fE = E - M - e * sin (E);
    double dfE = 1 - e * cos (E);
    double ddfE = e * sin (E);
    E = WrapTwoPI (E - fE / (dfE - 0.5 * fE * ddfE / dfE));
  }
  state[i].ta = WrapTwoPI ((2 * atan (sqrt ((1 + e) / (1 - e)) * tan (E / 2))));
  double tmp = 1 - e * cos (E);
  state[i].dta = param[sid].dma * sqrt (1 - e * e) / (tmp * tmp);
}

__global__ void CalcOrbitPosition (
  Vector*               position,
  adi::Satellite::Par*  param,
  adi::Satellite::Sta*  state,
  size_t                rows,
  size_t                cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  // Calculate the orbital radius
  state[i].radius = param[sid].slr / (1 + param[sid].elem.ecc * cos (state[i].ta));
  position[i].x = state[i].radius * cos (state[i].ta);
  position[i].y = state[i].radius * sin (state[i].ta);
  position[i].z = 0.0;
}

__global__ void CalcOrbitVelocity (
  Vector*               velocity,
  adi::Satellite::Par*  param,
  adi::Satellite::Sta*  state,
  size_t                rows,
  size_t                cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  double tmp = sqrt (K_MU / param[sid].slr);
  velocity[i].x = -tmp * sin (state[i].ta);
  velocity[i].y = tmp * (param[sid].elem.ecc + cos (state[i].ta));
  velocity[i].z = 0.0;
}

__global__ void CalcMatrixFromOrbToEci (
  Matrix*               matrix,
  adi::Satellite::Par*  param,
  adi::Satellite::Sta*  state,
  size_t                rows,
  size_t                cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  double cos_raan = cos (state[i].raan);
  double sin_raan = sin (state[i].raan);
  double cos_aop = cos (state[i].aop);
  double sin_aop = sin (state[i].aop);
  double cos_inc = cos (param[sid].elem.inc);
  double sin_inc = sin (param[sid].elem.inc);
  matrix[i].r1.x = cos_raan * cos_aop - sin_raan * cos_inc * sin_aop;
  matrix[i].r1.y = -cos_raan * sin_aop - sin_raan * cos_inc * cos_aop;
  matrix[i].r1.z = sin_raan * sin_inc;
  matrix[i].r2.x = sin_raan * cos_aop + cos_raan * cos_inc * sin_aop;
  matrix[i].r2.y = -sin_raan * sin_aop + cos_raan * cos_inc * cos_aop;
  matrix[i].r2.z = -cos_raan * sin_inc;
  matrix[i].r3.x = sin_inc * sin_aop;
  matrix[i].r3.y = sin_inc * cos_aop;
  matrix[i].r3.z = cos_inc;
}

__global__ void CalcMatrixFromEciToBody (
  Matrix*               matrix,
  adi::Satellite::Par*  param,
  adi::Satellite::Sta*  state,
  size_t                rows,
  size_t                cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  double cos_ta_plus_aop = cos (state[i].ta + state[i].aop);
  double sin_ta_plus_aop = sin (state[i].ta + state[i].aop);
  double cos_raan = cos (state[i].raan);
  double sin_raan = sin (state[i].raan);
  double cos_inc = cos (param[sid].elem.inc);
  double sin_inc = sin (param[sid].elem.inc);
  matrix[i].r1.x = -sin_ta_plus_aop * cos_raan - cos_ta_plus_aop * cos_inc * sin_raan;
  matrix[i].r1.y = -sin_ta_plus_aop * sin_raan + cos_ta_plus_aop * cos_inc * cos_raan;
  matrix[i].r1.z = cos_ta_plus_aop * sin_inc;
  matrix[i].r2.x = -sin_inc * sin_raan;
  matrix[i].r2.y = sin_inc * cos_raan;
  matrix[i].r2.z = -cos_inc;
  matrix[i].r3.x = -cos_ta_plus_aop * cos_raan + sin_ta_plus_aop * cos_inc * sin_raan;
  matrix[i].r3.y = -cos_ta_plus_aop * sin_raan - sin_ta_plus_aop * cos_inc * cos_raan;
  matrix[i].r3.z = -sin_ta_plus_aop * sin_inc;
}


__global__ void CalcDerivMatrixFromEciToBody (
  Matrix*               derivMatrix,
  Matrix*               matrix,
  adi::Satellite::Par*  param,
  adi::Satellite::Sta*  state,
  size_t                rows,
  size_t                cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  double dta_puls_daop = state[i].dta + param[sid].daop;
  double draan = param[sid].draan;
  derivMatrix[i].r1.x = matrix[i].r3.x * dta_puls_daop - matrix[i].r1.y * draan;
  derivMatrix[i].r1.y = matrix[i].r3.y * dta_puls_daop + matrix[i].r1.x * draan;
  derivMatrix[i].r1.z = matrix[i].r3.z * dta_puls_daop;
  derivMatrix[i].r2.x = -matrix[i].r2.y * draan;
  derivMatrix[i].r2.y = matrix[i].r2.x * draan;
  derivMatrix[i].r2.z = 0;
  derivMatrix[i].r3.x = -matrix[i].r1.x * dta_puls_daop - matrix[i].r3.y * draan;
  derivMatrix[i].r3.y = -matrix[i].r1.y * dta_puls_daop + matrix[i].r3.x * draan;
  derivMatrix[i].r3.z = -matrix[i].r1.z * dta_puls_daop;
  // printf ("|%3.8f, %3.8f, %3.8f|\n|%3.8f, %3.8f, %3.8f|\n|%3.8f, %3.8f, %3.8f|\n--------------------------------------------------------------------------------\n",
  //   derivMatrix[i].r1.x, derivMatrix[i].r1.y, derivMatrix[i].r1.z,
  //   derivMatrix[i].r2.x, derivMatrix[i].r2.y, derivMatrix[i].r2.z,
  //   derivMatrix[i].r3.x, derivMatrix[i].r3.y, derivMatrix[i].r3.z);
}
}
}
}