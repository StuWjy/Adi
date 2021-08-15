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

#include "adi-station-kernel.h"
#include "adi-constant.h"
#include "adi-util.h"

namespace adi {
namespace cuda {
namespace station {

__device__ double ToGreenwichSiderealTime(int64_t ticks)
{
  // julian date of previous midnight
  double jd0 = floor(ToJulian(ticks) + 0.5) - 0.5;
  // julian centuries since epoch
  double t   = (jd0 - 2451545.0) / 36525.0;
  double jdf = ToJulian(ticks) - jd0;

  double gt  = 24110.54841 + t * (8640184.812866 + t * (0.093104 - t * 6.2E-6));
  gt  += jdf * 1.00273790935 * 86400.0;
  double gst = WrapTwoPI ((gt / 240.0) * M_PI / 180.0);
  return gst; 
}

__global__ void CalcLocalGreenwichSiderealTime (
  Gst*                gst,
  adi::Station::Ele*  element,
  int64_t*            ticks,
  size_t              rows,
  size_t              cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  gst[i] = WrapTwoPI (ToGreenwichSiderealTime (ticks[tid]) + element[sid].lon);
}

__global__ void CalcEciPosition (
  Vector*             position,
  adi::Station::Ele*  element,
  Gst*                gst,
  size_t              rows,
  size_t              cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  double sin_theta = sin (gst[i]);
  double cos_theta = cos (gst[i]);
  double sin_lat = sin (element[sid].lat);
  double cos_lat = cos (element[sid].lat);
  double c = 1 / sqrt (1.0 + K_F * (K_F - 2.0) * sin_lat * sin_lat);
  double s = (1 - K_F) * (1 - K_F) * c;
  double achcp = (K_RE * c + element[sid].alt) * cos_lat;
  position[i].x = achcp * cos_theta;
  position[i].y = achcp * sin_theta;
  position[i].z = (K_RE * s + element[sid].alt) * sin_lat;
}


__global__ void CalcEciVelocity (
  Vector*       velocity,
  Vector*       position,
  size_t        rows,
  size_t        cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = tid + sid * cols;
  velocity[i].x = -K_EAR * position[i].y;
  velocity[i].y =  K_EAR * position[i].x;
  velocity[i].z = 0.0;
}

__global__ void CalcMatrixFromEciToBody (
  Matrix*             matrix,
  adi::Station::Ele*  element,
  Gst*                gst,
  size_t              rows,
  size_t              cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  double sin_lat = sin (element[sid].lat);
  double cos_lat = cos (element[sid].lat);
  double sin_theta = sin (gst[tid]);
  double cos_theta = cos (gst[tid]);
  uint32_t i = sid * cols + tid;
  matrix[i].r1.x = sin_lat * cos_theta;
  matrix[i].r1.y = sin_lat * sin_theta;
  matrix[i].r1.z = -cos_lat;
  matrix[i].r2.x = -sin_theta;
  matrix[i].r2.y = cos_theta;
  matrix[i].r2.z = 0.0;
  matrix[i].r3.x = cos_lat * cos_theta;
  matrix[i].r3.y = cos_lat * sin_theta;
  matrix[i].r3.z = sin_lat;
}

__global__ void CalcDerivMatrixFromEciToBody (
  Matrix*             derivMatrix,
  Matrix*             matrix,
  adi::Station::Ele*  element,
  Gst*                gst,
  size_t              rows,
  size_t              cols)
{
  uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t sid = blockIdx.y;
  if (tid >= cols || sid >= rows)
  {
    return;
  }
  uint32_t i = sid * cols + tid;
  derivMatrix[i].r1.x = -matrix[i].r1.y * K_EAR;
  derivMatrix[i].r1.y = matrix[i].r1.x * K_EAR;
  derivMatrix[i].r1.z = 0;
  derivMatrix[i].r2.x = -matrix[i].r2.y * K_EAR;
  derivMatrix[i].r2.y = matrix[i].r2.x * K_EAR;
  derivMatrix[i].r2.z = 0.0;
  derivMatrix[i].r3.x = -matrix[i].r3.y * K_EAR;
  derivMatrix[i].r3.y = matrix[i].r3.x * K_EAR;
  derivMatrix[i].r3.z = 0;
}

}
}
}