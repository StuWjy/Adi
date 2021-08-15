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

#include "adi-link-kernel.h"
#include "adi-util.h"
#include "adi-constant.h"

namespace adi {
namespace cuda {
namespace link {

__global__ void CalcLink(
  uint8_t*      state,
  Vector*       srcPos,
  Matrix*       srcMat,
  Vec2Vec       srcOp,
  PointingBound srcBound,
  Vector*       dstPos,
  Matrix*       dstMat,
  Vec2Vec       dstOp,
  PointingBound dstBound,
  Vector*       sunPos,
  double        maxDistance,
  PointingBound FOV,
  size_t        depth)
{
  uint32_t tid = threadIdx.x + blockIdx.z * blockDim.x;
  if (tid >= depth)
  {
    return;
  }
  state[tid] = 0;
  {
    // calculate whether satellite 1 is in daylight
    Vector _pos = sunPos[tid];
    Normalize (_pos);
    _pos = Cross (_pos, Cross (srcPos[tid], _pos));
    Normalize (_pos);
    double sin_lat = srcPos[tid].z / Norm (srcPos[tid]);
    double c = K_RE / sqrt (1.0 + K_F * (2.0 - K_F) * sin_lat * sin_lat);
    if (Dot (srcPos[tid], _pos) > c || Dot (srcPos[tid], sunPos[tid]) > 0)
    {
      state[tid] |= SRC_DAY;
    }
  }
  {
    // calculate whether satellite 2 is in daylight
    Vector _pos = sunPos[tid];
    Normalize (_pos);
    _pos = Cross (_pos, Cross (dstPos[tid], _pos));
    Normalize (_pos);
    double sin_lat = dstPos[tid].z / Norm (dstPos[tid]);
    double c = K_RE / sqrt (1.0 + K_F * (2.0 - K_F) * sin_lat * sin_lat);
    if (Dot (dstPos[tid], _pos) > c || Dot (dstPos[tid], sunPos[tid]) > 0)
    {
      state[tid] |= DST_DAY;
    }
  }
  bool inFov = false;
  {
    // calculate pointing from satellite 1 to satellite 2
    Vector pos = Sub (dstPos[tid], srcPos[tid]);
    if (Norm (pos) > maxDistance)
    {
      state[tid] |= BEYOND_DISTANCE;
    }
    Vector _pos = srcMat[tid] * pos;
    _pos = srcOp (_pos);
    Direction view = View (_pos);
    if (IsInside (view, FOV))
    {
      inFov = true;
    }
    if (IsInside (view, srcBound))
    {
      state[tid] |= SRC2DST;
    }
  }
  {
    // calculate pointing from satellite 2 to satellite 1
    Vector pos = Sub (srcPos[tid], dstPos[tid]);
    Vector _pos = dstMat[tid] * pos;
    _pos = dstOp (_pos);
    Direction view = View (_pos);
    if (IsInside (view, FOV) && inFov)
    {
      state[tid] |= IN_FOV;
    }
    if (IsInside (view, dstBound))
    {
      state[tid] |= DST2SRC;
    }
  }
}

__global__ void CalcLinkData (
  uint8_t*      state,
  Vector*       srcPos,
  Vector*       srcVel,
  Matrix*       srcMat,
  Matrix*       srcDMat,
  Vec2Vec       srcOp,
  PointingBound srcBound,
  Pointing*     srcView,
  Vector*       dstPos,
  Vector*       dstVel,
  Matrix*       dstMat,
  Matrix*       dstDMat,
  Vec2Vec       dstOp,
  PointingBound dstBound,
  Pointing*     dstView,
  double*       range,
  Vector*       sunPos,
  double        maxDistance,
  PointingBound FOV,
  size_t        depth)
{
  uint32_t tid = threadIdx.x + blockIdx.z * blockDim.x;
  if (tid >= depth)
  {
    return;
  }
  state[tid] = 0;
  // calculate pointing from satellite 1 to satellite 2
  Vector pos = Sub (dstPos[tid], srcPos[tid]);
  if (Norm (pos) > maxDistance)
  {
    state[tid] |= BEYOND_DISTANCE;
  }
  Vector vel = Sub (dstVel[tid], srcVel[tid]);
  Vector _pos = srcMat[tid] * pos;
  Vector _vel = srcDMat[tid] * pos + srcMat[tid] * vel;
  _pos = srcOp (_pos);
  _vel = srcOp (_vel);
  range[tid] = Norm (_pos);
  srcView[tid] = View (_pos, _vel);
  pos = -pos;
  vel = -vel;
  _pos = dstMat[tid] * pos;
  _vel = dstDMat[tid] * pos + dstMat[tid] * vel;
  _pos = dstOp (_pos);
  _vel = dstOp (_vel);
  dstView[tid] = View (_pos, _vel);
  if (IsInside (srcView[tid].angle, FOV) && IsInside (dstView[tid].angle, FOV))
  {
    state[tid] |= IN_FOV;
  }
  if (IsInside (srcView[tid].angle, srcBound))
  {
    state[tid] |= SRC2DST;
  }
  if (IsInside (dstView[tid].angle, dstBound))
  {
    state[tid] |= DST2SRC;
  }
  {
    // calculate whether satellite 1 is in daylight
    _pos = sunPos[tid];
    Normalize (_pos);
    _pos = Cross (_pos, Cross (srcPos[tid], _pos));
    Normalize (_pos);
    double sin_lat = srcPos[tid].z / Norm (srcPos[tid]);
    double c = K_RE / sqrt (1.0 + K_F * (2.0 - K_F) * sin_lat * sin_lat);
    if (Dot (srcPos[tid], _pos) > c || Dot (srcPos[tid], sunPos[tid]) > 0)
    {
      state[tid] |= SRC_DAY;
    }
  }
  {
    // calculate whether satellite 2 is in daylight
    _pos = sunPos[tid];
    Normalize (_pos);
    _pos = Cross (_pos, Cross (dstPos[tid], _pos));
    Normalize (_pos);
    double sin_lat = dstPos[tid].z / Norm (dstPos[tid]);
    double c = K_RE / sqrt (1.0 + K_F * (2.0 - K_F) * sin_lat * sin_lat);
    if (Dot (dstPos[tid], _pos) > c || Dot (dstPos[tid], sunPos[tid]) > 0)
    {
      state[tid] |= DST_DAY;
    }
  }
}


__global__ void CalcAccess (
  uint8_t*      state,
  Vector*       posSat,
  Matrix*       matSat,
  PointingBound boundSat,
  Vector*       posSta,
  Matrix*       matSta,
  PointingBound boundSta,
  Vector*       posSun,
  double        maxDistance,
  PointingBound FOV,
  size_t        rows,
  size_t        cols,
  size_t        depth)
{
  uint32_t tid = threadIdx.x + blockIdx.z * blockDim.x;
  if (tid >= depth || blockIdx.y >= rows || blockIdx.x >= cols)
  {
    return;
  }
  uint32_t i = tid + blockIdx.y * depth;
  uint32_t j = tid + blockIdx.x * depth;
  uint32_t k = tid + blockIdx.x * depth + blockIdx.y * depth * cols;
  double tmp = 0;
  state[k] = 0;
  // calculate the pointing from station to satellite
  Vector vec = Sub (posSat[i], posSta[j]);
  // if the distance exceeds the maximum distance
  if (Norm (vec) > maxDistance)
  {
    state[k] |= BEYOND_DISTANCE;
  }
  // satellite is in the back of station
  #if 0
  if (Dot (vec, posSta[j]) < 0)
  {
    state[k] |= DST_OUT_OF_VIEW;
  }
  #endif
  vec = matSta[j] * vec;
  // convert to the bottom side coordinate system
  // x-axis points to zenith
  // z-axis points to due north
  tmp = vec.x;
  vec.x = vec.z;
  vec.z = -tmp;
  Direction view = View (vec);
  if (IsInside (view, FOV))
  {
    state[k] |= IN_FOV;
  }
  if (IsInside (view, boundSta))
  {
    state[k] |= DST2SRC;
  }
  vec = Sub (posSta[j], posSat[i]);
  vec = matSat[i] * vec;
  // convert to the bottom side coordinate system
  // x-axis points to earth's center point
  // z-axis points to the opposite direction of advance
  tmp = vec.x;
  vec.x = vec.z;
  vec.z = -tmp;
  view = View (vec);
  if (IsInside (view, boundSat))
  {
    state[k] |= SRC2DST;
  }
  if (Dot (posSun[tid], posSta[j]) >= 0)
  {
    state[k] |= DST_DAY;
  }
  vec = posSun[tid];
  Normalize (vec);
  vec = Cross (vec, Cross (posSat[i], vec));
  Normalize (vec);
  // Considering the oblateness of the earth
  double sin_lat = posSat[i].z;
  // Approximate formula
  double c = K_RE / sqrt (1.0 + K_F * (K_F - 2.0) / (1 - K_F) * (1 - K_F) * sin_lat * sin_lat);
  if (Dot (posSat[i], vec) > c)
  {
    state[k] |= SRC_DAY;
  }
}

__global__ void CalcAccessData (
  uint8_t*      state,
  Vector*       posSat,
  Vector*       velSat,
  Matrix*       matSat,
  Matrix*       dMatSat,
  PointingBound boundSat,
  Pointing*     sat2sta,
  Vector*       posSta,
  Vector*       velSta,
  Matrix*       matSta,
  Matrix*       dMatSta,
  PointingBound boundSta,
  Pointing*     sta2sat,
  double*       range,
  Vector*       posSun,
  double        maxDistance,
  PointingBound FOV,
  size_t        rows,
  size_t        cols,
  size_t        depth)
{
  uint32_t tid = threadIdx.x + blockIdx.z * blockDim.x;
  if (tid >= depth || blockIdx.y >= rows || blockIdx.x >= cols)
  {
    return;
  }
  uint32_t i = tid + blockIdx.y * depth;
  uint32_t j = tid + blockIdx.x * depth;
  uint32_t k = tid + blockIdx.x * depth + blockIdx.y * depth * cols;
  double tmp = 0;
  state[k] = 0;
  // calculate pointing from station to satellite
  Vector pos = Sub (posSat[i], posSta[j]);
  Vector vel = Sub (velSat[i], velSta[j]);
  if (Norm (pos) > maxDistance)
  {
    state[k] |= BEYOND_DISTANCE;
  }
  // satellite is in the back of station
  #if 0
  if (Dot (pos, posSta[j]) < 0)
  {
    state[k] |= DST_OUT_OF_VIEW;
  }
  #endif
  Vector _pos = matSta[j] * pos;
  Vector _vel = dMatSta[j] * pos + matSta[j] * vel;
  // convert to the bottom side coordinate system
  // x-axis points to zenith
  // z-axis points to due north
  tmp = _pos.x;
  _pos.x = _pos.z;
  _pos.z = -tmp;
  tmp = _vel.x;
  _vel.x = _vel.z;
  _vel.z = -tmp;
  range[k] = Norm (_pos);
  sta2sat[k] = View (_pos, _vel);
  if (IsInside (sta2sat[k].angle, FOV))
  {
    state[k] |= IN_FOV;
  }
  if (IsInside (sta2sat[k].angle, boundSta))
  {
    state[k] |= DST2SRC;
  }
  pos = -pos;
  vel = -vel;
  _pos = matSat[i] * pos;
  _vel = dMatSat[i] * pos + matSat[i] * vel;
  // convert to the bottom side coordinate system
  // x-axis points to earth's center point
  // z-axis points to the opposite direction of advance
  tmp = _pos.x;
  _pos.x = _pos.z;
  _pos.z = -tmp;
  tmp = _vel.x;
  _vel.x = _vel.z;
  _vel.z = -tmp;
  sat2sta[k] = View (_pos, _vel);
  // printf (
  //   "(%f, %f, %f), (%f, %f, %f), (%f, %f), (%f, %f)\n"
  //   "|%3.8f, %3.8f, %3.8f|\n|%3.8f, %3.8f, %3.8f|\n|%3.8f, %3.8f, %3.8f|\n"
  //   "--------------------------------------------------------------------------------\n",
  //   _pos.x, _pos.y, _pos.z,
  //   _vel.x, _vel.y, _vel.z,
  //   sat2sta[k].azimuth * 180 / M_PI, sat2sta[k].pitch * 180 / M_PI,
  //   sat2staRate[k].azimuth * 180 / M_PI, sat2staRate[k].pitch * 180 / M_PI,
  //   dMatSat[i].r1.x, dMatSat[i].r1.y, dMatSat[i].r1.z,
  //   dMatSat[i].r2.x, dMatSat[i].r2.y, dMatSat[i].r2.z,
  //   dMatSat[i].r3.x, dMatSat[i].r3.y, dMatSat[i].r3.z);
  if (IsInside (sat2sta[k].angle, boundSat))
  {
    state[k] |= SRC2DST;
  }
  if (Dot (posSun[tid], posSta[j]) >= 0)
  {
    state[k] |= DST_DAY;
  }
  _pos = posSun[tid];
  Normalize (_pos);
  _pos = Cross (_pos, Cross (posSat[i], _pos));
  Normalize (_pos);
  // Considering the oblateness of the earth
  double sin_lat = posSat[i].z / Norm (posSat[i]);
  // Approximate formula
  double c = K_RE / sqrt (1.0 + K_F * (2.0 - K_F) * sin_lat * sin_lat);
  if (Dot (posSat[i], _pos) > c)
  {
    state[k] |= SRC_DAY;
  }
}

} // namespace link
} // namespace cuda
} // namespace adi