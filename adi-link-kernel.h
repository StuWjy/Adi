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

#ifndef ADI_LINK_KERNEL_H
#define ADI_LINK_KERNEL_H

#include "adi-type-define.h"
#include "adi-link.h"

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
  size_t        depth);

__global__ void CalcLinkData(
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
  size_t        depth);


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
  size_t        depth);

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
  size_t        depth);

} // namespace link
} // namespace cuda
} // namespace adi

#endif /* ADI_LINK_KERNEL_H */