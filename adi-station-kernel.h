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

#ifndef ADI_STATION_KERNEL_H
#define ADI_STATION_KERNEL_H

#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "adi-station.h"

namespace adi {
namespace cuda {
namespace station {

// /**
//  * \brief Calculate the loacl Greenwich sidereal time of given stations and time
//  * \param[out]  gst       the 1-dim array of local Greenwich sidereal time
//  * \param[in]   element   the 1-dim array of station's parameter
//  * \param[in]   t0_ticks  the ticks of initial time
//  * \param[in]   dt_ticks  the ticks of time step
//  * \param[in]   rows      the rows of threads
//  * \param[in]   cols      the columns of threads
//  */
// __global__ void CalcLocalGreenwichSiderealTime (
//   Gst*                gst,
//   adi::Station::Ele*  element,
//   int64_t             t0_ticks,
//   int64_t             dt_ticks,
//   size_t              rows,
//   size_t              cols);

/**
 * \brief Calculate the loacl Greenwich sidereal time of given stations and time
 * \param[in]   element   the 1-dim array of station's element
 * \param[in]   ticks     the 1-dim array of ticks of times
 * \param[in]   rows      the rows of threads
 * \param[in]   cols      the columns of threads
 */
__global__ void CalcLocalGreenwichSiderealTime (
  Gst*                gst,
  adi::Station::Ele*  element,
  int64_t*            ticks,
  size_t              rows,
  size_t              cols);

/**
 * \brief Calculate the earth-centered-inertial position of stations
 * \param[out]  position  the 1-dim array of eci position
 * \param[in]   element   the 1-dim array of station's element
 * \param[in]   gst       the 1-dim array of local Greenwich Sidereal time
 * \param[in]   rows      the rows of threads
 * \param[in]   cols      the columns of threads
 */
__global__ void CalcEciPosition (
  Vector*             position,
  adi::Station::Ele*  element,
  Gst*                gst,
  size_t              rows,
  size_t              cols);

/**
 * \brief Calculate the velocity of stations in Earth-centered-inertial coordinate system
 * \param[out]  velocity  the 1-dim array of velocity in eci coordinate system
 * \param[in]   position  the 1-dim array of position in eci coordinate system
 * \param[in]   rows      the rows of threads
 * \param[in]   cols      the columns of threads
 */
__global__ void CalcEciVelocity (
  Vector*       velocity,
  Vector*       position,
  size_t        rows,
  size_t        cols);

/**
 * \brief Calculate the transformation matrix from eci to body coordinate system
 * \param[out]  matrix    the transformation matrix
 * \param[in]   element   the 1-dim array of station's element
 * \param[in]   gst       the 1-dim array of local Greenwich Sidereal time
 * \param[in]   rows      the rows of threads
 * \param[in]   cols      the columns of threads
 */
__global__ void CalcMatrixFromEciToBody (
  Matrix*             matrix,
  adi::Station::Ele*  element,
  Gst*                gst,
  size_t              rows,
  size_t              cols);

/**
 * \brief Calculate the derivative of transformation matrix from eci to body coordinate system
 * \param[out]  derivMatrix the derivative of transformation matrix
 * \param[out]  matrix      the transformation matrix
 * \param[in]   element     the 1-dim array of station's element
 * \param[in]   gst         the 1-dim array of local Greenwich Sidereal time
 * \param[in]   rows        the rows of threads
 * \param[in]   cols        the columns of threads
 */
__global__ void CalcDerivMatrixFromEciToBody (
  Matrix*             derivMatrix,
  Matrix*             matrix,
  adi::Station::Ele*  element,
  Gst*                gst,
  size_t              rows,
  size_t              cols);

} // namespace station
} // namespace cuda
} // namespace adi
#endif /* STATION_KERNEL_H */