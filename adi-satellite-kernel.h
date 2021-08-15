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

#ifndef ADI_SATELLITE_KERNEL_H
#define ADI_SATELLITE_KERNEL_H

#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "adi-satellite.h"

namespace adi {
namespace cuda {
namespace satellite {
  /**
   * \brief Calculate the orbital parameters of satellites
   * \param[out] param    the 1-dim array of satellite's orbital parameters,
   *                      the parameter of satellite i is param[i]
   * \param[in] rows      the rows of param, i.e., the number of satellites
   */
  __global__ void CalcParam(
    adi::Satellite::Par*  param,
    size_t                rows);

  /**
   * \brief Calculate the orbital states of satellite
   * \param[out] state    the 1-dim array of satellite's orbital states,
   *                      the state of satellite i at time j is state[i * cols + j]
   * \param[in] param     the 1-dim vector of satellite's orbital parameters
   * \param[in] t0        the initial elapsed time
   * \param[in] dt        the time step
   * \param[in] rows      the rows of state, i.e., the number of satellites
   * \param[in] cols      the columns of state, i.e., the number of times
   */
  __global__ void CalcState (
    adi::Satellite::Sta*  state,
    adi::Satellite::Par*  param,
    double*               times,
    size_t                rows,
    size_t                cols);

    /**
   * \brief Calculate the ture anomaly
   * \param[out] state    the 1-dim array of satellite's orbital states
   * \param[in] param     the 1-dim vector of satellite's orbital parameters
   * \param[in] rows      the rows of state, i.e., the number of satellites
   * \param[in] cols      the columns of state, i.e., the number of times
   */
  __global__ void CalcTureAnomaly (
    adi::Satellite::Sta*  state,
    adi::Satellite::Par*  param,
    size_t                rows,
    size_t                cols);

  /**
   * \brief Calculate the position in orbital coordinate
   * \param[out] position     the 1-dim array of the position in orbital coordinate,
   *                          the position of satellite i at time j is position[i * cols + j]
   * \param[in] param         the 1-dim array of satellite's orbital parameters
   * \param[out] state        the 1-dim array of satellite's orbital states,
   * \param[in] rows          the rows of state, i.e., the number of satellites
   * \param[in] cols          the columns of state, i.e., the number of times
   */
  __global__ void CalcOrbitPosition (
    Vector*               position,
    adi::Satellite::Par*  param,
    adi::Satellite::Sta*  state,
    size_t                rows,
    size_t                cols);
    

  /**
   * \brief Calculate the velocity in orbital coordinate
   * \param[out] velocity     the 1-dim array of the velocity in orbital coordinate,
   *                          the velocity of satellite i at time j is position[i * cols + j]
   * \param[in] param         the 1-dim array of satellite's orbital parameters
   * \param[in] state         the 1-dim array of satellite's orbital states,
   * \param[in] rows          the rows of state, i.e., the number of satellites
   * \param[in] cols          the columns of state, i.e., the number of times
   */
  __global__ void CalcOrbitVelocity (
    Vector*               velocity,
    adi::Satellite::Par*  param,
    adi::Satellite::Sta*  state,
    size_t                rows,
    size_t                cols);

  /**
   * \brief Calculate the transform matrix from orbital to eci coordinate
   * \param[out] matrix       the 1-dim array of transform matrix,
   *                          the transform matrix of satellite i at time j is matrix[i * cols + j]
   * \param[in] param         the 1-dim array of satellite's orbital parameters
   * \param[in] state         the 1-dim array of satellite's orbital states
   * \param[in] rows          the rows of state, i.e., the number of satellites
   * \param[in] cols          the columns of state, i.e., the number of times
   */
  __global__ void CalcMatrixFromOrbToEci (
    Matrix*               matrix,
    adi::Satellite::Par*  param,
    adi::Satellite::Sta*  state,
    size_t                rows,
    size_t                cols);

  /**
   * \brief Calculate the transformation matrix from eci to body coordinate system
   * \param[out] matrix       the 1-dim array of transform matrix
   * \param[in] param         the 1-dim array of satellite's orbital parameters
   * \param[in] state         the 1-dim array of satellite's orbital states
   * \param[in] rows          the rows of state, i.e., the number of satellites
   * \param[in] cols          the columns of state, i.e., the number of times
   */
  __global__ void CalcMatrixFromEciToBody (
    Matrix*               matrix,
    adi::Satellite::Par*  param,
    adi::Satellite::Sta*  state,
    size_t                rows,
    size_t                cols);

  /**
   * \brief Calculate the derivative of transformation matrix from eci to body coordinate system
   * \param[out]  derivMatrix the 1-dim array of derivative of transform matrix
   * \param[in]   matrix      the 1-dim array of transform matrix
   * \param[in]   param       the 1-dim array of satellite's orbital parameters
   * \param[in]   state       the 1-dim array of satellite's orbital states
   * \param[in]   rows        the rows of state, i.e., the number of satellites
   * \param[in]   cols        the columns of state, i.e., the number of times
   */
  __global__ void CalcDerivMatrixFromEciToBody (
    Matrix*               derivMatrix,
    Matrix*               matrix,
    adi::Satellite::Par*  param,
    adi::Satellite::Sta*  state,
    size_t                rows,
    size_t                cols);
} // namespace satellite
} // namespace cuda
} // namespace adi
#endif /* CUDA_SATELLITE_KERNEL_H */