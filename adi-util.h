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

#ifndef ADI_UTIL_H
#define ADI_UTIL_H

#include <cuda_runtime.h>
#include "adi-type-define.h"

namespace adi {
#define cudaCheck \
{ \
  cudaError_t cudaStatus; \
  cudaStatus = cudaGetLastError (); \
  if (cudaSuccess != cudaStatus)  \
  { \
    std::cout << __FILE__ << ": " << __LINE__ << " Error code: " << cudaGetErrorString (cudaStatus) << std::endl; \
  } \
}

extern __host__ __device__ double Wrap(const double x, const double y);

extern __host__ __device__ double WrapTwoPI (double x);

extern __host__ __device__ double Wrap360 (double x);

extern __host__ __device__ double DegreesToRadians (double x);
extern __host__ __device__ double RadiansToDegrees (double x);

extern __host__ __device__ double ToJulian (int64_t ticks);

extern __host__ __device__ double ToJ2000(int64_t ticks);

extern __host__ __device__ Vector Add (Vector v1, Vector v2);

extern __host__ __device__ Vector Sub (Vector v1, Vector v2);

extern __host__ __device__ void Scale (Vector& v, double scale);

extern __host__ __device__ double Dot (Vector v1, Vector v2);

extern __host__ __device__ Vector operator* (Matrix& m, Vector& v);

extern __host__ __device__ Vector operator+ (const Vector& v1, const Vector& v2);
extern __host__ __device__ Vector operator- (const Vector& v1, const Vector& v2);
extern __host__ __device__ Vector operator- (const Vector& v);

extern __host__ __device__ Vector Cross (Vector v1, Vector v2);

extern __host__ __device__ void Normalize (Vector& v);

extern __host__ __device__ double Norm (const Vector& v);

extern __host__ __device__ Direction View (Vector v);

extern __host__ __device__ Pointing View (const Vector& pos, const Vector& vel);

// extern __host__ __device__ Pointing ViewRate (const Vector& pos, const Vector& vel);

extern __host__ __device__ bool IsInside (const Direction& d, const PointingBound& b);

extern __host__ __device__ void Swap (double& a, double& b);

// template<class T>
// T* MemcpyDeviceToHost (T* ptrDevice, size_t bytes);

// template<class T>
// void PrintDevice (T* ptrDevice, size_t bytes);

// Vectors MemcpyDeviceToHost (Vectors device);
// Trajectory MemcpyDeviceToHost (Trajectory device);
// Trajectories MemcpyDeviceToHost (Trajectories device);

namespace cuda {

/**
 * \brief Calculate the multiplication of square matrix to vector b[i] = A[i] * b[i]
 * \param[in]     A     the square matrix A
 * \param[inout]  b     the vector b
 * \param[in]     rows  the rows of threads
 * \param[in]     cols  the columns of threads
 */
__global__ void CalcMatsMulVecs (
  Matrix* A,
  Vector* b,
  size_t  rows,
  size_t  cols);

/**
 * \brief Calculate the multiplication of square matrix to vector c[i] = A[i] * b[i]
 * \param[out]  c     the vector c
 * \param[in]   A     the square matrix A
 * \param[in]   b     the vector b
 * \param[in]   rows  the rows of threads
 * \param[in]   cols  the columns of threads
 */
__global__ void CalcVecsEqMatsMulVecs (
  Vector* c,
  Matrix* A,
  Vector* b,
  size_t  rows,
  size_t  cols);

/**
 * \brief Calculate the multiplication of square matrix to vector c[i] = A * b[i]
 * \param[out]  c     the vector c
 * \param[in]   A     the square matrix A
 * \param[in]   b     the vector b
 * \param[in]   rows  the rows of threads
 * \param[in]   cols  the columns of threads
 */
__global__ void CalcVecsEqMatMulVecs (
  Vector* c,
  Matrix  A,
  Vector* b,
  size_t  rows,
  size_t  cols);

/**
 * \brief Subtraction of vectors,c = a - b
 * \param[out]  c     the vector c
 * \param[in]   a     the vector a
 * \param[in]   b     the vector b
 * \param[in]   rows  the rows of threads
 * \param[in]   cols  the columns of threads
 */
__global__ void CalcVectorSubtraction (
  Vector* c,
  Vector* a,
  Vector* b,
  size_t  rows,
  size_t  cols);

/**
 * \brief Normalization of vectors,a = a / |a|
 * \param[in]   a     the vector a
 * \param[in]   b     the vector b
 * \param[in]   rows  the rows of threads
 * \param[in]   cols  the columns of threads
 */
__global__ void CalcVectorNormalization (
  Vector* a,
  size_t  rows,
  size_t  cols);

/**
 * \brief Calculate a and b, c = a and b
 * \param[out]  c     the 1-dim array of bool
 * \param[in]   a     the 1-dim array of bool
 * \param[in]   b     the 1-dim array of bool
 * \param[in]   b     the vector b
 * \param[in]   rows  the rows of threads
 * \param[in]   cols  the columns of threads
 */
__global__ void CalcAnd (
  bool*   c,
  bool*   a,
  bool*   b,
  size_t  rows,
  size_t  cols);
} // namespace cuda
} // namespace adi
#endif /* UTIL_H */