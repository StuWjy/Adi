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

#include "adi-util.h"
#include "adi-time-span.h"

namespace adi {
__host__ __device__ double Wrap(const double x, const double y)
{
  if (y == 0.0)
  {
      return x;
  }
  return x - y * floor(x / y);
}

__host__ __device__ double WrapTwoPI (double x)
{
  return Wrap (x, M_PI * 2.0);
}

__host__ __device__ double Wrap360 (double x)
{
  return Wrap (x, 360.0);
}

__host__ __device__ double DegreesToRadians (double x)
{
  return x * M_PI / 180.0;
}

__host__ __device__ double RadiansToDegrees (double x)
{
  return x * 180.0 / M_PI;
}

__host__ __device__ double ToJulian (int64_t ticks)
{
  return static_cast<double>(ticks) / TicksPerDay + 1721425.5;
}

__host__ __device__ double ToJ2000(int64_t ticks)
{
    return ToJulian(ticks) - 2415020.0;
}

__host__ __device__ Vector Add (Vector v1, Vector v2)
{
  Vector v;
  v.x = v1.x + v2.x;
  v.y = v1.y + v2.y;
  v.z = v1.z + v2.z;
  return v;
}

__host__ __device__ Vector Sub (Vector v1, Vector v2)
{
  Vector v;
  v.x = v1.x - v2.x;
  v.y = v1.y - v2.y;
  v.z = v1.z - v2.z;
  return v;
}

__host__ __device__ void Scale (Vector& v, double scale)
{
  v.x *= scale;
  v.y *= scale;
  v.z *= scale;
}

__host__ __device__ double Dot (Vector v1, Vector v2)
{
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

__host__ __device__ Vector operator* (Matrix& m, Vector& v)
{
  Vector r;
  r.x = m.r1.x * v.x + m.r1.y * v.y + m.r1.z * v.z;
  r.y = m.r2.x * v.x + m.r2.y * v.y + m.r2.z * v.z;
  r.z = m.r3.x * v.x + m.r3.y * v.y + m.r3.z * v.z;
  return r;
}

__host__ __device__ Vector operator+ (const Vector& v1, const Vector& v2)
{
  Vector r;
  r.x = v1.x + v2.x;
  r.y = v1.y + v2.y;
  r.z = v1.z + v2.z;
  return r;
}

__host__ __device__ Vector operator- (const Vector& v1, const Vector& v2)
{
  Vector r;
  r.x = v1.x - v2.x;
  r.y = v1.y - v2.y;
  r.z = v1.z - v2.z;
  return r;
}

__host__ __device__ Vector operator- (const Vector& v)
{
  Vector r;
  r.x = -v.x;
  r.y = -v.y;
  r.z = -v.z;
  return r;
}

__host__ __device__ Vector Cross (Vector v1, Vector v2)
{
  Vector v;
  v.x = v1.y * v2.z - v1.z * v2.y;
  v.y = v1.z * v2.x - v1.x * v2.z;
  v.z = v1.x * v2.y - v1.y * v2.x;
  return v;
}

__host__ __device__ void Normalize (Vector& v)
{
  double r = Norm (v);
  v.x /= r;
  v.y /= r;
  v.z /= r;
}

__host__ __device__ double Norm (const Vector& v)
{
  return sqrt (v.x * v.x + v.y * v.y + v.z * v.z);
}

__host__ __device__ Direction View (Vector v)
{
  Direction direction;
  direction.azimuth = atan2 (v.y, v.x);
  direction.pitch = -atan (v.z / sqrt (v.x * v.x + v.y * v.y));
  return direction;
}

__host__ __device__ Pointing View (const Vector& pos, const Vector& vel)
{
  Pointing pointing;
  pointing.angle.azimuth = atan2 (pos.y, pos.x);
  pointing.angle.pitch = -atan (pos.z / sqrt (pos.x * pos.x + pos.y * pos.y));
  pointing.rate.azimuth = (pos.x * vel.y - pos.y * vel.x) / (pos.x * pos.x + pos.y * pos.y);
  double r = Norm (pos);
  double d2 = pos.x * pos.x + pos.y * pos.y;
  pointing.rate.pitch =  (pos.z * (pos.x * vel.x + pos.y * vel.y) - d2 * vel.z) / (r * r * sqrt (d2));
  return pointing;
}

// __host__ __device__ Pointing ViewRate (const Vector& pos, const Vector& vel)
// {
//   Pointing rate;
//   rate.azimuth = (pos.x * vel.y - pos.y * vel.x) / (pos.x * pos.x + pos.y * pos.y);
//   double r = Norm (pos);
//   double d2 = pos.x * pos.x + pos.y * pos.y;
//   rate.pitch =  (pos.z * (pos.x * vel.x + pos.y * vel.y) - d2 * vel.z) / (r * r * sqrt (d2));
//   return rate;
// }

__host__ __device__ bool IsInside (const Direction& d, const PointingBound& b)
{
  if (d.azimuth > b.azimuth.max || d.azimuth < b.azimuth.min)
  {
    return false;
  }
  if (d.pitch > b.pitch.max || d.pitch < b.pitch.min)
  {
    return false;
  }
  return true;
}

__host__ __device__ void Swap (double& a, double& b)
{
  a = a + b;
  b = a - b;
  a = a - b;
}

// template<class T>
// T* MemcpyDeviceToHost (T* ptrDevice, size_t bytes)
// {
//   T* ptrHost = new T[bytes];
//   cudaMemcpy (ptrHost, ptrDevice, bytes, cudaMemcpyDeviceToHost);
//   return ptrHost;
// }

// template<class T>
// void PrintDevice (T* ptrDevice, size_t bytes)
// {
//   T* ptrHost = MemcpyDeviceToHost (ptrDevice, bytes);
//   for (size_t i = 0;i < bytes / sizeof (T);i++)
//     std::cout << ptrHost[i] << std::endl;
// }

// Vectors MemcpyDeviceToHost (Vectors device)
// {
//   Vectors host;
//   host.ptr = (Vector*)malloc (sizeof (Vector) * device.num);
//   host.num = device.num;
//   cudaMemcpy (host.ptr, device.ptr, sizeof (Vector) * device.num, cudaMemcpyDeviceToHost);
//   return host;
// }

// Trajectory MemcpyDeviceToHost (Trajectory device)
// {
//   Trajectory host;
//   host.position.ptr = (Vector*)malloc (sizeof (Vector) * device.position.num);
//   host.position.num = device.position.num;
//   host.velocity.ptr = (Vector*)malloc (sizeof (Vector) * device.velocity.num);
//   host.velocity.num = device.velocity.num;
//   cudaMemcpy (host.position.ptr, device.position.ptr, sizeof (Vector) * device.position.num, cudaMemcpyDeviceToHost);
//   cudaMemcpy (host.velocity.ptr, device.velocity.ptr, sizeof (Vector) * device.velocity.num, cudaMemcpyDeviceToHost);
//   return host;
// }

// Trajectories MemcpyDeviceToHost (Trajectories device)
// {
//   Trajectories host;
//   for (size_t i = 0;i < device.size ();i++)
//   {
//     host.push_back (MemcpyDeviceToHost (device[i]));
//   }
//   return host;
// }

namespace cuda {

__global__ void CalcMatsMulVecs (
  Matrix* A,
  Vector* b,
  size_t  rows,
  size_t  cols)
{
  uint32_t j = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t i = threadIdx.y + blockIdx.y * blockDim.y;
  if (j >= cols || i >= rows)
  {
    return;
  }
  uint32_t k = j + i * cols;
  double x = b[k].x;
  double y = b[k].y;
  double z = b[k].z;
  b[k].x = A[k].r1.x * x + A[k].r1.y * y + A[k].r1.z * z;
  b[k].y = A[k].r2.x * x + A[k].r2.y * y + A[k].r2.z * z;
  b[k].z = A[k].r3.x * x + A[k].r3.y * y + A[k].r3.z * z;
}

__global__ void CalcVecsEqMatsMulVecs (
  Vector* c,
  Matrix* A,
  Vector* b,
  size_t  rows,
  size_t  cols)
{
  uint32_t j = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t i = threadIdx.y + blockIdx.y * blockDim.y;
  if (j >= cols || i >= rows)
  {
    return;
  }
  uint32_t k = j + i * cols;
  c[k].x = A[k].r1.x * b[k].x + A[k].r1.y * b[k].y + A[k].r1.z * b[k].z;
  c[k].y = A[k].r2.x * b[k].x + A[k].r2.y * b[k].y + A[k].r2.z * b[k].z;
  c[k].z = A[k].r3.x * b[k].x + A[k].r3.y * b[k].y + A[k].r3.z * b[k].z;
}

__global__ void CalcVecsEqMatMulVecs (
  Vector* c,
  Matrix  A,
  Vector* b,
  size_t  rows,
  size_t  cols)
{
  uint32_t j = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t i = threadIdx.y + blockIdx.y * blockDim.y;
  if (j >= cols || i >= rows)
  {
    return;
  }
  uint32_t k = j + i * cols;
  c[k].x = A.r1.x * b[i].x + A.r1.y * b[i].y + A.r1.z * b[i].z;
  c[i].y = A.r2.x * b[i].x + A.r2.y * b[i].y + A.r2.z * b[i].z;
  c[k].z = A.r3.x * b[i].x + A.r3.y * b[i].y + A.r3.z * b[i].z;
}

__global__ void CalcVectorSubtraction (
  Vector* c,
  Vector* a,
  Vector* b,
  size_t  rows,
  size_t  cols)
{
  uint32_t j = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t i = threadIdx.y + blockIdx.y * blockDim.y;
  if (j >= cols || i >= rows)
  {
    return;
  }
  uint32_t k = j + i * cols;
  c[k].x = a[k].x - b[k].x;
  c[k].y = a[k].y - b[k].y;
  c[k].z = a[k].z - b[k].z;
}

__global__ void CalcVectorNormalization (
  Vector* a,
  size_t  rows,
  size_t  cols)
{
  uint32_t j = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t i = threadIdx.y + blockIdx.y * blockDim.y;
  if (j >= cols || i >= rows)
  {
    return;
  }
  uint32_t k = j + i * cols;
  double r = a[k].x * a[k].x + a[k].y * a[k].y + a[k].z * a[k].z;
  r = sqrt (r);
  a[k].x = a[k].x / r;
  a[k].y = a[k].y / r;
  a[k].z = a[k].z / r;
}

__global__ void CalcAnd (
  bool* c,
  bool* a,
  bool* b,
  size_t  rows,
  size_t  cols)
{
  uint32_t j = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t i = threadIdx.y + blockIdx.y * blockDim.y;
  if (j >= cols || i >= rows)
  {
    return;
  }
  uint32_t k = j + i * cols;
  c[k] = a[k] && b[k];
}
}
}