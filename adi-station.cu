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
#include "adi-station.h"
#include "adi-station-kernel.h"
#include "adi-station-list.h"

namespace adi {
std::string
Station::Element::ToString () const
{
  std::stringstream ss;
  ss.setf (std::ios::fixed);
  ss.precision (4);
  ss.width (10);
  ss << " altitude: " << alt                    << std::endl
      << " latitude: " << RadiansToDegrees (lat) << std::endl
      << "longitude: " << RadiansToDegrees (lon) << std::endl;
  return ss.str ();
}

Station::Station ()
: Object  ()
// , m_uid   (0)
, m_ele   (Ele ())
, d_ele   (NULL)
, d_gst   (NULL)
{
}

Station::Station (const Station& sta)
: Object  (sta)
// , m_uid   (sta.m_uid)
, m_ele   (sta.m_ele)
, d_ele   (sta.d_ele)
, d_gst   (sta.d_gst)
{
}

Station::~Station ()
{
}

void
Station::SetElement (Ele ele)
{
  m_ele = ele;
  Construct ();
}

void
Station::SetElement (
  double lat,
  double lon,
  double alt,
  std::string name,
  bool isDegree)
{
  if (isDegree)
  {
    lat = DegreesToRadians (lat);
    lon = DegreesToRadians (lon);
  }
  lat = WrapTwoPI (lat);
  lon = WrapTwoPI (lon);
  m_ele = Ele {alt, lat, lon};
  Construct ();
  m_name = name;
}

// void
// Station::SetElement (
//   double lat,
//   double lon,
//   double alt,
//   std::string name,
//   bool isDegree)
// {
//   if (isDegree)
//   {
//     lat = DegreesToRadians (lat);
//     lon = DegreesToRadians (lon);
//   }
//   lat = WrapTwoPI (lat);
//   lon = WrapTwoPI (lon);
//   m_ele = Ele {alt, lat, lon};
//   Construct ();
//   m_name = name;
// }

std::string
Station::GetName () const
{
  return m_name;
}

size_t
Station::GetId () const
{
  return m_uid;;
}

Object::Type
Station::GetType () const
{
  return STATION;
}

Trajectory
Station::CalcTrajectory (Matrix** dmat)
{
  Initialize ();
  size_t N_TIME = GetIntervalTicks ();
  std::vector<int64_t> h_times = Interval::CreateTicks (m_intervals, m_epoch, m_step);
  int64_t* d_ticks = NULL;
  cudaMalloc ((void**)&d_ticks, sizeof (int64_t) * N_TIME);
  cudaMemcpy (d_ticks, &h_times[0], sizeof (int64_t) * N_TIME, cudaMemcpyHostToDevice);
  int minGridSize, BlockSize;
  dim3 gridSize, blockSize;
  // calculate the Greenwich Sidereal time
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::station::CalcLocalGreenwichSiderealTime, 0, N_TIME);
  cudaCheck;
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  // std::cout << d_gst << std::endl;
  // std::cout << d_ele << std::endl;
  // std::cout << blockSize << std::endl;
  // std::cout << gridSize << std::endl;
  cuda::station::CalcLocalGreenwichSiderealTime <<< gridSize, blockSize >>> (d_gst, d_ele, d_ticks, 1, N_TIME);
  cudaCheck;
  cudaFree (d_ticks);
  cudaCheck;
  // calculate position in Earth-centered-inertial coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::station::CalcEciPosition, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::station::CalcEciPosition <<< gridSize, blockSize >>> (d_pos, d_ele, d_gst, 1, N_TIME);
  cudaCheck;
  // {
  //   size_t bytes = sizeof (Vector) * N_TIME;
  //   Vector* h_pos = new Vector[N_TIME];
  //   cudaMemcpy (h_pos, d_pos, bytes, cudaMemcpyDeviceToHost);
  //   for (size_t i = 0;i < N_TIME;i++)
  //     std::cout << h_pos[i] << std::endl;
  // }
  // calculate velocity in Earth-centered-inertial coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::station::CalcEciVelocity, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::station::CalcEciVelocity <<< gridSize, blockSize >>> (d_vel, d_pos, 1, N_TIME);
  // cudaCheck;
  // {
  //   size_t bytes = sizeof (Vector) * N_TIME;
  //   Vector* h_vel = new Vector[N_TIME];
  //   cudaMemcpy (h_vel, d_vel, bytes, cudaMemcpyDeviceToHost);
  //   for (size_t i = 0;i < N_TIME;i++)
  //     std::cout << h_vel[i] << std::endl;
  // }
  // calculate the transform matrix from Earth-centered-inertial to body coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::station::CalcMatrixFromEciToBody, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::station::CalcMatrixFromEciToBody <<< gridSize, blockSize >>> (d_mat, d_ele, d_gst, 1, N_TIME);
  cudaCheck;
  if (dmat != NULL)
  {
    Matrix* d_dmat;
    cudaMalloc ((void**)&d_dmat, sizeof (Matrix) * N_TIME);
    // calculate the derivative of transform matrix from Earth-centered-inertial to body coordinate system
    cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::station::CalcDerivMatrixFromEciToBody, 0, N_TIME);
    blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
    gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
    cuda::station::CalcDerivMatrixFromEciToBody <<< gridSize, blockSize >>> (d_dmat, d_mat, d_ele, d_gst, 1, N_TIME);
    cudaCheck;
    *dmat = d_dmat;
  }
  Trajectory tra;
  tra.pos = d_pos;
  tra.vel = d_vel;
  tra.mat = d_mat;
  tra.num = N_TIME;
  return tra;
}

Turntable*
Station::GetTurntable (Face face)
{
  if (face == Top)
    return m_turntables[face];
  else
    return NULL;
}

void
Station::Construct ()
{
  if (m_uid == 0)
    m_uid = StationList::Add (this);
  cudaFree (d_ele);
  cudaCheck;
  cudaMalloc ((void**)&d_ele, sizeof (Station::Ele));
  // memory copy
  cudaMemcpy ((void*)d_ele, (void*)&m_ele, sizeof (Station::Ele), cudaMemcpyHostToDevice);
  cudaCheck;
}

void
Station::Initialize ()
{
  Release ();
  size_t N_TIME = GetIntervalTicks ();
  // allocating memory 
  cudaMalloc ((void**)&d_gst, sizeof (Gst) * N_TIME);
  cudaCheck;
  cudaMalloc ((void**)&d_pos, sizeof (Vector) * N_TIME);
  cudaCheck;
  cudaMalloc ((void**)&d_vel, sizeof (Vector) * N_TIME);
  cudaCheck;
  cudaMalloc ((void**)&d_mat, sizeof (Matrix) * N_TIME);
  cudaCheck;
}

void
Station::Release ()
{
  cudaFree (d_gst);
  cudaCheck;
  cudaFree (d_pos);
  cudaCheck;
  cudaFree (d_vel);
  cudaCheck;
  cudaFree (d_mat);
  cudaCheck;
}

std::ostream& operator<< (std::ostream& os, Station::Ele& ele)
{
  os << ele.ToString ();
  return os;
}

}