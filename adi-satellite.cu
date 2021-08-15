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
#include "adi-satellite.h"
#include "adi-satellite-kernel.h"
#include "adi-satellite-list.h"

extern const double K_RE;

namespace adi {

std::string
Satellite::Element::ToString () const
{
  std::stringstream ss;
  ss.setf (std::ios::fixed);
  ss.precision (4);
  ss.width (10);
  ss << "                  semi-major axis: " << std::setw (10) << sma                    << std::endl
     << "                     eccentricity: " << std::setw (10) << ecc                    << std::endl
     << "                      inclination: " << std::setw (10) << RadiansToDegrees (inc) << std::endl
     << "right ascension of ascending node: " << std::setw (10) << RadiansToDegrees (raan)<< std::endl
     << "              argument of perigee: " << std::setw (10) << RadiansToDegrees (aop) << std::endl
     << "                     mean anomaly: " << std::setw (10) << RadiansToDegrees (ma)  << std::endl;
  return ss.str ();
}

bool
Satellite::Ele::operator== (Satellite::Ele& ele) const
{
  return sma == ele.sma && ecc == ele.ecc && inc == ele.inc
    && raan == ele.raan && aop == ele.aop && ma == ele.ma;
}

std::string
Satellite::Parameter::ToString () const
{
  std::stringstream ss;
  ss.setf (std::ios::fixed);
  ss.precision (4);
  ss.width (10);
  ss << elem.ToString ();
  ss << "                     angular rate: " << std::setw (10) << ar                       << std::endl
     << "                semi-latus rectum: " << std::setw (10) << slr                      << std::endl
     << "rate of                            "                                               << std::endl
     << "right ascension of ascending node: " << std::setw (10) << RadiansToDegrees (draan) << std::endl
     << "              argument of perigee: " << std::setw (10) << RadiansToDegrees (daop)  << std::endl
     << "                     mean anomaly: " << std::setw (10) << RadiansToDegrees (dma)   << std::endl;
  return ss.str ();
}

    
std::string
Satellite::State::ToString () const
{
  std::stringstream ss;
  ss.setf (std::ios::fixed);
  ss.precision (4);
  ss.width (10);
  ss << "                           current radius: " << std::setw (10) << RadiansToDegrees (radius)<< std::endl
     << "current right ascension of ascending node: " << std::setw (10) << RadiansToDegrees (raan)  << std::endl
     << "              current argument of perigee: " << std::setw (10) << RadiansToDegrees (aop)   << std::endl
     << "                     current mean anomaly: " << std::setw (10) << RadiansToDegrees (ma)    << std::endl
     << "                     current true anomaly: " << std::setw (10) << RadiansToDegrees (ta)    << std::endl;
  return ss.str ();
}

Satellite::Satellite ()
: Object  ()
, m_par   (Par ())
, d_par   (NULL)
, d_sta   (NULL)
{
}

Satellite::Satellite (const Satellite& sat)
: Object  (sat)
, m_par   (sat.m_par)
, d_par   (sat.d_par)
, d_sta   (sat.d_sta)
{
}

Satellite::~Satellite ()
{
}

void
Satellite::SetElement (Ele ele)
{
  m_par.elem = ele;
  Construct ();
}

void
Satellite::SetElement (
  double sma,
  double ecc,
  double inc,
  double raan,
  double aop,
  double ma,
  bool isDegree)
{
  if (isDegree)
  {
    inc = DegreesToRadians (inc);
    raan = DegreesToRadians (raan);
    aop = DegreesToRadians (aop);
    ma = DegreesToRadians (ma);
  }
  inc = WrapTwoPI (inc);
  raan = WrapTwoPI (raan);
  aop = WrapTwoPI (aop);
  ma = WrapTwoPI (ma);
  m_par.elem = Ele {sma, ecc, inc, raan, aop, ma};
  Construct ();
}

std::string
Satellite::GetName () const
{
  return m_name;
}

size_t
Satellite::GetId () const
{
  return m_uid;
}

Object::Type
Satellite::GetType () const
{
  return SATELLITE;
}

Trajectory
Satellite::CalcTrajectory (Matrix** dmat)
{
  Initialzie ();
  size_t N_TIME = GetIntervalTicks ();
  std::vector<double> h_times = Interval::CreateSeconds (m_intervals, m_epoch, m_step);
  double* d_times = NULL;
  cudaMalloc ((void**)&d_times, sizeof (double) * N_TIME);
  cudaMemcpy (d_times, &h_times[0], sizeof (double) * N_TIME, cudaMemcpyHostToDevice);
  int minGridSize, BlockSize;
  dim3 gridSize, blockSize;
  // calculate states of satellite
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::satellite::CalcState, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  // std::cout << blockSize << std::endl;
  // std::cout << gridSize << std::endl;
  cuda::satellite::CalcState <<< gridSize, blockSize >>> (d_sta, d_par, d_times, 1, N_TIME);
  cudaCheck;
  cudaFree (d_times);
  cudaCheck;
  // calculate true anomaly
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::satellite::CalcTureAnomaly, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::satellite::CalcTureAnomaly <<< gridSize, blockSize >>> (d_sta, d_par, 1, N_TIME);
  cudaCheck;
  // {
  //   size_t bytes = sizeof (Satellite::Sta) * N_TIME;
  //   Satellite::Sta* h_sta = new Satellite::Sta[N_TIME];
  //   cudaMemcpy (h_sta, d_sta, bytes, cudaMemcpyDeviceToHost);
  //   for (size_t i = 0;i < N_TIME;i++)
  //     std::cout << h_sta[i] << std::endl;
  // }
  // calculate position in orbital coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::satellite::CalcOrbitPosition, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::satellite::CalcOrbitPosition <<< gridSize, blockSize >>> (d_pos, d_par, d_sta, 1, N_TIME);
  cudaCheck;
  // calculate velocity in orbital coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::satellite::CalcOrbitVelocity, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::satellite::CalcOrbitVelocity <<< gridSize, blockSize >>> (d_vel, d_par, d_sta, 1, N_TIME);
  cudaCheck;
  // calculate transform matrix from orbital to Earth-centered-inertial coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::satellite::CalcMatrixFromOrbToEci, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::satellite::CalcMatrixFromOrbToEci <<< gridSize, blockSize >>> (d_mat, d_par, d_sta, 1, N_TIME);
  cudaCheck;
  // calculate position in Earth-centered-inertial coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::CalcMatsMulVecs, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::CalcMatsMulVecs <<< gridSize, blockSize >>> (d_mat, d_pos, 1, N_TIME);
  cudaCheck;
  // {
  //   size_t bytes = sizeof (Vector) * N_TIME;
  //   Vector* h_pos = new Vector[N_TIME];
  //   cudaMemcpy (h_pos, d_pos, bytes, cudaMemcpyDeviceToHost);
  //   for (size_t i = 0;i < N_TIME;i++)
  //     std::cout << h_pos[i] << std::endl;
  // }
  // calculate velocity in Earth-centered-inertial coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::CalcMatsMulVecs, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::CalcMatsMulVecs <<< gridSize, blockSize >>> (d_mat, d_vel, 1, N_TIME);
  cudaCheck;
  // {
  //   size_t bytes = sizeof (Vector) * N_TIME;
  //   Vector* h_vel = new Vector[N_TIME];
  //   cudaMemcpy (h_vel, d_vel, bytes, cudaMemcpyDeviceToHost);
  //   for (size_t i = 0;i < N_TIME;i++)
  //     std::cout << h_vel[i] << std::endl;
  // }
  // calculate the transform matrix from Earth-centered-inertial to body coordinate system
  cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::satellite::CalcMatrixFromEciToBody, 0, N_TIME);
  blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
  gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
  cuda::satellite::CalcMatrixFromEciToBody <<< gridSize, blockSize >>> (d_mat, d_par, d_sta, 1, N_TIME);
  cudaCheck;
  if (dmat != NULL)
  {
    Matrix* d_dmat;
    cudaMalloc ((void**)&d_dmat, sizeof (Matrix) * N_TIME);
    cudaCheck;
    // calculate the derivative transform matrix from Earth-centered-inertial to body coordinate system if need
    cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::satellite::CalcDerivMatrixFromEciToBody, 0, N_TIME);
    blockSize = dim3 (std::min (N_TIME, (size_t)BlockSize));
    gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);
    cuda::satellite::CalcDerivMatrixFromEciToBody <<< gridSize, blockSize >>> (d_dmat, d_mat, d_par, d_sta, 1, N_TIME);
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
Satellite::GetTurntable (Face face)
{
  if (m_turntables.find (face) == m_turntables.end ())
  {
    return NULL;
  }
  return m_turntables[face];
}

void
Satellite::Construct ()
{
  if (m_uid == 0)
    m_uid = SatelliteList::Add (this);
  cudaFree (d_par);
  cudaCheck;
  cudaMalloc ((void**)&d_par, sizeof (Satellite::Par));
  cudaCheck;
  // memory copy
  cudaMemcpy ((void*)d_par, (void*)&m_par, sizeof (Satellite::Par), cudaMemcpyHostToDevice);
  cudaCheck;
  // calculate the parameters of satellite
  dim3 blockSize = dim3 (1);
  dim3 gridSize = dim3 (1);
  cuda::satellite::CalcParam <<< gridSize, blockSize >>> (d_par, 1);
  cudaCheck;
  m_name = "Sat-" + std::to_string (m_uid);
}

void
Satellite::Initialzie ()
{
  Release ();
  size_t N_TIME = GetIntervalTicks ();
  // allocating memory 
  cudaMalloc ((void**)&d_sta, sizeof (Sta)    * N_TIME);
  cudaMalloc ((void**)&d_pos, sizeof (Vector) * N_TIME);
  cudaMalloc ((void**)&d_vel, sizeof (Vector) * N_TIME);
  cudaMalloc ((void**)&d_mat, sizeof (Matrix) * N_TIME);
}

void
Satellite::Release ()
{
  cudaFree (d_sta);
  cudaCheck;
  cudaFree (d_pos);
  cudaCheck;
  cudaFree (d_vel);
  cudaCheck;
  cudaFree (d_mat);
  cudaCheck;
}

std::ostream& operator<< (std::ostream& os, const Satellite::Ele& ele)
{
  os << ele.ToString ();
  return os;
}

std::ostream& operator<< (std::ostream& os, const Satellite::Par& par)
{
  os << par.ToString ();
  return os;
}

std::ostream& operator<< (std::ostream& os, const Satellite::Sta& sta)
{
  os << sta.ToString ();
  return os;
}

}