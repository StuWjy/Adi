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

#include <queue>
#include "adi-link.h"
#include "adi-satellite.h"
#include "adi-turntable.h"
#include "adi-link-kernel.h"
#include "adi-sun.h"
#include "adi-util.h"

namespace adi {

__host__ __device__ Vector Body2Left (const Vector& body)
{
  return Vector {-body.y, body.x, body.z};
}

__host__ __device__ Vector Body2Right (const Vector& body)
{
  return Vector {body.y, -body.x, body.z};
}

__host__ __device__ Vector Body2Front (const Vector& body)
{
  return Vector {body.x, body.y, body.z};
}

__host__ __device__ Vector Body2Back (const Vector& body)
{
  return Vector {-body.x, -body.y, body.z};
}

__host__ __device__ Vector Body2Top (const Vector& body)
{
  return Vector {-body.z, -body.y, body.x};
}

__host__ __device__ Vector Body2Bottom (const Vector& body)
{
  return Vector {body.z, body.y, -body.x};
}

// __host__ __device__ Vector StationBody2Top (const Vector& body)
// {
//   return Vector {body.z, body.y, -body.x};
// }

const __device__ Vec2Vec d_body2left = adi::Body2Left;
const __device__ Vec2Vec d_body2right= adi::Body2Right;
const __device__ Vec2Vec d_body2front= adi::Body2Front;
const __device__ Vec2Vec d_body2back   = adi::Body2Back;
const __device__ Vec2Vec d_body2top    = adi::Body2Top;
const __device__ Vec2Vec d_body2bottom = adi::Body2Bottom;

Link::Link ()
: m_epoch         (J2000)
, m_step          (Second)
, m_intervals     (Intervals ())
, m_src           (NULL)
, m_dst           (NULL)
, m_maxDistance   (K_MAX_DISTANCE)
, m_includedState (SRC2DST | DST2SRC)
, m_excludedState (SRC_DAY | DST_DAY)
{
}

Link::Link (
  Turntable*  src,
  Turntable*  dst,
  uint8_t     included,
  uint8_t     excluded
)
: m_epoch         (J2000)
, m_step          (Second)
, m_intervals     (Intervals ())
, m_src           (src)
, m_dst           (dst)
, m_maxDistance   (K_MAX_DISTANCE)
, m_includedState (included)
, m_excludedState (excluded)
{
}

Link::~Link ()
{
}

void
Link::SetIntervalList (const Intervals& intervals)
{
  m_intervals = intervals;;
}

const Intervals&
Link::AddInterval (const Interval& interval)
{
  Interval _interval = interval;
  _interval.InsertToList (m_intervals);
  return m_intervals;
}

const Intervals&
Link::GetIntervalList () const
{
  return m_intervals;
}

void
Link::SetEpoch (const DateTime& epoch)
{
  m_epoch = epoch;
}

const DateTime&
Link::GetEpoch () const
{
  return m_epoch;
}

void
Link::SetStep (const TimeSpan& step)
{
  m_step = step;
}

const TimeSpan&
Link::GetStep () const
{
  return m_step;
}

void
Link::SetMaxDistance (const double& maxDistance)
{
  m_maxDistance = maxDistance;
}

const double&
Link::GetMaxDistance () const
{
  return m_maxDistance;
}

Intervals
Link::CalcLinkInterval ()
{
  int minGridSize, BlockSize;
  size_t N_TIME = Interval::GetTotalTicks (m_intervals, m_step);
  Object* src = m_src->GetObject ();
  Object* dst = m_dst->GetObject ();
  Object::Type srcType = src->GetType ();
  Object::Type dstType = dst->GetType ();
  src->SetEpoch (m_epoch);
  src->SetStep (m_step);
  src->SetIntervalList (m_intervals);
  dst->SetEpoch (m_epoch);
  dst->SetStep (m_step);
  dst->SetIntervalList (m_intervals);
  Trajectory srcTra = src->CalcTrajectory ();
  // {
  //   std::vector<Vector> h_pos = std::vector<Vector> (N_TIME);
  //   cudaMemcpy (&h_pos[0], srcTra.pos, sizeof (adi::Vector) * N_TIME, cudaMemcpyDeviceToHost);
  //   for (const Vector& pos : h_pos)
  //   {
  //     std::cout << pos << std::endl;
  //   }
  //   std::vector<Vector> h_vel = std::vector<Vector> (N_TIME);
  //   cudaMemcpy (&h_vel[0], srcTra.vel, sizeof (adi::Vector) * N_TIME, cudaMemcpyDeviceToHost);
  //   for (const Vector& vel : h_vel)
  //   {
  //     std::cout << vel << std::endl;
  //   }
  // }
  Trajectory dstTra = dst->CalcTrajectory ();
  // {
  //   std::vector<Vector> h_pos = std::vector<Vector> (N_TIME);
  //   cudaMemcpy (&h_pos[0], dstTra.pos, sizeof (adi::Vector) * N_TIME, cudaMemcpyDeviceToHost);
  //   for (const Vector& pos : h_pos)
  //   {
  //     std::cout << pos << std::endl;
  //   }
  //   std::vector<Vector> h_vel = std::vector<Vector> (N_TIME);
  //   cudaMemcpy (&h_vel[0], dstTra.vel, sizeof (adi::Vector) * N_TIME, cudaMemcpyDeviceToHost);
  //   for (const Vector& vel : h_vel)
  //   {
  //     std::cout << vel << std::endl;
  //   }
  // }
  Vector* d_solPos = Sun::CalcPosition (m_intervals, m_epoch, m_step);
  // {
  //   std::vector<Vector> h_pos = std::vector<Vector> (N_TIME);
  //   cudaMemcpy (&h_pos[0], d_solPos, sizeof (adi::Vector) * N_TIME, cudaMemcpyDeviceToHost);
  //   for (const Vector& pos : h_pos)
  //   {
  //     std::cout << pos << std::endl;
  //   }
  // }
  dim3 blockSize, gridSize;
  uint8_t* d_state = NULL;
  cudaMalloc ((void**)&d_state, sizeof (uint8_t) * N_TIME);
  if (srcType == Object::SATELLITE && dstType == Object::SATELLITE)
  {
    Vec2Vec h_vec2vec_1, h_vec2vec_2;
    FaceTransform (h_vec2vec_1, m_src->GetFace ());
    FaceTransform (h_vec2vec_2, m_dst->GetFace ());
    cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::link::CalcLink, 0, N_TIME);
    blockSize = dim3 (BlockSize);
    gridSize = dim3 (1, 1, (N_TIME + blockSize.x - 1) / blockSize.x);
    cuda::link::CalcLink <<< gridSize, blockSize >>> (
      d_state,
      srcTra.pos, srcTra.mat, h_vec2vec_1, m_src->GetPointingBound (),
      dstTra.pos, dstTra.mat, h_vec2vec_2, m_dst->GetPointingBound (),
      d_solPos, m_maxDistance, K_FOV,
      N_TIME);
    cudaCheck;
  }
  else if (srcType == Object::SATELLITE && dstType == Object::STATION)
  {
    cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::link::CalcAccess, 0, N_TIME);
    blockSize = dim3 (BlockSize);
    gridSize = dim3 (1, 1, (N_TIME + blockSize.x - 1) / blockSize.x);
    cuda::link::CalcAccess <<< gridSize, blockSize >>> (
      d_state,
      srcTra.pos, srcTra.mat, m_src->GetPointingBound (),
      dstTra.pos, dstTra.mat, m_dst->GetPointingBound (),
      d_solPos, m_maxDistance, K_FOV,
      1, 1, N_TIME);
    cudaCheck;
  }
  std::vector<uint8_t> h_state = std::vector<uint8_t> (N_TIME);
  cudaMemcpy ((void*)&h_state[0], (void*)d_state, sizeof (uint8_t) * N_TIME, cudaMemcpyDeviceToHost);
  cudaCheck;
  cudaFree (d_state);
  cudaCheck;
  Intervals intervals;
  bool prevState = false;
  bool currState;
  DateTime startTime = m_intervals.front ().GetStart ();
  DateTime stopTime  = m_intervals.back ().GetStop ();
  uint32_t k = 0;
  for (uint32_t i = 0;i < m_intervals.size ();++i)
  {
    uint32_t Ticks = m_intervals[i].GetTicks (m_step);
    for (uint32_t j = 0;j < Ticks;++j)
    {
      currState = IsAcceptedState (h_state[k]);
      if (currState && !prevState)
      {
        startTime = std::max (
          m_intervals[i].GetStart ().AddTicks (j * m_step.Ticks ()),
          m_intervals[i].GetStart ()
        );
      }
      if ((!currState && prevState) || (i == m_intervals.size () - 1 && j == Ticks - 1 && currState))
      {
        stopTime = std::min (
          m_intervals[i].GetStart ().AddTicks (j * m_step.Ticks ()),
          m_intervals[i].GetStop ()
        );
        Interval (startTime, stopTime).InsertToList (intervals);
      }
      prevState = currState;
      k++;
    }
  }
  return intervals;
}

LinkInfoList
Link::CalcLinkData ()
{
  LinkInfoList result;
  int minGridSize, BlockSize;
  size_t N_TIME = Interval::GetTotalTicks (m_intervals, m_step);
  Object* src = m_src->GetObject ();
  Object* dst = m_dst->GetObject ();
  Object::Type srcType = src->GetType ();
  Object::Type dstType = dst->GetType ();
  src->SetEpoch (m_epoch);
  src->SetStep (m_step);
  src->SetIntervalList (m_intervals);
  dst->SetEpoch (m_epoch);
  dst->SetStep (m_step);
  dst->SetIntervalList (m_intervals);
  Matrix** d_srcDMat = new Matrix*;
  Matrix** d_dstDMat = new Matrix*;
  Trajectory srcTra = src->CalcTrajectory (d_srcDMat);
  Trajectory dstTra = dst->CalcTrajectory (d_dstDMat);
  Sun::CalcPosition (m_intervals, m_epoch, m_step);
  dim3 blockSize, gridSize;
  uint8_t* d_state = NULL;
  cudaMalloc ((void**)&d_state, sizeof (uint8_t) * N_TIME);
  cudaCheck;
  Pointing* d_src2dst = NULL;
  Pointing* d_dst2src = NULL;
  cudaMalloc ((void**)&d_src2dst, sizeof (Pointing) * N_TIME);
  cudaCheck;
  cudaMalloc ((void**)&d_dst2src, sizeof (Pointing) * N_TIME);
  cudaCheck;
  double* d_distance = NULL;
  cudaMalloc ((void**)&d_distance, sizeof (double) * N_TIME);
  cudaCheck;
  if (srcType == Object::SATELLITE && dstType == Object::SATELLITE)
  {
    Vec2Vec h_vec2vec_1, h_vec2vec_2;
    FaceTransform (h_vec2vec_1, m_src->GetFace ());
    FaceTransform (h_vec2vec_2, m_dst->GetFace ());
    cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::link::CalcLinkData, 0, N_TIME);
    cudaCheck;
    blockSize = dim3 (BlockSize);
    gridSize = dim3 (1, 1, (N_TIME + blockSize.x - 1) / blockSize.x);
    cuda::link::CalcLinkData <<< gridSize, blockSize >>> (
      d_state,
      srcTra.pos, srcTra.vel, srcTra.mat,
      *d_srcDMat, h_vec2vec_1, m_src->GetPointingBound (),
      d_src2dst,
      dstTra.pos, dstTra.vel, dstTra.mat,
      *d_dstDMat, h_vec2vec_2, m_dst->GetPointingBound (),
      d_dst2src,
      d_distance, Sun::d_pos, m_maxDistance, K_FOV,
      N_TIME);
    cudaCheck;
  }
  else if (srcType == Object::SATELLITE && dstType == Object::STATION)
  {
    cudaOccupancyMaxPotentialBlockSize (&minGridSize, &BlockSize, cuda::link::CalcAccessData, 0, N_TIME);
    cudaCheck;
    blockSize = dim3 (BlockSize);
    gridSize = dim3 (1, 1, (N_TIME + blockSize.x - 1) / blockSize.x);
    cuda::link::CalcAccessData <<< gridSize, blockSize >>> (
      d_state,
      srcTra.pos, srcTra.vel, srcTra.mat,
      *d_srcDMat, m_src->GetPointingBound (),
      d_src2dst,
      dstTra.pos, dstTra.vel, dstTra.mat,
      *d_dstDMat, m_dst->GetPointingBound (),
      d_dst2src,
      d_distance, Sun::d_pos, m_maxDistance, K_FOV,
      1, 1, N_TIME);
    cudaCheck;
  }
  uint8_t* h_state = new uint8_t[N_TIME];
  Pointing* h_src2dst = new Pointing[N_TIME];
  Pointing* h_dst2src = new Pointing[N_TIME];
  double* h_distance = new double[N_TIME];
  cudaMemcpy (h_state, d_state, sizeof (uint8_t) * N_TIME, cudaMemcpyDeviceToHost);
  cudaCheck;
  cudaMemcpy (h_src2dst, d_src2dst, sizeof (Pointing) * N_TIME, cudaMemcpyDeviceToHost);
  cudaCheck;
  cudaMemcpy (h_dst2src, d_dst2src, sizeof (Pointing) * N_TIME, cudaMemcpyDeviceToHost);
  cudaCheck;
  cudaMemcpy (h_distance, d_distance, sizeof (double) * N_TIME, cudaMemcpyDeviceToHost);
  cudaCheck;
  cudaFree (d_state);
  cudaCheck;
  cudaFree (d_src2dst);
  cudaCheck;
  cudaFree (d_dst2src);
  cudaCheck;
  cudaFree (d_distance);
  cudaCheck;
  cudaFree (*d_srcDMat);
  cudaCheck;
  cudaFree (*d_dstDMat);
  cudaCheck;
  LinkDatas linkDatas;
  uint32_t k = 0;
  for (uint32_t i = 0;i < m_intervals.size ();++i)
  {
    uint32_t Ticks = m_intervals[i].GetTicks (m_step);
    for (uint32_t j = 0;j < Ticks;++j)
    {
      if (IsAcceptedState (h_state[k]))
      {
        DateTime t = m_intervals[i].GetStart ().AddTicks (j * m_step.Ticks ());
        linkDatas.push_back (LinkData {h_state[k], t, h_src2dst[k], h_dst2src[k], h_distance[k]});
      }
      if (!IsAcceptedState (h_state[k]) && linkDatas.size () > 0)
      {
        result.push_back (LinkInfo {m_src, m_dst, linkDatas});
        linkDatas.clear ();
      }
      if (j == Ticks - 1)
      {
        if (!linkDatas.empty ())
        {
          result.push_back (LinkInfo {m_src, m_dst, linkDatas});
          linkDatas.clear ();
        }
      }
      k++;
    }
  }
  delete h_state;
  delete h_src2dst;
  delete h_dst2src;
  delete h_distance;
  return result;
}

Satellite*
Link::GetSatellite (size_t i)
{
  if (i == 0)
    return dynamic_cast<Satellite*> (m_src->GetObject ());
  else if (i == 1)
    return dynamic_cast<Satellite*> (m_dst->GetObject ());
  else
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": input a number either 0 or 1";
    return NULL;
  }
}

Turntable*
Link::GetTurntable (size_t i)
{
  if (i == 0)
    return m_src;
  else if (i == 1)
    return m_dst;
  else
    return NULL;
}

Face
Link::GetFace (size_t i) const
{
  if (i == 0)
    return m_src->GetFace ();
  else if (i == 1)
    return m_dst->GetFace ();
  else
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": i < 2";
    return ErrFace;
  }
}

void
Link::FaceTransform (Vec2Vec& host, Face face)
{
  switch (face)
  {
    case Left:
      cudaMemcpyFromSymbol(&host, d_body2left, sizeof(Vec2Vec));
      break;
    case Right:
      cudaMemcpyFromSymbol(&host, d_body2right, sizeof(Vec2Vec));
      break;
    case Front:
      cudaMemcpyFromSymbol(&host, d_body2front, sizeof(Vec2Vec));
      break;
    case Back:
      cudaMemcpyFromSymbol(&host, d_body2back, sizeof(Vec2Vec));
      break;
    case Top:
      cudaMemcpyFromSymbol(&host, d_body2top, sizeof(Vec2Vec));
      break;
    case Bottom:
      cudaMemcpyFromSymbol(&host, d_body2bottom, sizeof(Vec2Vec));
      break;
    default:
      std::cerr << __FILE__ << " " << __LINE__ << ": link only use left, right, front and back face";
      return;
  }
}

bool
Link::IsAcceptedState (const uint8_t& state) const
{
  return ((state & m_includedState) == m_includedState) && ((state & m_excludedState) == 0);
}

bool
Link::operator== (const Link& link)
{
  return (m_src == link.m_src && m_dst == link.m_dst);
}

}