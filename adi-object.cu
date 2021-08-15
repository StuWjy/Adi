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

#include "adi-object.h"
#include "adi-turntable.h"
#include "adi-util.h"
#include "adi-constant.h"

namespace adi {

Object::Object ()
: m_ns3         (NULL)
, m_epoch       (J2000)
, m_step        (Second)
, m_intervals   (Intervals ())
, m_uid         (0)
, d_pos         (NULL)
, d_vel         (NULL)
, d_mat         (NULL)
{
}

Object::Object (const Object& obj)
: m_ns3         (obj.m_ns3)
, m_epoch       (obj.m_epoch)
, m_step        (obj.m_step)
, m_intervals   (obj.m_intervals)
, m_uid         (obj.m_uid)
, d_pos         (obj.d_pos)
, d_vel         (obj.d_vel)
, d_mat         (obj.d_mat)
, m_turntables  (obj.m_turntables)
{
}

Object::~Object ()
{
  cudaFree (d_pos);
  cudaCheck;
  cudaFree (d_vel);
  cudaCheck;
  cudaFree (d_mat);
  cudaCheck;
}

void
Object::SetEpoch (const DateTime& epoch)
{
  m_epoch = epoch;
}

const DateTime&
Object::GetEpoch () const
{
  return m_epoch;
}

void
Object::SetStep (const TimeSpan& step)
{
  m_step = step;
}

const TimeSpan&
Object::GetStep () const
{
  return m_step;
}

void
Object::SetIntervalList (const Intervals& intervals)
{
  m_intervals = intervals;
}

const Intervals&
Object::AddInterval (const Interval& interval)
{
  interval.InsertToList (m_intervals);
  return m_intervals;
}

const Intervals&
Object::GetIntervalList () const
{
  return m_intervals;
}

Turntable*
Object::Install (const Face& face, const PointingBound& bound)
{
  if (m_turntables.find (face) == m_turntables.end ())
  {
    Turntable* turntable = new Turntable;
    turntable->SetFace (face);
    turntable->SetObject (this);
    turntable->SetPointingBound (bound);
    m_turntables[face] = turntable;
    return turntable;
  }
  return NULL;
}

Turntable*
Object::GetTurntable (const Face& face)
{
  return m_turntables[face];
}

void
Object::SetNs3 (void* ns3)
{
  m_ns3 = ns3;
}

void*
Object::GetNs3 () const
{
  return m_ns3;
}

size_t
Object::GetIntervalTicks () const
{
  size_t N_TIME = 0;
  for (const Interval& interval : m_intervals)
  {
    N_TIME += interval.GetTicks (m_step);
  }
  return N_TIME;
}

}