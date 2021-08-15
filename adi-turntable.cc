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

#include "adi-turntable.h"
#include "adi-turntable-list.h"
#include "adi-util.h"
#include "adi-constant.h"

namespace adi {

Turntable::Turntable ()
: m_ns3   (NULL)
, m_face  (ErrFace)
, m_obj   (NULL)
, m_bound (K_FOV)
{
  m_uid = TurntableList::Add (this);
}

Turntable::Turntable (const Turntable& turntable)
: m_ns3   (turntable.m_ns3)
, m_face  (turntable.m_face)
, m_obj   (turntable.m_obj)
, m_bound (turntable.m_bound)
{
}

Turntable::~Turntable ()
{
}

void
Turntable::SetFace (const Face& face)
{
  m_face = face;
}

const Face&
Turntable::GetFace () const
{
  return m_face;
}

void
Turntable::SetObject (Object* object)
{
  m_obj = object;
}

Object*
Turntable::GetObject () const
{
  return m_obj;
}

void
Turntable::SetPointingBound (const PointingBound& bound)
{
  m_bound = bound;
  if (m_bound.azimuth.max > M_PI_2 || m_bound.azimuth.max < -M_PI_2)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": the azimuth boundary is incorrect";
  }
  if (m_bound.pitch.max > M_PI_2 || m_bound.pitch.max < -M_PI_2)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": the pitch boundary is incorrect";
  }
}

const PointingBound&
Turntable::GetPointingBound () const
{
  return m_bound;
}

uint32_t
Turntable::GetId () const
{
  return m_uid;
}

void
Turntable::SetNs3 (void* ns3)
{
  m_ns3 = ns3;
}

void*
Turntable::GetNs3 (void) const
{
  return m_ns3;
}

}