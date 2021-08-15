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

#include <algorithm>
#include "adi-satellite-container.h"

namespace adi {

size_t
SatelliteContainer::GetN () const
{
  return m_container.size ();
}

Satellite*
SatelliteContainer::operator[] (size_t i)
{
  if (i >= m_container.size ())
  {
    return NULL;
  }
  return m_container[i];
}

SatelliteContainer
SatelliteContainer::operator() (size_t i)
{
  SatelliteContainer c;
  c << m_container[i];
  return c;
}

SatelliteContainer&
SatelliteContainer::operator<< (Satellite* sat)
{
  m_container.push_back (sat);
  return *this;
}

SatelliteContainer&
SatelliteContainer::operator<< (Satellite& sat)
{
  m_container.push_back (&sat);
  return *this;
}

SatelliteContainer
SatelliteContainer::DistributeEvenly (adi::Satellite::Ele ele, size_t N)
{
  SatelliteContainer container;
  Satellite* sats = new Satellite[N] ;
  for (uint32_t i = 0;i < N;i++)
  {
    sats[i].SetElement (ele);
    container << sats[i];
    ele.ma += M_PI * 2.0 / N;
  }
  return container;
}

void
SatelliteContainer::Install (uint32_t faces, const PointingBound& bound)
{
  for (size_t i = 0;i < m_container.size ();i++)
  {
    if ((faces & Top) != 0)
    {
      Install (i, Top, bound);
    }
    if ((faces & Bottom) != 0)
    {
      Install (i, Bottom, bound);
    }
    if ((faces & Left) != 0)
    {
      Install (i, Left, bound);
    }
    if ((faces & Right) != 0)
    {
      Install (i, Right, bound);
    }
    if ((faces & Front) != 0)
    {
      Install (i, Front, bound);
    }
    if ((faces & Back) != 0)
    {
      Install (i, Back, bound);
    }
  }
}

void
SatelliteContainer::Install (size_t i, const Face& face, const PointingBound& bound)
{
  m_container[i]->Install (face, bound);
}

}