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

#include "adi-station-container.h"

namespace adi {
size_t
StationContainer::GetN () const
{
  return m_container.size ();
}

StationContainer
StationContainer::Create (size_t n)
{
  StationContainer c;
  Station* stations = new Station[n];
  for (uint32_t i = 0;i < n;i++)
  {
    stations[i].SetElement (Station::Ele {0.0, 0.0, 0.0});
    c << stations + i;
  }
  return c;
}

Station*
StationContainer::operator[] (size_t i)
{
  if (i >= m_container.size ())
  {
    return NULL;
  }
  return m_container[i];
}

StationContainer
StationContainer::operator() (size_t i)
{
  StationContainer c;
  c << m_container[i];
  return c;
}

StationContainer&
StationContainer::operator<< (Station* sta)
{
  m_container.push_back (sta);
  return *this;
}

StationContainer&
StationContainer::operator<< (Station& sta)
{
  m_container.push_back (&sta);
  return *this;
}

void
StationContainer::Install (uint32_t faces, const PointingBound& bound)
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
StationContainer::Install (size_t i, const Face& face, const PointingBound& bound)
{
  m_container[i]->Install (face, bound);
}

}