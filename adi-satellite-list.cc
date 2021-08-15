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

#include "adi-satellite-list.h"

namespace adi {

std::map<uint32_t, Satellite*>  SatelliteList::m_satellites = std::map<uint32_t, Satellite*> ();
size_t                          SatelliteList::m_count      = 0;

uint32_t
SatelliteList::Add (Satellite* sat)
{
  m_count++;
  m_satellites[m_count] = sat;
  // std::cout << "A new satellite added, whose id is " << m_count << std::endl;
  return m_count;
}

Satellite*
SatelliteList::Find (uint32_t id)
{
  return m_satellites[id];
}

uint32_t
SatelliteList::GetN ()
{
  return m_satellites.size ();
}

SatelliteContainer
SatelliteList::GetAll ()
{
  SatelliteContainer c;
  for (uint32_t i = 0;i < m_count;++i)
  {
    c << m_satellites[i];
  }
  return c;
}

// void
// SatelliteList::InstallAll ()
// {
//   for (uint32_t i = 0;i < m_satellites.size ();i++)
//   {
//     m_satellites[i]->InstallAllFace ();
//   }
// }

}