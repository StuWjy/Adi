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

#include "adi-station-list.h"

namespace adi {

std::map<uint32_t, Station*> StationList::m_stations = std::map<uint32_t, Station*> ();
uint32_t                     StationList::m_count    = 0;

uint32_t
StationList::Add (Station* sta)
{
  m_count++;
  m_stations[m_count] = sta;
  // std::cout << "A new station added, whose id is " << m_count << std::endl;
  return m_count;
}

Station*
StationList::Find (uint32_t id)
{
  return m_stations[id];
}

uint32_t
StationList::GetN ()
{
  return m_stations.size ();
}

StationContainer
StationList::GetAll ()
{
  StationContainer c;
  for (uint32_t i = 0;i < m_count;++i)
  {
    c << m_stations[i];
  }
  return c;
}

// void
// StationList::InstallAll ()
// {
//   for (uint32_t i = 0;i < m_stations.size ();i++)
//   {
//     m_stations[i]->InstallAllFace ();
//   }
// }

}