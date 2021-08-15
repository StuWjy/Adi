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

#include "adi-turntable-list.h"

namespace adi {

std::map<uint32_t, Turntable*>  TurntableList::m_turntables = std::map<uint32_t, Turntable*> ();
uint32_t                        TurntableList::m_count      = 0;

uint32_t
TurntableList::Add (Turntable* t)
{
  m_count++;
  m_turntables[m_count] = t;
  // std::cout << "A new turntable installed, whose id is " << m_count << std::endl;
  return m_count;
}

Turntable*
TurntableList::Find (uint32_t id)
{
  return m_turntables[id];
}

uint32_t
TurntableList::GetN ()
{
  return m_turntables.size ();
}

}