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

#ifndef ADI_TURNTABLE_LIST_H
#define ADI_TURNTABLE_LIST_H

#include <map>
#include "adi-turntable.h"

namespace adi {

class TurntableList
{
public:
  static uint32_t Add (Turntable* t);
  static Turntable* Find (uint32_t id);
  static uint32_t GetN ();
private:
  static std::map<uint32_t, Turntable*> m_turntables;
  static uint32_t                       m_count;
};

}

#endif /* ADI_TURNTABLE_LIST_H */