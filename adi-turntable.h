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

#ifndef ADI_TURNTABLE_H
#define ADI_TURNTABLE_H

#include "adi-object.h"
#include "adi-type-define.h"

namespace adi {

class Turntable
{
public:
  Turntable ();
  Turntable (const Turntable& turntable);
  ~Turntable ();
  void SetFace (const Face& face);
  const Face& GetFace () const;
  void SetObject (Object* obj);
  Object* GetObject () const;
  void SetPointingBound (const PointingBound& bound);
  const PointingBound& GetPointingBound () const;
  uint32_t GetId () const;
  void SetNs3 (void* ns3);
  void* GetNs3 (void) const;
private:
  void*   m_ns3;
  Face    m_face;
  Object* m_obj;
  PointingBound m_bound;
  uint32_t m_uid;
};

}

#endif /* ADI_TURNTABLE_H */