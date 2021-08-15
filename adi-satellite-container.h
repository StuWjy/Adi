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

#ifndef ADI_SATELLITE_CONTAINER_H
#define ADI_SATELLITE_CONTAINER_H

#include <vector>
#include "adi-satellite.h"

namespace adi {

class SatelliteContainer
{
public:
  /**
   * \brief Default constructor
   */
  SatelliteContainer (){m_container.clear ();}

  /**
   * \brief Copy constructor
   */
  SatelliteContainer (const SatelliteContainer& c)
  {
    m_container.insert (m_container.begin (), c.m_container.begin (), c.m_container.end ());
  }

  /**
   * \brief Deconstructor
   */
  ~SatelliteContainer (){}
  size_t GetN () const;
  Satellite* operator[] (size_t i);
  SatelliteContainer operator() (size_t i);
  SatelliteContainer& operator<< (Satellite* sat);
  SatelliteContainer& operator<< (Satellite& sta);
  static SatelliteContainer DistributeEvenly (adi::Satellite::Ele ele, size_t N);

  void Install (uint32_t faces, const PointingBound& bound);
  void Install (size_t i, const Face& face, const PointingBound& bound);
private:
  std::vector<Satellite*> m_container;
};

}

#endif /* ADI_SATELLITE_CONTAINER_H */