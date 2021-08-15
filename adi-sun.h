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

#ifndef ADI_SUN_H
#define ADI_SUN_H

#include "adi-type-define.h"
#include "adi-interval.h"

namespace adi {

class Sun
{
public:
  friend class Access;
  friend class Link;
  Sun (){};
  ~Sun (){};

  /**
   * \brief Calculate the position of the sun with given interval
   * \param[in] interval the given interval
   * \return the pointer of vector in device
   */
  static Vector* CalcPosition (const Intervals& intervals, const DateTime& epoch, const TimeSpan& step);

private:
  static Vector*  d_pos;
};

}

#endif /* ADI_SUN_H */