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

#ifndef ADI_CONSTANT_H
#define ADI_CONSTANT_H

#include "adi-type-define.h"
#include "adi-time-span.h"

namespace adi {

const double K_MU = 398600.5;                   //!< Gravitational constant
const double K_RE = 6378.137;                   //!< Earth radius
const double K_J2 = 1082.63e-6;                 //!< J2 constant
const double K_EAV = 1.00273779174;             //!< Earth rotation angle per sidereal day
const double K_AU = 1.49597870691e8;            //!< Astronomical Unit
const double K_F = 1.0 / 298.257223563;         //!< Earth flattening
const double K_EAR = M_PI * (K_EAV / 86400.0);  //!< rate of Earth spin
const double K_MAX_DISTANCE = 1500.0;


const DateTime J2000  = DateTime (2000, 1, 1, 12, 0, 0);
const TimeSpan Day    = TimeSpan (1, 0, 0, 0);
const TimeSpan Hour   = TimeSpan (1, 0, 0);
const TimeSpan Minute = TimeSpan (0, 1, 0);
const TimeSpan Second = TimeSpan (0, 0, 1);

#define D2R(x) (x * M_PI / 180.0)
#define R2D(x) (x * 180.0 / M_PI)

const PointingBound K_FOV = PointingBound {
  Bound {D2R (-85.0), D2R (85.0)},
  Bound {D2R (-85.0), D2R (85.0)}
};

#define BEYOND_DISTANCE 0x01
#define SRC_DAY         0x80
#define SRC2DST         0x40
#define DST_DAY         0x08
#define DST2SRC         0x04
#define IN_FOV          0x10
}
#endif /* ADI_CONSTANT_H */