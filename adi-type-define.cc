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

#include <cuda_runtime.h>
#include "adi-type-define.h"
#include "adi-turntable.h"
#include "adi-constant.h"

namespace adi {

extern __host__ __device__ double RadiansToDegrees (double x);

std::string
Vector::ToString () const
{
  std::stringstream ss;
  ss.setf (std::ios::fixed);
  ss.precision(4);
  ss << "("
      << std::setw (10) << x << ", "
      << std::setw (10) << y << ", "
      << std::setw (10) << z
      << ")";
  return ss.str ();
}


std::string
Matrix::ToString () const
{
  std::stringstream ss;
  ss.setf (std::ios::fixed);
  ss.precision(4);
  ss << " --------------------------------" << std::endl;
  ss << "|" << std::setw (10) << r1.x << " " << std::setw (10) << r1.y << " " << std::setw (10) << r1.z << "|" << std::endl;
  ss << "|" << std::setw (10) << r2.x << " " << std::setw (10) << r2.y << " " << std::setw (10) << r2.z << "|" << std::endl;
  ss << "|" << std::setw (10) << r3.x << " " << std::setw (10) << r3.y << " " << std::setw (10) << r3.z << "|" << std::endl;
  ss << " --------------------------------";
  return ss.str ();
}

std::string
Trajectory::ToString () const
{
  std::stringstream ss;
  for (size_t i = 0;i < num;i++)
  {
    ss << std::setw (8) << i << " | Position: " << pos[i].ToString () << " Velocity: " << vel[i].ToString () << std::endl;
  }
  return ss.str ();
}

// std::string
// Direction::ToString () const
// {
//   std::stringstream ss;
//   ss.precision (4);
//   ss.setf(std::ios::fixed);
//   ss << "Azimuth: " << std::setw (10) << RadiansToDegrees (azimuth) << std::endl
//      << "  Pitch: " << std::setw (10) << RadiansToDegrees (pitch)   << std::endl;
//   return ss.str ();
// }

std::string
PointingBound::ToString () const
{
  std::stringstream ss;
  ss.precision (4);
  ss.setf (std::ios::fixed);
  ss << "Azimuth min: " << azimuth.min * 180.0 / M_PI  << " max: " << azimuth.max * 180.0 / M_PI << std::endl
     << "  Pitch min: " << pitch.min * 180.0 / M_PI   << " max: " << pitch.max * 180.0 / M_PI   << std::endl;
    return ss.str ();
}

std::string FaceToString (Face face)
{
  switch (face)
  {
    case Left:    return "Left";
    case Right:   return "Right";
    case Front:   return "Front";
    case Back:    return "Back";
    case Top:     return "Top";
    case Bottom:  return "Bottom";
    default:      return "Error face";
  }
}

std::string TypeToString (Object::Type type)
{
  switch (type)
  {
    case Object::STATION:   return "Station";
    case Object::SATELLITE: return "Satellite";
    default:                return "Error type";
  }
}

std::string StateToString (uint8_t state)
{
  std::string s;
  if (state & IN_FOV)
  {
    s = s + "IN_FOV|";
  }
  if (state & BEYOND_DISTANCE)
  {
    s = s + "BEYOND|";
  }
  if (state & SRC2DST)
  {
    s = s + "SRC2DST|";
  }
  if (state & DST2SRC)
  {
    s = s + "DST2SRC|";
  }
  if (state & SRC_DAY)
  {
    s = s + "SRC_DAY|";
  }
  if (state & DST_DAY)
  {
    s = s + "DST_DAY|";
  }
  return s;
}

std::string
LinkInfo::ToString () const
{
  if (linkDatas.size () ==0)
  {
    return "";
  }
  std::stringstream ss;
  ss.precision (4);
  ss.setf(std::ios::fixed);
  ss << std::setw (30) << " "
     << "|"
     << std::setw (35) << "Start"
     << "|"
     << std::setw (35) << "Stop"
     << "|"
     << std::endl
     << std::setw (30) << " "
     << "|"
     << std::setw (35) << linkDatas.front ().time
     << "|"
     << std::setw (35) << linkDatas.back ().time
     << "|"
     << std::endl
     << std::setw (30) << linkDatas.size ()
     << "|"
     << std::setw (35) << "Source:"
     << "|"
     << std::setw (35) << "Destination:"
     << "|"
     << std::endl
     << std::setw (30) << "Type:"
     << "|"
     << std::setw (35) << TypeToString (src->GetObject ()->GetType ())
     << "|"
     << std::setw (35) << TypeToString (dst->GetObject ()->GetType ())
     << "|"
     << std::endl
     << std::setw (30) << "Id:"
     << "|"
     << std::setw (35) << src->GetObject ()->GetId ()
     << "|"
     << std::setw (35) << dst->GetObject ()->GetId ()
     << "|"
     << std::endl
     << std::setw (30) << "Direction: "
     << "|"
     << std::setw (35) << FaceToString (src->GetFace ())
     << "|"
     << std::setw (35) << FaceToString (dst->GetFace ())
     << "|"
     << std::endl
     << "Time                          | Azimuth|    rate|   Pitch|    rate| Azimuth|    rate|   Pitch|    rate|  distance|state:"
     << std::endl;
  for (size_t i = 0;i < linkDatas.size ();++i)
  {
    ss << linkDatas[i].time << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromSrc.angle.azimuth)  << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromSrc.rate.azimuth)   << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromSrc.angle.pitch)    << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromSrc.rate.pitch)     << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromDst.angle.azimuth)  << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromDst.rate.azimuth)   << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromDst.angle.pitch)    << "|"
       << std::setw (8) << RadiansToDegrees (linkDatas[i].fromDst.rate.pitch)     << "|"
       << std::setw (10) << linkDatas[i].distance << "|"
       << StateToString (linkDatas[i].state)
       << std::endl;
  }
  return ss.str ();
}


std::ostream& operator<< (std::ostream& os, const Vector& v)
{
  os << v.ToString ();
  return os;
}

std::ostream& operator<< (std::ostream& os, const Matrix& m)
{
  os << m.ToString ();
  return os;
}

std::ostream& operator<< (std::ostream& os, const Trajectory& t)
{
  os << t.ToString ();
  return os;
}

// std::ostream& operator<< (std::ostream& os, const Pointing& p)
// {
//   os << p.ToString ();
//   return os;
// }

std::ostream& operator<< (std::ostream& os, const LinkInfo& l)
{
  os << l.ToString ();
  return os;
}

}

std::ostream& operator<< (std::ostream& os, const dim3& d)
{
  os << "(" << d.x << ", " << d.y << ", " << d.z << ")";
  return os;
}
