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

#ifndef ADI_TYPE_DEFINE_H
#define ADI_TYPE_DEFINE_H

#include <vector>
#include <sstream>
#include <iomanip>
#include "adi-date-time.h"

namespace adi {

struct Vector
{
  double x;
  double y;
  double z;
  std::string ToString () const;
};

struct Matrix
{
  Vector r1;
  Vector r2;
  Vector r3;
  std::string ToString () const;
};

template<class T>
struct Array
{
  T*        ptr;
  uint32_t  num;
};

typedef Array<Vector> Vectors;
typedef Array<Matrix> Matrices;
typedef Vectors Position;
typedef Vectors Velocity;
typedef Array<bool> Shadows;

struct Trajectory
{
  Vector* pos;
  Vector* vel;
  Matrix* mat;
  size_t  num;
  std::string ToString () const;
};

typedef std::vector<Trajectory> Trajectories;

struct Direction
{
  double azimuth;
  double pitch;
  // std::string ToString () const;
};

struct Pointing
{
  Direction angle;
  Direction rate;
  // std::string ToString () const;
};

typedef std::vector<Pointing> Pointings;

struct Bound
{
  double min;
  double max;
};

struct PointingBound
{
  Bound azimuth;
  Bound pitch;
  std::string ToString () const;
};

typedef double Gst;

enum Face
{
  ErrFace = 0,
  Left    = 1 << 0,
  Right   = 1 << 1,
  Front   = 1 << 2,
  Back    = 1 << 3,
  Top     = 1 << 4,
  Bottom  = 1 << 5
};

std::string FaceToString (Face face);

class Turntable;

typedef std::vector<DateTime> DateTimes;

typedef std::vector<uint8_t> States;

struct LinkData
{
  uint8_t  state;
  DateTime time;
  Pointing fromSrc;
  Pointing fromDst;
  double   distance;
};

typedef std::vector<LinkData> LinkDatas;

struct LinkInfo
{
  Turntable*  src;
  Turntable*  dst;
  LinkDatas   linkDatas;
  std::string ToString () const;
};

typedef std::vector<LinkInfo> LinkInfoList;

class Satellite;
typedef std::vector<Satellite> Satellites;
class Station;
typedef std::vector<Station> Stations;
class Link;
typedef std::vector<Link> Links;

typedef Vector (*Vec2Vec)(const Vector&);

// #define D2R DegreesToRadians
// #define R2D RadiansToDegrees

std::ostream& operator<< (std::ostream& os, const Vector& v);

std::ostream& operator<< (std::ostream& os, const Matrix& v);

std::ostream& operator<< (std::ostream& os, const Trajectory& t);

// std::ostream& operator<< (std::ostream& os, const Pointing& p);

std::ostream& operator<< (std::ostream& os, const LinkInfo& l);

}

struct dim3;
std::ostream& operator<< (std::ostream& os, const dim3& d);

#endif /* ADI_TYPE_DEFINE_H */