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

#ifndef ADI_STATION_H
#define ADI_STATION_H

#include <vector>
#include <sstream>
#include <memory>
#include "adi-type-define.h"
#include "adi-object.h"
#include "adi-turntable.h"

namespace adi {
class Station : public Object
{
public:
  friend class StationHelper;
  struct Element
  {
    double alt; //!< altitude, in km
    double lat; //!< latitude, in radian
    double lon; //!< longitude, in radian
    std::string ToString () const;
  };

  typedef Element Ele;

  /**
   * \brief Default constructor
   */
  Station ();
  
  /**
   * \brief Copy constructor
   */
  Station (const Station& sta);

  /**
   * \brief Deconstructor
   */
  ~Station ();

  /**
   * \brief Set the element of station
   * \param[in] element the element of station
   */
  void SetElement (Ele ele);

  /**
   * \brief Set the element of station
   * \param[in] lat the latitude, default in degree
   * \param[in] lon the longitude, default in degree
   * \param[in] alt the altitude, in kilometer
   * \param[in] name the name of station
   * \param[in] isDegree whether the latitude or longitude is in degree or radian
   */
  void SetElement (
    double lat,
    double lon,
    double alt,
    std::string name = "NoName",
    bool isDegree = true);

  // void SetElement (
  //   double lat,
  //   double lon,
  //   double alt,
  //   bool isDegree = true);

  virtual std::string GetName () const;

  /**
   * \return the global id of station
   */
  size_t GetId () const;

  Type GetType () const;
  
  virtual Trajectory CalcTrajectory (Matrix** dmat = NULL);
  virtual Turntable* GetTurntable (Face face = Top);
private:
  virtual void Construct ();
  void Initialize ();
  void Release ();
  // datas in host
  Ele   m_ele;   //!< orbital element of station
  // datas in device
  Ele*  d_ele;
  Gst*  d_gst;
};

std::ostream& operator<< (std::ostream& os, Station::Ele& elem);

}
#endif /* ADI_STATION_H */