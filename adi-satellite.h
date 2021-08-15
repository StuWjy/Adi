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

#ifndef ADI_SATELLITE_H
#define ADI_SATELLITE_H

#include <vector>
#include <sstream>
#include <memory>
#include "adi-type-define.h"
#include "adi-object.h"
#include "adi-turntable.h"

namespace adi {

class SatelliteHelper;

class Satellite : public Object
{
public:
  friend class SatelliteHelper;
  struct Element
  {
    double sma; //!< semi-major axis, in km
    double ecc; //!< eccentricity
    double inc; //!< inclination, in radian
    double raan;//!< right ascension of ascending node, in radian
    double aop; //!< argument of perigee, in radian
    double ma;  //!< mean anomaly, in radian
    std::string ToString () const;
    bool operator== (Element& ele) const;
  };

  struct Parameter
  {
    Element elem;
    double ar;    //!< angular rate, in radian per second
    double slr;   //!< semi-latus rectum, in km
    double draan; //!< rate of right ascension of ascending node, in radian per second
    double daop;  //!< rate of argument of perigee, in radian per second
    double dma;   //!< rate of mean anomaly, in radian per second
    std::string ToString () const;
  };

  struct State
  {
    double radius;  //!< current radius, in km
    double raan;    //!< current right ascension of ascending node, in radian
    double aop;     //!< current argument of perigee, in radian
    double ma;      //!< current mean anomaly, in radian
    double ta;      //!< current true anomaly, in radian
    double dta;     //!< current true anomaly rate, in radian per second
    std::string ToString () const;
  };

  typedef State Sta;
  typedef Parameter Par;
  typedef Element Ele;

  /**
   * \brief Default constructor
   */
  Satellite ();

  /**
   * \brief Copy constructor
   */
  Satellite (const Satellite& sat);

  /**
   * \brief Deconstructor
   */
  ~Satellite ();

  /**
   * \brief Set the element of satellite
   * \param[in] element the element of satellite
   */
  void SetElement (Ele ele);

  /**
   * \brief Set the element of satellite
   * \param[in] sma       the semi-major axis, in kilometer
   * \param[in] ecc       the eccentricity, default in degree
   * \param[in] inc       the inclination, default in degree
   * \param[in] raan      the right ascension of ascending node, default in degree
   * \param[in] aop       the argument of perigee, default in degree
   * \param[in] ma        the mean anomaly, default in degree
   * \param[in] isDegree  whether the inc, raan, aop and ma are in degree or radian
   */
  void SetElement (
    double sma,
    double ecc,
    double inc,
    double raan,
    double aop,
    double ma,
    bool   isDegree = true);

  virtual std::string GetName () const;

  /**
   * \return the global id of satellite
   */
  size_t GetId () const;

  /**
   * \return the type of object
   */
  virtual Type GetType () const;

  /**
   * \brief Calculate the trajectory of the satellite with given interval
   * \param[in] interval the given interval
   */
  virtual Trajectory CalcTrajectory (Matrix** dmat = NULL);

  /**
   * \brief Get the pointer of turntable with given face
   * \param[in] face the given face that the turntable is installed on
   * \return the pointer of turntable, NULL if no turntable is installed on given face
   */
  virtual Turntable* GetTurntable (Face face);
private:
  virtual void Construct ();
  void Initialzie ();
  void Release ();
  // datas in host
  Par   m_par;  //!< orbital element of satellite
  //datas in device
  Par*  d_par;
  Sta*  d_sta;
};

std::ostream& operator<< (std::ostream& os, const Satellite::Ele& ele);
std::ostream& operator<< (std::ostream& os, const Satellite::Par& par);
std::ostream& operator<< (std::ostream& os, const Satellite::Sta& sta);

}
#endif /* ADI_SATELLITE_H */