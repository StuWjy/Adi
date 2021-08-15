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

#ifndef ADI_OBJECT_H
#define ADI_OBJECT_H

#include <map>
#include "adi-type-define.h"
#include "adi-interval.h"

namespace adi {

class Object;
class Turntable;
typedef std::vector<Object*> Objects;

class Object {
public:
  enum Type
  {
    STATION,
    SATELLITE
  };

  /**
   * \brief Default constructor
   */
  Object ();

  /**
   * \brief Copy constructor
   */
  Object (const Object& obj);

  /**
   * \brief Deconstructor
   */
  ~Object ();

  /**
   * \brief Set the epoch
   * \param[in] epoch the epoch
   */
  void SetEpoch (const DateTime& epoch);

  /**
   * \return the epoch
   */
  const DateTime& GetEpoch () const;

  /**
   * \brief Set the time step
   * \param[in] step the time step
   */
  void SetStep (const TimeSpan& step);

  /**
   * \return the time step
   */
  const TimeSpan& GetStep () const;

  /**
   * \brief Set the time interval list to be calculated
   * \param[in] intervals the time interval list
   */
  void SetIntervalList (const Intervals& intervals);

  /**
   * \brief Add the time interval to be calculated
   * \param[in] interval the time interval
   * \return the interval list after added
   */
  const Intervals& AddInterval (const Interval& interval);

  /**
   * \return the interval list
   */
  const Intervals& GetIntervalList () const;

  /**
   * \brief Get the string of name
   * \return the name
   */
  virtual std::string GetName () const = 0;

  /**
   * \return Get the object id
   */
  virtual size_t GetId () const = 0;

  /**
   * \return the type of object
   */
  virtual Type GetType () const = 0;

  /**
   * \brief Calculate the trajectory of the object with given interval
   * \param[in] interval the given interval
   */
  virtual Trajectory CalcTrajectory (Matrix** dmat = NULL) = 0;

  /**
   * \brief install the turntable with given face and bound of pointing
   * \param[in] face  the given face to be installed
   * \param[in] bound the given bound of pointing
   * \return the pointer of turntable, NULL -- if has been installed
   */
  virtual Turntable* Install (const Face& face, const PointingBound& bound);
  
  /**
   * \return the pointer of turntable with given face
   */
  Turntable* GetTurntable (const Face& face);

  /**
   * \brief Set the ns3 object
   * \param[in] ns3 the ns3 object
   */
  void SetNs3 (void* ns3);

  /**
   * \return the ns3 object
   */
  void* GetNs3 () const;
protected:
  typedef std::map<size_t, Turntable*> TurntableList;
  virtual void Construct () = 0;
  size_t GetIntervalTicks () const;
  void*           m_ns3;          //!< the object pointer of ns3
  DateTime        m_epoch;        //!< the epoch
  TimeSpan        m_step;         //!< the time step
  Intervals       m_intervals;    //!< the time interval list to be calculated
  size_t          m_uid;          //!< the object id
  TurntableList   m_turntables;   //!< the turntables list
  std::string     m_name;
  // datas in device
  Vector* d_pos;
  Vector* d_vel;
  Matrix* d_mat;
};

}
#endif /* ADI_OBJECT_H */