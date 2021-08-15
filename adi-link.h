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

#ifndef ADI_LINK_H
#define ADI_LINK_H

#include "adi-type-define.h"
#include "adi-interval.h"
#include "adi-constant.h"

namespace adi {

class Turntable;

class Link 
{
private:
  /**
   * \brief Default constructor
   */
  Link ();

public:
  /**
   * \brief Constructor
   * \param[in] src       the pointer of source's turntable
   * \param[in] dst       the pointer of destination's turntable
   * \param[in] included  the state that must be included
   * \param[in] excluded  the state that must be excluded
   */
  Link (
    Turntable* src,
    Turntable* dst,
    uint8_t included = SRC2DST | DST2SRC,
    uint8_t excluded = SRC_DAY | DST_DAY
  );

  /**
   * \brief Deconstructor
   */
  ~Link ();

  /**
   * \brief Set the time interval list to be calculated
   * \param[in] intervals the time interval list
   */
  void SetIntervalList (const Intervals& intervals);

  /**
   * \brief Add the time interval to be calculated
   * \param[in] interval the time interval
   */
  const Intervals& AddInterval (const Interval& interval);

  /**
   * \return the interval list
   */
  const Intervals& GetIntervalList () const;

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
   * \brief Set the maximum distance for calculating link data
   * \param[in] maxDistance the maximum distance
   */
  void SetMaxDistance (const double& maxDistance);

  /**
   * \return the maximum distance
   */
  const double& GetMaxDistance () const;

  /**
   * \brief Calculate the interval list of the link
   * \return the interval list
   */
  Intervals CalcLinkInterval ();
  LinkInfoList CalcLinkData ();
  Satellite* GetSatellite (size_t i);
  Turntable* GetTurntable (size_t i);
  Face GetFace (size_t i) const;
  bool operator== (const Link& link);
  void SetIncludedState (uint8_t included);
  void SetExcludedState (uint8_t excluded);
private:
  void FaceTransform (Vec2Vec& host, Face face);
  bool IsAcceptedState (const uint8_t& state) const;
  DateTime    m_epoch;
  TimeSpan    m_step;
  Intervals   m_intervals;
  Turntable*  m_src;
  Turntable*  m_dst;
  double      m_maxDistance;
  uint8_t     m_includedState;
  uint8_t     m_excludedState;
};

}

#endif /* ADI_LINK_H */