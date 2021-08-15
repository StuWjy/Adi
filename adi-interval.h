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

#ifndef ADI_INTERVAL_H
#define ADI_INTERVAL_H

#include <vector>
#include "adi-date-time.h"

namespace adi {

class Interval;
typedef std::vector<Interval> Intervals;

class Interval
{
public:
  /**
   * \brief Default constructor
   */
  Interval ();

  /**
   * \brief Copy constructor
   */
  Interval (const Interval& interval);

  /**
   * \brief Constructor
   * \param[in] _start the start of interval 
   * \param[in] _stop  the stop of interval
   */
  Interval (DateTime _start, DateTime _stop);

  /**
   * \brief Deconstructor
   */
  ~Interval (){}

  /**
   * \brief Set the start of interval
   * \param[in] start the start of interval
   */
  void SetStart (const DateTime& start);

  /**
   * \return the start of interval
   */
  const DateTime& GetStart () const;

  /**
   * \brief Set the stop of interval
   * \param[in] stop the stop of interval
   */
  void SetStop (const DateTime& stop);

  /**
   * \return the stop of interval
   */
  const DateTime& GetStop () const;

  /**
   * \brief Shift the interval by given delay
   * \param[in] delay the given delay
   * \return the shifted interval
   */
  Interval& Shift (const TimeSpan& delay);

  /**
   * \brief Convert to string
   * \return the string of interval
   */
  std::string ToString () const;

  /**
   * \return true if start and stop of intervals are same
   */
  bool operator== (const Interval& interval) const;

  /**
   * \return true if start and stop of intervals are not totally same
   */
  bool operator!= (const Interval& interval) const;

  /**
   * \return true if this start is less than the start of interval
   */
  bool operator<  (const Interval& interval) const;

  /**
   * \brief Whether the interval intersects with the given interval
   * \param[in] interval the given interval
   * \return true if they intersect with each other
   */
  bool Intersect  (const Interval& interval) const;

  /**
   * \brief Whether the given date time is between the interval
   * \param[in] dateTime the given date time
   * \return true if the interval intersects with the given date time
   */
  bool Intersect  (const DateTime& dateTime) const;

  /**
   * \brief Merge the interval with given interval
   * \param[in] interval the given interval
   */
  void Merge (const Interval& interval);

  /**
   * \brief Insert the interval into the interval list
   * \param[inout] intervals the interval list
   */
  void InsertToList (Intervals& intervals) const;

  /**
   * \brief Get the ticks of the interval
   * \param[in] step the time step
   * \return the ticks of the interval
   */
  size_t GetTicks (const TimeSpan& step) const;

  /**
   * \brief Get total ticks of the interval list with given time step
   * \param[in] intervals the interval list
   * \param[in] step      the time step
   * \return the total ticks
   */
  static size_t GetTotalTicks (const Intervals& intervals, const TimeSpan& step);

  /**
   * \brief Create the elapsed seconds according the given interval list and time step, 0 is the epoch
   * \param[in] intervals the interval list
   * \param[in] epoch     the epoch
   * \param[in] step      the time step
   * \return    the vector of elapsed seconds
   */
  static std::vector<double> CreateSeconds (const Intervals& intervals, const DateTime& epoch, const TimeSpan& step);

  /**
   * \brief Create the ticks according the given interval list and time step
   * \param[in] intervals the interval list
   * \param[in] epoch     the epoch
   * \param[in] step      the time step
   * \return    the vector of ticks
   */
  static std::vector<int64_t> CreateTicks (const Intervals& intervals, const DateTime& epoch, const TimeSpan& step);
private:
  DateTime m_start; //!< the start of interval
  DateTime m_stop;  //!< the stop of interval
};

std::ostream& operator<< (std::ostream& os, const Interval& i);

} //namespace adi

#endif  /* ADI_INTERVAL_H */