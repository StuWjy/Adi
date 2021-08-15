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

#include "adi-interval.h"

namespace adi {

Interval::Interval ()
: m_start (0)
, m_stop  (0)
{
}

Interval::Interval (const Interval& interval)
: m_start (interval.m_start)
, m_stop  (interval.m_stop)
{
}

Interval::Interval (DateTime _start, DateTime _stop)
: m_start (_start)
, m_stop  (_stop)
{
  try
  {
    if (_start > _stop)
    {
      throw std::exception ();
    }
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    m_start = _stop;
    m_stop  = _start;
  }
}

void
Interval::SetStart (const DateTime& start)
{
  m_start = start;
}

const DateTime&
Interval::GetStart () const
{
  return m_start;
}

void
Interval::SetStop (const DateTime& stop)
{
  m_stop = stop;
}

const DateTime&
Interval::GetStop () const
{
  return m_stop;
}

Interval&
Interval::Shift (const TimeSpan& delay)
{
  m_start = m_start + delay;
  m_stop  = m_stop + delay;
  return *this;
}

std::string
Interval::ToString () const
{
  std::stringstream ss;
  ss << "Interval: "            << std::endl
     << "    from: " << m_start << std::endl
     << "      to: " << m_stop  << std::endl;
  return ss.str ();
}

bool
Interval::operator== (const Interval& interval) const
{
  return m_start == interval.m_start && m_stop == interval.m_stop;
}

bool
Interval::operator!= (const Interval& interval) const
{
  return !(*this == interval);
}


bool
Interval::operator<  (const Interval& interval) const
{
  return m_start < interval.m_start;
}

bool
Interval::Intersect  (const Interval& interval) const
{
  return Intersect (interval.m_start) || Intersect (interval.m_stop);
}

bool
Interval::Intersect  (const DateTime& dateTime) const
{
  return dateTime <= m_stop && dateTime >= m_start;
}

void
Interval::Merge (const Interval& interval)
{
  m_start = std::min (m_start, interval.m_start);
  m_stop = std::max (m_stop, interval.m_stop);
}

void
Interval::InsertToList (Intervals& intervals) const
{
  Interval interval = *this;
  // the vector of intervals' indices that intersect with this interval
  std::vector<uint32_t> intersectedIndex;
  // find all the satisfied indices
  for (uint32_t i = 0;i < intervals.size ();++i)
  {
    if (interval.Intersect (intervals[i]))
    {
      intersectedIndex.push_back (i);
    }
  }
  /**
   * if there are no intersected indices
   * insert this interval and sort the vector
   */
  if(intersectedIndex.empty ())
  {
    intervals.push_back (interval);
    std::sort (intervals.begin (), intervals.end ());
    return;
  }
  else
  {
    // merge all the intersected intervals with this interval
    for (uint32_t i = 0;i < intersectedIndex.size ();++i)
    {
      interval.Merge (intervals[intersectedIndex[i]]);
    }
    // erase all the intersected intervals
    for (uint32_t i = 0;i < intersectedIndex.size ();++i)
    {
      intervals.erase (intervals.begin () + intersectedIndex[i]);
      // for the remain indices, decrease the indices by 1
      for (uint32_t j = i + 1;j < intersectedIndex.size ();++j)
      {
        intersectedIndex[j]--;
      }
    }
    // insert this interval and sort the vector
    intervals.push_back (interval);
    std::sort (intervals.begin (), intervals.end ());
    return;
  }
}

size_t
Interval::GetTicks (const TimeSpan& step) const
{
  return (m_stop - m_start).TotalSeconds () / step.TotalSeconds () + 1;
}

size_t
Interval::GetTotalTicks (const Intervals& intervals, const TimeSpan& step)
{
  size_t ticks = 0;
  for (const Interval& interval : intervals)
  {
    ticks += interval.GetTicks (step);
  }
  return ticks;
}

std::vector<double>
Interval::CreateSeconds (const Intervals& intervals, const DateTime& epoch, const TimeSpan& step)
{
  std::vector<double> seconds;
  uint32_t k = 0;
  for (uint32_t i = 0;i < intervals.size ();++i)
  {
    for (uint32_t j = 0;j < intervals[i].GetTicks (step);++j)
    {
      seconds.push_back ((intervals[i].m_start.AddTicks (step.Ticks () * j) - epoch).TotalSeconds ());
    }
  }
  return seconds;
}

std::vector<int64_t>
Interval::CreateTicks (const Intervals& intervals, const DateTime& epoch, const TimeSpan& step)
{
  std::vector<int64_t> ticks;
  uint32_t k = 0;
  for (uint32_t i = 0;i < intervals.size ();++i)
  {
    for (uint32_t j = 0;j < intervals[i].GetTicks (step);++j)
    {
      ticks.push_back ((intervals[i].m_start.AddTicks (step.Ticks () * j)).Ticks ());
    }
  }
  return ticks;
}

std::ostream& operator<< (std::ostream& os, const Interval& i)
{
  os << i.ToString ();
  return os;
}

} //namespace adi