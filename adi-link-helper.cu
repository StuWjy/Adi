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

#include "adi-link-helper.h"
#include "adi-constant.h"
#include "adi-turntable.h"

namespace adi {

void
LinkHelper::SetInterval (const Interval& interval)
{
  m_intervals.clear ();
  m_intervals.push_back (interval);
}

void
LinkHelper::SetIntervalList (const Intervals& intervals)
{
  m_intervals = intervals;
}

const Intervals&
LinkHelper::GetIntervalList () const
{
  return m_intervals;
}

void
LinkHelper::AddLink (const Link& link)
{
  if (!IsExist (link))
  {
    m_links.push_back (link);
  }
}

LinkInfoList&
LinkHelper::CalcLink ()
{
  m_linkDatas.clear ();
  for (Link& link : m_links)
  {
    link.SetStep (Minute);
    link.SetIntervalList (m_intervals);
    Intervals intervals = link.CalcLinkInterval ();
    if (intervals.empty ())
    {
      continue;
    }
    for (Interval& interval : intervals)
    {
      interval.SetStart (std::max (m_intervals.front ().GetStart (), interval.GetStart () - Minute));
      interval.SetStop (std::min (m_intervals.back ().GetStop (), interval.GetStop () + Minute));
    }
    link.SetStep (Second);
    link.SetIntervalList (intervals);
    LinkInfoList linkData = link.CalcLinkData ();
    m_linkDatas.insert (m_linkDatas.end (), linkData.begin (), linkData.end ());
  }
  return m_linkDatas;
}

bool
LinkHelper::IsExist (Link link)
{
  for (Link& _link : m_links)
  {
    if (_link == link)
    {
      return true;
    }
  }
  return false;
}

}