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

#ifndef ADI_LINK_HELPER_H
#define ADI_LINK_HELPER_H

#include "adi-link.h"
#include "adi-interval.h"

namespace adi {

class LinkHelper
{
public:

  /**
   * \brief Default constructor
   */
  LinkHelper (){}

  /**
   * \brief Deconstructor
   */
  ~LinkHelper (){}

  /**
   * \brief Set the interval to be calculated
   * \param[in] interval the interval
   */
  void SetInterval (const Interval& interval);

  /**
   * \brief Set the interval list to be calculated
   * \param[in] intervals the interval list
   */
  void SetIntervalList (const Intervals& intervals);

  /**
   * \brief Get the interval list
   * \return the interval list
   */
  const Intervals& GetIntervalList () const;

  /**
   * \brief Add the link to be calculated
   * \param[in] link the link
   */
  void AddLink (const Link& link);

  /**
   * \brief Calculate the link data
   * \return the link data
   */
  LinkInfoList& CalcLink ();
private:

  /**
   * \brief Check whether the link exists in the link list
   * \param[in] link the link
   * \return true if the link has existed, false otherwise
   */
  bool IsExist (Link link);

  Intervals m_intervals;  //!< the interval list to be calculated
  Links     m_links;      //!< the link list
  LinkInfoList m_linkDatas;
};
}

#endif /* ADI_LINK_HELPER_H */