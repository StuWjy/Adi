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

#include "adi-sun.h"
#include "adi-sun-kernel.h"
#include "adi-util.h"

namespace adi {

Vector*   Sun::d_pos       = NULL;

Vector*
Sun::CalcPosition (const Intervals& intervals, const DateTime& epoch, const TimeSpan& step)
{
  uint32_t N_TIME = Interval::GetTotalTicks (intervals, step);
  std::vector<int64_t> h_ticks = Interval::CreateTicks (intervals, epoch, step);
  int64_t* d_ticks = NULL;
  cudaMalloc ((void**)&d_ticks, sizeof (int64_t) * N_TIME);
  cudaMemcpy (d_ticks, &h_ticks[0], sizeof (int64_t) * N_TIME, cudaMemcpyHostToDevice);
  cudaFree (d_pos);
  cudaMalloc ((void**)&d_pos, sizeof (Vector) * N_TIME);
  dim3 blockSize = dim3 (std::min (N_TIME, 256u));
  dim3 gridSize = dim3 ((N_TIME + blockSize.x - 1) / blockSize.x);

  cuda::sun::CalcEciPosition <<< gridSize, blockSize >>> (d_pos, d_ticks, N_TIME);
  // {
  //   std::vector<Vector> h_pos = std::vector<Vector> (N_TIME);
  //   cudaMemcpy (&h_pos[0], d_pos, sizeof (Vector) * N_TIME, cudaMemcpyDeviceToHost);
  //   for (const Vector& pos : h_pos)
  //   {
  //     std::cout << pos << std::endl;
  //   }
  // }
  cudaFree (d_ticks);
  cudaCheck;
  return d_pos;
}

}