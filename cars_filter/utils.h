// Copyright (C) 2024 Centre National d'Etudes Spatiales (CNES).
//
// This file is part of cars-filter
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>
#include <tuple>

namespace cars_filter
{
using PointType = std::array<double, 3>;

/*
* Check if a point is NaN
* 
* \param point Input point
* \return true if any coordinate of the point is nan
*
*/
inline bool isnan(const PointType& point)
{
  return std::isnan(point[0]) || std::isnan(point[1]) || std::isnan(point[2]);
}

// Squared euclidian distance (single coordinates version)
inline double squared_euclidian_distance(const double x1,
                                         const double y1,
                                         const double z1,
                                         const double x2,
                                         const double y2,
                                         const double z2)
{
  const double dx = x1-x2;
  const double dy = y1-y2;
  const double dz = z1-z2;
  return dx * dx + dy * dy + dz * dz;
}


// Squared euclidian distance (Point version)
inline double squared_euclidian_distance(const PointType& point1,
                                         const PointType& point2)
{
  return squared_euclidian_distance(point1[0],
                            point1[1],
                            point1[2],
                            point2[0],
                            point2[1],
                            point2[2]);
}

/*
*
* \brief Compute first quartile, median and third quartile of input vector
* 
* This function returns the closest element (rounded down) to the quantile.
* For example the median of a vector of 30 elements is the 15th element,
* not the mean between the 15th and 16th.
* 
* Note: the vector is copied because nth_element modifies the input container
*/
template<typename T>
std::tuple<T, T, T> compute_approximate_quantiles(std::vector<T> input_data)
{

  // clean list
  std::vector<T> cleaned_data;
  for (const T& value : input_data) {
      if (!std::isnan(value)) {
          cleaned_data.push_back(value);
      }
  }
  unsigned int size = cleaned_data.size();
  unsigned int first_quarter_pos = static_cast<unsigned int>(size/4);
  unsigned int half_pos = static_cast<unsigned int>(size/2);
  unsigned int third_quarter_pos = first_quarter_pos + half_pos;

  if (size == 0) {
        return {std::nan(""), std::nan(""), std::nan("")};
  }

  // Split the container at the median
  std::nth_element(cleaned_data.begin(),
                   cleaned_data.begin() + half_pos,
                   cleaned_data.end());

  // Split the first half of the container at the 0.25 quantile
  std::nth_element(cleaned_data.begin(),
                   cleaned_data.begin() + first_quarter_pos,
                   cleaned_data.begin() + half_pos);

  // Split the other half of the container at the 0.75 quantile
  std::nth_element(cleaned_data.begin() + half_pos+1,
                   cleaned_data.begin() + third_quarter_pos,
                   cleaned_data.end());

  return {cleaned_data[first_quarter_pos],
          cleaned_data[half_pos],
          cleaned_data[third_quarter_pos]};
}


} // namespace cars_filter