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

#include "statistical_filtering.h"

namespace cars_filter
{

std::vector<unsigned int> statistical_filtering(double* x_coords,
                                                double* y_coords,
                                                double* z_coords,
                                                const unsigned int num_elem,
                                                const double filtering_constant,
                                                const double mean_factor,
                                                const double dev_factor,
                                                const unsigned int k,
                                                const bool use_median)
{
  auto input_point_cloud = PointCloud(x_coords, y_coords, z_coords, num_elem);

  auto tree = KDTree(input_point_cloud, 40);

  std::vector<double> mean_distances(num_elem, 0);
  for (auto& node: tree.getNodes())
  {
    // Leaf Node case
    if (node.m_indices.empty())
    {
      auto k_neighbors = tree.findKNNIterative(input_point_cloud.m_coords[0][node.m_idx],
                                               input_point_cloud.m_coords[1][node.m_idx],
                                               input_point_cloud.m_coords[2][node.m_idx],
                                               k+1,
                                               &node);
      double acc =0;
      for (const auto& elem : k_neighbors)
      {
        acc += std::sqrt(elem.distance);
      }
      mean_distances[node.m_idx] = acc/k;
    }
    else
    {
      for (auto idx: node.m_indices)
      {
        //const PointType point = {input_point_cloud.m_x[idx], input_point_cloud.m_y[idx], input_point_cloud.m_z[idx]};
        auto k_neighbors = tree.findKNNIterative(input_point_cloud.m_coords[0][idx],
                                                 input_point_cloud.m_coords[1][idx],
                                                 input_point_cloud.m_coords[2][idx],
                                                 k+1,
                                                 &node);
        double acc =0;
        for (const auto& elem : k_neighbors)
        {
          acc += std::sqrt(elem.distance);
        }
        mean_distances[idx] = acc/k;
      }
    }
  }

  double dist_thresh;

  if (use_median)
  {
    auto [percentile_25, median, percentile_75] = compute_approximate_quantiles(mean_distances);
    double interquartile_distance = percentile_75 - percentile_25;
    dist_thresh = filtering_constant + mean_factor * median + dev_factor * interquartile_distance;
  }
  else
  {
    // mean and variance mode
    const double mean = std::accumulate(mean_distances.begin(), mean_distances.end(), 0.0) / num_elem;

    // variance (mean)
    auto variance_lambda= [mean](double acc, double input) {
        return acc + (input - mean)*(input - mean);
    };

    const double var = std::accumulate(mean_distances.begin(), mean_distances.end(), 0.0, variance_lambda)/num_elem;
    const double stddev = std::sqrt(var);
    dist_thresh = filtering_constant + mean_factor * mean + dev_factor * stddev;
  }

  std::vector<unsigned int> result;

  // list version
  for (unsigned int i = 0; i<num_elem; i++)
  {
    if (mean_distances[i] > dist_thresh)
    {
      result.push_back(i);
    }
  }

  return result;
}


void epipolar_statistical_filtering(Image<double>& x_coords,
                                    Image<double>& y_coords,
                                    Image<double>& z_coords,
                                    Image<double>& outlier_array,
                                    const unsigned int k,
                                    const unsigned int half_window_size,
                                    const double filtering_constant,
                                    const double mean_factor,
                                    const double dev_factor,
                                    const double use_median)
{
  InMemoryImage<double> mean_distance_image(x_coords.number_of_rows(), x_coords.number_of_cols());

  for (unsigned int row=0; row < x_coords.number_of_rows(); row++)
  {
    for (unsigned int col=0; col < x_coords.number_of_cols(); col++)
    {
      auto mean_dist = epipolar_knn(x_coords, y_coords, z_coords, row, col, k+1, half_window_size);
      mean_distance_image.get(row, col) = mean_dist;
    }
  }

  double distance_threshold;

  if (use_median)
  {
    std::vector<double> mean_distances_no_nan;
    std::copy_if(mean_distance_image.getData(), 
                mean_distance_image.getData() + mean_distance_image.size(),
                std::back_inserter(mean_distances_no_nan),
                [](double elem){return !std::isnan(elem);});
    auto [percentile_25, median, percentile_75] = compute_approximate_quantiles(mean_distances_no_nan);
    double interquartile_distance = percentile_75 - percentile_25;
;
    distance_threshold = filtering_constant + mean_factor * median + dev_factor * interquartile_distance;
  }
  else
  {
    // mean with number of valid pixel counting
    unsigned int num_valid = 0;
    auto mean_lambda = [&num_valid](double acc, double input)
    {
      if (std::isnan(input))
      {
        return acc;
      }
      else
      {
        num_valid++;
        return acc + input;
      }
    };

    auto mean = std::accumulate(mean_distance_image.getData(), mean_distance_image.getData()+mean_distance_image.size(), 0., mean_lambda);
    mean/=num_valid;

    // variance (ignoring nan)
    auto variance_lambda = [mean](double acc, double input)
    {
      if (std::isnan(input))
      {
        return acc;
      }
      else
      {
        return acc + (input - mean)*(input - mean);
      }
    };

    const double var = std::accumulate(mean_distance_image.getData(),
                                       mean_distance_image.getData() + mean_distance_image.size(),
                                       0.,
                                       variance_lambda) / num_valid;
    const double stddev = std::sqrt(var);

    distance_threshold = filtering_constant + mean_factor * mean + dev_factor * stddev;
  }

  for (unsigned int row=0; row < x_coords.number_of_rows(); row++)
  {
    for (unsigned int col=0; col < x_coords.number_of_cols(); col++)
    {
      if (mean_distance_image.get(row, col) > distance_threshold)
      {
        outlier_array.get(row, col) = true;
        x_coords.get(row, col) = std::numeric_limits<double>::quiet_NaN();
        y_coords.get(row, col) = std::numeric_limits<double>::quiet_NaN();
        z_coords.get(row, col) = std::numeric_limits<double>::quiet_NaN();
      }
      else
      {
        outlier_array.get(row, col) = false;
      }
    }
  }
}

} // namespace cars_filter
