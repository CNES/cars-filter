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

#include "small_component_filtering.h"
#include <unordered_set>

namespace cars_filter
{

std::vector<unsigned int> point_cloud_small_component_filtering(
    double* x_coords,
    double* y_coords,
    double* z_coords,
    const unsigned int num_elem,
    const double radius,
    const int min_cluster_size,
    const double clusters_distance_threshold
)
{
  auto input_point_cloud = PointCloud(x_coords, y_coords, z_coords, num_elem);

  auto tree = KDTree(input_point_cloud, 40);

  std::vector<unsigned int> visited_points(num_elem, 0);

  std::vector<unsigned int> result;

  // Iterate on the point cloud nodes, and get the index of a point that have 
  // not been processed yet
  for (const auto& node: tree.getNodes())
  {
    std::vector<unsigned int> indices;
    if (node.m_indices.empty()) 
    {
      indices = {node.m_idx};
    }
    else
    {
      indices = node.m_indices;
    }
    
    for (auto idx: indices)
    {
      if (visited_points[idx])
      {
        continue;
      }

      // Get the neighbors of the point
      auto neighbor_list = tree.neighbors_in_ball(input_point_cloud.m_coords[0][idx], 
                     input_point_cloud.m_coords[1][idx],
                     input_point_cloud.m_coords[2][idx],
                     radius);

      // Add the neighbors to the current cluster
      std::unordered_set<unsigned int> seed;
      for (auto& elem: neighbor_list)
      {
        visited_points[elem] = true;
        seed.insert(elem);
      }

      while (!neighbor_list.empty())
      {
        auto current_idx = neighbor_list.back();
        neighbor_list.pop_back();

        // visited_points[current_idx] = true;
        auto new_neighbors = tree.neighbors_in_ball(input_point_cloud.m_coords[0][current_idx], 
                     input_point_cloud.m_coords[1][current_idx],
                     input_point_cloud.m_coords[2][current_idx],
                     radius);

        for (auto elem: new_neighbors)
        {
          auto [it, is_new_elem] = seed.insert(elem);
          if (is_new_elem)
          {
            visited_points[elem] = true;
            neighbor_list.push_back(elem);
          }
        }

      }

      if (seed.size() < min_cluster_size)
      {
        if (std::isnan(clusters_distance_threshold))
        {
          for (auto elem: seed)
          {
            result.push_back(elem);
          }
        }
        else
        {
          // search if the current cluster has any neighbors in the 
          // clusters_distance_threshold radius
          bool neighbor_found = false;
          for (auto current_idx_it = seed.begin(); current_idx_it != seed.end() && !neighbor_found; current_idx_it++)
          {
            // Optimization note: it would probably be more efficient to mask
            // the pixels in the epipolar neighborhood that also are in the 
            // current cluster, to avoid computing some euclidian distances.
            // This have not been implemented because only a few points should 
            // be processed here, in comparison to the full algorithm and 
            // therefore should not be that imapctful.
            auto new_neighbors = tree.neighbors_in_ball(input_point_cloud.m_coords[0][*current_idx_it],
                                                        input_point_cloud.m_coords[1][*current_idx_it],
                                                        input_point_cloud.m_coords[2][*current_idx_it],
                                                        clusters_distance_threshold);
            
            // Check neighbors of current point and check the set to know if a 
            // new neighbour is found
            for (auto elem: new_neighbors)
            {
              if (!seed.count(elem))
              {
                neighbor_found = true;
                break;
              }
            }
          }
          if (!neighbor_found)
          {
            for (auto elem: seed)
            {
              result.push_back(elem);
            }
          }
        }
      }
    }
  }
  return result;
}

std::vector<std::pair<unsigned int, unsigned int>> epipolar_small_component_filtering(
    Image<double>& x_coords,
    Image<double>& y_coords,
    Image<double>& z_coords,
    Image<double>& outlier_array,
    const unsigned int min_cluster_size,
    const double radius,
    const unsigned int half_window_size,
    const double clusters_distance_threshold)
{
  InMemoryImage<double> visited_pixels(x_coords.number_of_rows(), x_coords.number_of_cols());

  std::vector<std::pair<unsigned int, unsigned int>> clusters;
  std::vector<std::pair<unsigned int, unsigned int>> result;


  // outlier_array initialization
  for (unsigned int row=0; row < x_coords.number_of_rows(); row++)
  {
    for (unsigned int col=0; col < x_coords.number_of_cols(); col++)
    {
      outlier_array.get(row, col) = false;
    }
  }

  for (unsigned int row=0; row < x_coords.number_of_rows(); row++)
  {
    for (unsigned int col=0; col < x_coords.number_of_cols(); col++)
    {
      if (std::isnan(z_coords.get(row, col)))
      {
        continue;
      }

      clusters.push_back({row,col});

      // unordered set might be more efficient, but std pair is not hashable
      // so it would require to work with one index instead of row/col
      std::vector<std::pair<unsigned int, unsigned int>> current_cluster;
      while (!clusters.empty())
      {
        auto [ref_row, ref_col] = clusters.back();
        clusters.pop_back();
        auto& is_ref_visted = visited_pixels.get(ref_row, ref_col);
        if (is_ref_visted)
        {
          continue;
        }
        current_cluster.push_back({ref_row, ref_col});
        is_ref_visted = true;
        auto neighbors = epipolar_neighbors_in_ball(x_coords, y_coords, z_coords, ref_row, ref_col, radius, half_window_size);

        //clusters.reserve(clusters.size() + neighbors.size());
        for (const auto& [elem_row, elem_col] :neighbors)
        {
          if (!visited_pixels.get(elem_row, elem_col))
          {
            clusters.push_back({elem_row, elem_col});
          }
        }
      }

      if (!current_cluster.empty() && current_cluster.size() < min_cluster_size)
      {
        if (std::isnan(clusters_distance_threshold))
        {
          for (const auto& [elem_row, elem_col]:current_cluster)
          {
            result.emplace_back(elem_row, elem_col);
          }
        }
        else
        {
          bool neighbor_found = false;
          for (auto current_idx_it = current_cluster.begin(); current_idx_it != current_cluster.end() && !neighbor_found; current_idx_it++)
          {
            auto& [row, col] = *current_idx_it;
            auto neighbors = epipolar_neighbors_in_ball(x_coords, y_coords, z_coords, row, col, clusters_distance_threshold, half_window_size);
            for (auto& item: neighbors)
            {
              // Check if the new neighbor is already in the current neighborhood
              if (std::find(current_cluster.begin(), current_cluster.end(), item) == current_cluster.end())
              {
                neighbor_found = true;
                break;
              }


            }

          }
          if (!neighbor_found)
          {
            for (const auto& [elem_row, elem_col]:current_cluster)
            {
                result.emplace_back(elem_row, elem_col);
            }
          }
        }
      }
    }
  }
  return result;
}

}