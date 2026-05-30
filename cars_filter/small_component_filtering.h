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


#include "KDTree.h"
#include "image.h"
#include "epipolar_utils.h"

namespace cars_filter
{

/*
*
* \brief Filter a point cloud using the small component method
*
* \param x_coords array containing the x triangulated coordinates
* \param y_coords array containing the y triangulated coordinates
* \param z_coords array containing the altitudes
* \param num_elem number of points in the cloud
* \param radius points are considered connected if their euclidian distance is less than radius
* \param min_cluster_size cluster with less than min_cluster_size will be marked as outliers
* \param clusters_distance_threshold a cluster will not be removed if it has a neighbor within this distance
*
*/
std::vector<unsigned int> point_cloud_small_component_filtering(
    double* x_coords,
    double* y_coords,
    double* z_coords,
    const unsigned int num_elem,
    const double radius = 3,
    const int min_cluster_size = 15,
    const double clusters_distance_threshold = std::numeric_limits<double>::quiet_NaN()
);


/*
*
* \brief Filter an epipolar depth map using the small component method
*
* \param x_coords Image in epipolar geometry containing the x triangulated coordinates
* \param y_coords Image in epipolar geometry containing the y triangulated coordinates
* \param z_coords Image in epipolar geometry containing the z triangulated coordinates
* \param outlier_array output outlier mask (true = outlier)
* \param min_cluster_size cluster with less than min_cluster_size will be marked as outliers
* \param radius points are considered connected if their euclidian distance is less than radius
* \param half window size half size of the epipolar search window (in rows and columns)
* \param clusters_distance_threshold a cluster will not be removed if it has a neighbor within this distance
*
*/
std::vector<std::pair<unsigned int, unsigned int>> epipolar_small_component_filtering(
    Image<double>& x_coords,
    Image<double>& y_coords,
    Image<double>& z_coords,
    Image<double>& outlier_array,
    const unsigned int min_cluster_size = 15,
    const double radius = 10,
    const unsigned int half_window_size = 5,
    const double clusters_distance_threshold = std::numeric_limits<double>::quiet_NaN());

}