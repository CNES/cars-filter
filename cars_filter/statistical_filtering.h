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

#include <vector>
#include "KDTree.h"
#include "epipolar_utils.h"

namespace cars_filter
{

/*
*
* \brief Filter a point cloud using the statistical method
*
* \param x_coords array containing the x triangulated coordinates
* \param y_coords array containing the y triangulated coordinates
* \param z_coords array containing the altitudes
* \param num_elem number of points in the cloud
* \param filtering_constant constant added to the distance threshold distance threshold
* \param mean_factor ratio applied to mean or median distance in distance threshold
* \param dev_factor ratio applied to stddev or interquartile distance in distance threshold
* \param k number of neighbors in the KNN algorithm
* \param use_median Use median+interquartile distance (true) or mean+stddev (false) to compute distance threshold
*
*/
std::vector<unsigned int> statistical_filtering(double* x_coords,
                                                double* y_coords,
                                                double* z_coords,
                                                const unsigned int num_elem,
                                                const double filtering_constant = 0,
                                                const double mean_factor = 1,
                                                const double dev_factor = 1.,
                                                const unsigned int k = 50,
                                                const bool use_median = false);

/*
*
* \brief Filter an epipolar depth map using the statistical method
*
* \param x_coords Image in epipolar geometry containing the x triangulated coordinates
* \param y_coords Image in epipolar geometry containing the y triangulated coordinates
* \param z_coords Image in epipolar geometry containing the z triangulated coordinates
* \param outlier_array output outlier mask (true = outlier)
* \param k number of neighbors in the KNN algorithm
* \param half window size half size of the epipolar search window (in rows and columns)
* \param filtering_constant constant added to the distance threshold distance threshold
* \param mean_factor ratio applied to mean or median distance in distance threshold
* \param dev_factor ratio applied to stddev or interquartile distance in distance threshold
* \param use_median Use median+interquartile distance (true) or mean+stddev (false) to compute distance threshold
*
*/
void epipolar_statistical_filtering(Image<double>& x_coords,
                                    Image<double>& y_coords,
                                    Image<double>& z_coords,
                                    Image<double>& outlier_array,
                                    const unsigned int k = 50,
                                    const unsigned int half_window_size = 15,
                                    const double filtering_constant = 0,
                                    const double mean_factor = 1,
                                    const double dev_factor = 1,
                                    const double use_median = false);


} // namespace cars_filter