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
#include "statistical_filtering.h"
#include "image.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "KDTree.h"

namespace py = pybind11;

template<typename T>
cars_filter::Image<T> pyarray_to_image(py::array_t<T> input_array)
{
  py::buffer_info info = input_array.request();
  return cars_filter::Image<T>(static_cast<T*>(info.ptr), info.shape[0], info.shape[1]);
}


py::array_t<double, py::array::c_style> pyEpipolarStatisticalOutlierFiltering(
        py::array_t<double, py::array::c_style>& x_values,
        py::array_t<double, py::array::c_style>& y_values,
        py::array_t<double, py::array::c_style>& z_values,
        const unsigned int k,
        const unsigned int half_window_size,
        const double filtering_constant,
        const double mean_factor,
        const double dev_factor,
        const double use_median
)
{
  // Build image wrappers around numpy arrays
  auto x_image = pyarray_to_image<double>(x_values);
  auto y_image = pyarray_to_image<double>(y_values);
  auto z_image = pyarray_to_image<double>(z_values);

  // // Create a Python object
  py::array_t<double, py::array::c_style> outlier_array({x_image.number_of_rows(),  x_image.number_of_cols()},
                                                        {x_image.number_of_cols() * sizeof(double),  sizeof(double)});
  auto outlier_image = pyarray_to_image<double>(outlier_array);

  epipolar_statistical_filtering(x_image,
                                 y_image,
                                 z_image,
                                 outlier_image,
                                 k,
                                 half_window_size,
                                 filtering_constant,
                                 mean_factor,
                                 dev_factor,
                                 use_median);

  return outlier_array;
}



py::array_t<double, py::array::c_style> pyEpipolarSmallComponentOutlierFiltering(
        py::array_t<double, py::array::c_style>& x_values,
        py::array_t<double, py::array::c_style>& y_values,
        py::array_t<double, py::array::c_style>& z_values,
        const unsigned int min_cluster_size = 15,
        const double radius = 10,
        const unsigned int half_window_size = 5,
        const double clusters_distance_threshold = 4
)
{
  // Build image wrappers around numpy arrays
  auto x_image = pyarray_to_image<double>(x_values);
  auto y_image = pyarray_to_image<double>(y_values);
  auto z_image = pyarray_to_image<double>(z_values);

  // Create a Python object
  py::array_t<double, py::array::c_style> outlier_array({x_image.number_of_rows(),  x_image.number_of_cols()},
                                                        {x_image.number_of_cols() * sizeof(double),  sizeof(double)});
  auto outlier_image = pyarray_to_image<double>(outlier_array);

  epipolar_small_component_filtering(x_image, y_image, z_image, outlier_image, min_cluster_size, radius, half_window_size, clusters_distance_threshold);

  return outlier_array;
}




py::list pyPointCloudStatisticalOutlierFiltering(py::array_t<double, py::array::c_style> x_array,
                                                 py::array_t<double, py::array::c_style> y_array,
                                                 py::array_t<double, py::array::c_style> z_array,
                                                 const double filtering_constant,
                                                 const double mean_factor,
                                                 const double dev_factor,
                                                 const unsigned int k,
                                                 const bool use_median)
{
  /* Request a buffer descriptor from Python */
  py::buffer_info x_info = x_array.request();
  auto x_coords = static_cast<double *>(x_info.ptr);

  py::buffer_info y_info = y_array.request();
  auto y_coords = static_cast<double *>(y_info.ptr);

  py::buffer_info z_info = z_array.request();
  auto z_coords = static_cast<double *>(z_info.ptr);

  auto result = cars_filter::statistical_filtering(x_coords, y_coords, z_coords, x_info.shape[0], filtering_constant, mean_factor, dev_factor, k, use_median);

  // Copy C++ vector to Python list
  // As result has variable length, I am not sure if the copy can be avoided here in pyBind
  return py::cast(result);
}


py::list pyPointCloudSmallComponentOutlierFiltering(py::array_t<double, py::array::c_style> x_array,
                                                     py::array_t<double, py::array::c_style> y_array,
                                                     py::array_t<double, py::array::c_style> z_array,
                                                     const double radius = 3,
                                                     const unsigned int min_cluster_size = 15,
                                                     const double clusters_distance_threshold = 4)
{
  /* Request a buffer descriptor from Python */
  py::buffer_info x_info = x_array.request();
  auto x_coords = static_cast<double *>(x_info.ptr);

  py::buffer_info y_info = y_array.request();
  auto y_coords = static_cast<double *>(y_info.ptr);

  py::buffer_info z_info = z_array.request();
  auto z_coords = static_cast<double *>(z_info.ptr);

  auto result = cars_filter::point_cloud_small_component_filtering(
                                        x_coords,
                                        y_coords,
                                        z_coords,
                                        x_info.shape[0],
                                        radius,
                                        min_cluster_size,
                                        clusters_distance_threshold
                                        );

  // Copy C++ vector to Python list
  // As result has variable length, I am not sure if the copy can be avoided here in pyBind
  return py::cast(result);
}




// wrap as Python module
PYBIND11_MODULE(outlier_filter, m)
{
  m.doc() = "Wrapper module of cars-filter, a c++ module for 3D point filtering";
  m.def("pc_small_component_outlier_filtering",
        &pyPointCloudSmallComponentOutlierFiltering,
        "Filter outliers from point cloud using small components method",
        py::arg("x_array"),
        py::arg("y_array"),
        py::arg("z_array"),
        py::arg("radius") = 3.,
        py::arg("min_cluster_size") = 50,
        py::arg("clusters_distance_threshold") = std::numeric_limits<double>::quiet_NaN()
  );

  m.def("pc_statistical_outlier_filtering",
        &pyPointCloudStatisticalOutlierFiltering,
        "Filter outliers from point cloud using statistical method",
        py::arg("x_array"),
        py::arg("y_array"),
        py::arg("z_array"),
        py::arg("filtering_constant") = 0.,
        py::arg("mean_factor") = 1.,
        py::arg("dev_factor") = 1.,
        py::arg("k") = 15,
        py::arg("use_median") = false
  );

  m.def("epipolar_small_component_outlier_filtering",
        &pyEpipolarSmallComponentOutlierFiltering,
        "Filter outliers from depth map in epipolar geometry using small components method",
        py::arg("x_values"),
        py::arg("y_values"),
        py::arg("z_values"),
        py::arg("min_cluster_size") = 50,
        py::arg("radius") = 3.,
        py::arg("half_window_size") = 10,
        py::arg("clusters_distance_threshold") = std::numeric_limits<double>::quiet_NaN(),
        py::return_value_policy::take_ownership
  );

  m.def("epipolar_statistical_outlier_filtering",
        &pyEpipolarStatisticalOutlierFiltering,
        "Filter outliers from depth map in epipolar geometry using statistical method",
        py::arg("x_values"),
        py::arg("y_values"),
        py::arg("z_values"),
        py::arg("k") = 15,
        py::arg("half_window_size") = 10,
        py::arg("filtering_constant") = 0.,
        py::arg("mean_factor") = 1.,
        py::arg("dev_factor") = 1.,
        py::arg("use_median") = false,
        py::return_value_policy::take_ownership
  );
}