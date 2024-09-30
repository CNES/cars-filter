#!/usr/bin/env python
# coding: utf8
#
# Copyright (C) 2024 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of cars-filter
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""Tests for `cars-filter` package."""

import outlier_filter


def test_synthetic():
    """
    Test outlier filtering on small synthetic data
    """
    synthetic_array = [
        [10.2, 10.5, 10.4, 10.8, 10.9, 10.2, 10.3, 10.1],
        [20.2, 20.5, 20.4, 20.8, 20.9, 20.2, 20.3, 20.1],
        [30.2, 30.5, 30.4, 30.8, 30.9, 30.2, 30.3, 30.1],
    ]
    outlier_filter.pc_outlier_filtering(synthetic_array)
