///  @file  eigen_manipulation.cpp
///  @author  Daiki Watanabe
///  @date  November 20, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once

#include <Eigen/Core>
#include <problem/base.hpp>

namespace icarat
{
/// This is a reference to numpy.
inline Eigen::MatrixXd hstack(Eigen::MatrixXd const &a,
                              Eigen::MatrixXd const &b)
{
  Eigen::MatrixXd c(a.rows(), a.cols() + b.cols());
  c << a, b;
  return c;
}

/// This is a reference to numpy.
inline Eigen::MatrixXd vstack(Eigen::MatrixXd const &a,
                              Eigen::MatrixXd const &b)
{
  Eigen::MatrixXd c(a.rows() + b.rows(), a.cols());
  c << a, b;
  return c;
}

} // namespace icarat