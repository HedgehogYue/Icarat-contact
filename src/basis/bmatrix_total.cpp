///  @file  bmatrix_total.cpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "bmatrix_total.hpp"
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/LU>
#include <iostream>
#include <problem/base.hpp>

namespace icarat
{
void BmatrixTotal::make(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                        Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &F,
                        int ndim, int ne, int ipmax, Eigen::MatrixXd &X,
                        Eigen::MatrixXd &U, int ip)
{
  // init
  BL.setZero();
  BNL.setZero();
  Ne.setZero();
  jac = 0.0;
  F.setZero();

  if(ndim == 2)
  {
    if(ne == 3)
      bmatrix2dTotal<ShapeFunction3Triangle, Gauss1Triangle>(BL, BNL, Ne, jac,
                                                             F, X, U, ip);
    else if(ne == 6)
      bmatrix2dTotal<ShapeFunction6Triangle, Gauss3Triangle>(BL, BNL, Ne, jac,
                                                             F, X, U, ip);
    else if(ne == 4)
    {
      if(ipmax == 1)
        bmatrix2dTotal<ShapeFunction4Square, Gauss1Square>(BL, BNL, Ne, jac, F,
                                                           X, U, ip);
      else if(ipmax == 4)
        bmatrix2dTotal<ShapeFunction4Square, Gauss4Square>(BL, BNL, Ne, jac, F,
                                                           X, U, ip);
      else
      {
        std::cerr << "Ipmax is invalid in bmatrix. " << std::endl;
        exit(1);
      }
    }
    else if(ne == 8)
      bmatrix2dTotal<ShapeFunction8Square, Gauss9Square>(BL, BNL, Ne, jac, F, X,
                                                         U, ip);
    else
    {
      std::cerr << "The combination of ne and ndim is invalid in bmatrix. "
                << std::endl;
      exit(1);
    }
  }
  else if(ndim == 3)
  {
    if(ne == 4)
      bmatrix3dTotal<ShapeFunction4Tetrahedron, Gauss1Tetrahedron>(
          BL, BNL, Ne, jac, F, X, U, ip);
    else if(ne == 8)
      bmatrix3dTotal<ShapeFunction8Cubic, Gauss8Cubic>(BL, BNL, Ne, jac, F, X,
                                                       U, ip);
    else if(ne == 20)
      bmatrix3dTotal<ShapeFunction20Cubic, Gauss27Cubic>(BL, BNL, Ne, jac, F, X,
                                                         U, ip);
    else
    {
      std::cerr << "The combination of ne and ndim is invalid in bmatrix. "
                << std::endl;
      exit(1);
    }
  }
}

} // namespace icarat
