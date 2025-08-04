///  @file  saint_venant.hpp
///  @author  Daiki Watanabe
///  @date  June 4, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "elasticity.hpp"

namespace icarat
{
/// Saint-Venant material model class in Total Lagrangian Method
class SaintVenantTotal : public LinearElastic
{
public:
  SaintVenantTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : LinearElastic(voigt_, ne_, ipmax_, numdof_)
  {
  }

  void totalKirchhoff(FEM &fem, Eigen::VectorXd &sigmaVector,
                      Eigen::MatrixXd &F, Eigen::MatrixXd &sigmaMatrix);
};

class FiniteSaintVenantTotal : public SaintVenantTotal
{
public:
  FiniteSaintVenantTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : SaintVenantTotal(voigt_, ne_, ipmax_, numdof_)
  {
  }

  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};

/// micro scale elastic model class based on homogenization theory
///@sa  寺田賢二郎ら，「数値材料試験」，丸善出版，2021.
class mHomoSaintVenantTotal : public SaintVenantTotal
{
public:
  mHomoSaintVenantTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : SaintVenantTotal(voigt_, ne_, ipmax_, numdof_), sEnergy_(voigt_, voigt_)
  {
  }

  /// used for homogenization in micro analysis
  void makeKeHomo(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF, FEM &fem,
                  std::vector<HomoNode> &nodes, int dir);

  /// generate the homogenized elasticity tensor in micro analysis
  void makeMicroValues(FEM &fem, std::vector<HomoNode> &nodes);

  /// getter
  const Eigen::MatrixXd &sEnergy() const { return sEnergy_; }

protected:
  Eigen::MatrixXd sEnergy_; ///< strain energy
};

} // namespace icarat