///  @file  multibody_presetting.hpp
///  @author  Takeshi Chang
///  @date  November 24, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <problem/elasticity.hpp>

namespace icarat
{
namespace multibody
{

struct NodeMB : Node
{
  /// w=0 design domain do not care about the strain energy
  /// w=1 non-design domain need to calculate the strain energy
  double w; /// distinguish design domain
  // double wdof[3];

  double isConstraint; /// 0 is not constraint, 1 is constraint
  double constvalue;   /// the constraint value of this node

  int nwd; /// distinguish dirichlet condition
           /// nwd=0 two degree of freedom is fixed
           /// nwd=1 x degree of freedom is fixed
           /// nwd=-1 y degree of freedom is fixed
           /// nwd=2 no dirichlet condition
};

class ElementMB : public LinearElastic
{
public:
  int nonID;         // the ID of the nondesign domain
  bool isDesignable; /// 0 is non-design domain, 1 is design domain
  const double &strainenergy() const noexcept { return strainenergy_; }
  double strainenergy_;

public:
  ElementMB(int voigt_, int ne_, int ipmax_, int numdof_)
      : LinearElastic(voigt_, ne_, ipmax_, numdof_)
  {
  }
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<NodeMB> &node);

  /// make mises stress in this elmenet (mises_)
  void makeMisesStress(FEM &fem, std::vector<NodeMB> &node);

  /// set De parameters to this element form input file
  void setParameter(const toml::value &config);
};

} // namespace multibody
} // namespace icarat