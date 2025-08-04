///  @file  multibody_element_presetting.cpp
///  @author  Takeshi Chang
///  @date  November 24, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "multibody_presetting.hpp"
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_shell.hpp>
#include <iostream>
#include <misc/assembling.hpp>

namespace icarat
{
namespace multibody
{
using namespace std;
using namespace Eigen;

void ElementMB::makeKeFint(vector<Triplet<double>> &tripletsK,
                           vector<Triplet<double>> &tripletsF,
                           vector<Triplet<double>> &tripletsR, FEM &fem,
                           vector<NodeMB> &node)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      disp.coeffRef(fem.dofnp * i + j) = node[nodeID[i]].val[j];
      idof.coeffRef(fem.dofnp * i + j) = node[nodeID[i]].dof[j];
    }
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = node[nodeID[i]].x[j];
  }

  young_ = SIMP(young1_, young2_, design_s_, pp_);

  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  // gauss loop
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    fe_ += Be.transpose() * De_ * Be * disp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

void ElementMB::makeMisesStress(FEM &fem, vector<NodeMB> &node)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  MatrixXd V = MatrixXd::Zero(fem.voigt, fem.voigt);

  if(fem.ndim == 2)
    V << 1.0, -0.5, 0.0, //
        -0.5, 1.0, 0.0,  //
        0.0, 0.0, 3.0;
  else if(fem.ndim == 3)
    V << 1.0, -0.5, -0.5, 0.0, 0.0, 0.0, //
        -0.5, 1.0, -0.5, 0.0, 0.0, 0.0,  //
        -0.5, -0.5, 1.0, 0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0, 3.0, 0.0, 0.0,    //
        0.0, 0.0, 0.0, 0.0, 3.0, 0.0,    //
        0.0, 0.0, 0.0, 0.0, 0.0, 3.0;

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      disp.coeffRef(fem.ndim * i + j) = node[nodeID[i]].val[j];
      X.coeffRef(i, j) = node[nodeID[i]].x[j];
    }
  }

  strain_.setZero();
  stress_.setZero();

  double fac = (double)(1.0 / ipmax);
  double dummy;
  Bmatrix bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, dummy, eType, ipmax, X, ip);
    strain_ += Be * disp * fac;
    stress_ += De_ * Be * disp * fac;
  }

  mises_ = sqrt(stress_.transpose() * V * stress_);
  // strainenergy_ = sqrt(strain_.transpose() * V * strain_);
  strainenergy_ = 0.5 * stress_.dot(strain_);
}

void ElementMB::setParameter(const toml::value &config)
{
  /// reading the design material
  const auto &mat = toml::find(config, "material");
  assert(toml::find<string>(mat, "type") == "Elastic");
  pp_ = toml::find<double>(mat, "penalty");
  mattype_ = toml::find<string>(mat, "state");
  if(isDesignable == 1)
  {
    auto young = toml::find<array<double, 2>>(mat, "young");
    design_s_ = toml::find<double>(mat, "designs_0");
    young1_ = get<0>(young);
    young2_ = get<1>(young);
    poisson_ = toml::find<double>(mat, "poisson");
  }
  else if(isDesignable == 0)
  {
    /// reading the non-design material
    // using material_type = std::tuple<std::array<double, 1>,
    // std::array<double, 1>,std::array<double, 2>>;
    const auto &non_mat =
        toml::find<toml::array>(config, "non_design", "material");
    design_s_ = toml::get<double>(non_mat[nonID][0]);
    poisson_ = toml::get<double>(non_mat[nonID][1]);
    young1_ = toml::get<double>(non_mat[nonID][2][0]);
    young2_ = toml::get<double>(non_mat[nonID][2][1]);
    // cout << design_s_ << " " << poisson_ << " " << young1_ << " " << young2_
    //      << endl;
  }
  if(toml::find<int>(config, "mesh", "dimension") == 3)
    mattype_ = "3D";
}
} // namespace multibody
} // namespace icarat