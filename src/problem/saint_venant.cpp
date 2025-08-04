///  @file  saint_venant.cpp
///  @author  Daiki Watanabe
///  @date  June 4, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "saint_venant.hpp"
#include <Eigen/LU>
#include <basis/bmatrix_total.hpp>
#include <basis/bmatrix_updated.hpp>
#include <misc/assembling.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace icarat
{
using namespace std;
using namespace Eigen;

void SaintVenantTotal::totalKirchhoff(FEM &fem, VectorXd &sigmaVector,
                                      MatrixXd &F, MatrixXd &sigmaMatrix)
{
  MatrixXd I = MatrixXd::Identity(fem.ndim, fem.ndim);

  // Cauchy Greeen tensor E
  MatrixXd E = 0.5 * (F.transpose() * F - I);

  if(fem.ndim == 2)
  {
    // convert to voigt expression
    VectorXd Ev(fem.voigt);
    Ev << E(0, 0), E(1, 1), E(0, 1) + E(1, 0);

    // Create elasticity tensor.
    double mu0 = 0.5 * young_ / (1 + poisson_);
    double lambda0 = 2.0 * mu0 * poisson_ / (1.0 - 2.0 * poisson_);

    lameDe(De_, mu0, lambda0);

    // Saint-Venent Kirchhoff
    sigmaVector = De_ * Ev;

    // make 9*9 matrix (refer bonet book(5.43))
    MatrixXd sigma(fem.ndim, fem.ndim);

    sigma << sigmaVector(0), sigmaVector(2), //
        sigmaVector(2), sigmaVector(1);

    // matrix expression
    sigmaMatrix = kroneckerProduct(sigma, I);
  }
  else
  {
    // convert to voigt expression
    VectorXd Ev(fem.voigt);
    Ev << E(0, 0), E(1, 1), E(2, 2), E(0, 1) + E(1, 0), E(1, 2) + E(2, 1),
        E(2, 0) + E(0, 2);

    // Create elasticity tensor.
    young_ = SIMP(young1_, young2_, design_s_, pp_);
    double mu0 = 0.5 * young_ / (1 + poisson_);
    double lambda0 = 2.0 * mu0 * poisson_ / (1.0 - 2.0 * poisson_);

    lameDe(De_, mu0, lambda0);

    // Saint-Venent Kirchhoff
    sigmaVector = De_ * Ev;

    // make 9*9 matrix (refer bonet book(5.43))
    MatrixXd sigma(fem.ndim, fem.ndim);

    sigma << sigmaVector(0), sigmaVector(3), sigmaVector(5), //
        sigmaVector(3), sigmaVector(1), sigmaVector(4),      //
        sigmaVector(5), sigmaVector(4), sigmaVector(2);

    // matrix expression
    sigmaMatrix = kroneckerProduct(sigma, I);
  }
}

void FiniteSaintVenantTotal::makeKeFint(vector<Triplet<double>> &tripletsK,
                                        vector<Triplet<double>> &tripletsF,
                                        vector<Triplet<double>> &tripletsR,
                                        FEM &fem, vector<Node> &nodes)
{
  /// initial coordinate
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  /// displacement in this increment
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  /// array of DOF
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i * j) = nodes[nodeID[i]].dof[j];
    }
  }

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  BmatrixTotal bmatrix;
  /// gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    /// second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bmatrix.make(BeL, BeNL, Ne, jac, F, fem.ndim, ne, ipmax, X, U, ip);

    young_ = SIMP(young1_, young2_, design_s_, pp_);
    totalKirchhoff(fem, sigmaVector, F, Sbar);

    /// inner stress
    fe_ += BeL.transpose() * sigmaVector * jac;

    /// element stiffness matrix
    /// refer bonet's book formula (8.10)
    Ke_ += BeL.transpose() * De_ * BeL * jac;
    Ke_ += BeNL.transpose() * Sbar * BeNL * jac;

    volume += jac;
  } // gauss point loop end

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

void mHomoSaintVenantTotal::makeKeHomo(vector<Triplet<double>> &tripletsK,
                                       vector<Triplet<double>> &tripletsF,
                                       FEM &fem, vector<HomoNode> &nodes,
                                       int dir)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
      U.coeffRef(i, j) = nodes[nodeID[i]].nmtval[dir][j];
    }

    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  young_ = SIMP(young1_, young2_, design_s_, pp_);
  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  BmatrixTotal bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    // stess for internal force
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bmatrix.make(BeL, BeNL, Ne, jac, F, fem.ndim, ne, ipmax, X, U, ip);

    totalKirchhoff(fem, sigmaVector, F, Sbar);

    // inner stress(Bonet's book (9.15b))
    fe_ += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    Ke_ += BeL.transpose() * De_ * BeL * jac +
           BeNL.transpose() * Sbar * BeNL * jac;

    volume += jac;
  }

  assemblingHomo(ne, numdof, idof, Ke_, fe_, tripletsK, tripletsF, fem, nodeID,
                 nodes);
}

void mHomoSaintVenantTotal::makeMicroValues(FEM &fem, vector<HomoNode> &nodes)
{
  VectorXi idof = VectorXi::Zero(numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd disps = MatrixXd::Zero(numdof, fem.voigt); ///< for nmt

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      idof.coeffRef(i * fem.dofnp + j) = nodes[nodeID[i]].dof[j];

      for(int k = 0; k < fem.voigt; k++)
        disps(fem.dofnp * i + j, k) = nodes[nodeID[i]].nmtval[k][j];
    }
  }

  sEnergy_.setZero();
  strain_.setZero();
  stress_.setZero();
  volume = 0.0;
  double fac = (double)(1.0 / ipmax);
  BmatrixTotal bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    // stess for internal force
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    bmatrix.make(BeL, BeNL, Ne, jac, F, fem.ndim, ne, ipmax, X, U, ip);

    totalKirchhoff(fem, sigmaVector, F, Sbar);

    MatrixXd mstrain = BeL * disps;
    volume += jac;
    // strain_ += BeL * disp * fac;
    //  stress_ += De_ * Be * disp * fac;

    sEnergy_ += mstrain.transpose() * De_ * mstrain * jac;
  }

  // make mises stress
  // MatrixXd V = MatrixXd::Zero(fem.voigt, fem.voigt);
  // if(fem.ndim == 2)
  //   V << 1.0, -0.5, 0.0, //
  //       -0.5, 1.0, 0.0,  //
  //       0.0, 0.0, 3.0;
  // else if(fem.ndim == 3)
  //   V << 1.0, -0.5, -0.5, 0.0, 0.0, 0.0, //
  //       -0.5, 1.0, -0.5, 0.0, 0.0, 0.0,  //
  //       -0.5, -0.5, 1.0, 0.0, 0.0, 0.0,  //
  //       0.0, 0.0, 0.0, 3.0, 0.0, 0.0,    //
  //       0.0, 0.0, 0.0, 0.0, 3.0, 0.0,    //
  //       0.0, 0.0, 0.0, 0.0, 0.0, 3.0;

  // mises_ = sqrt(stress_.transpose() * V * stress_);
}

} // namespace icarat