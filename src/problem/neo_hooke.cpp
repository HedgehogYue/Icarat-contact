///  @file  neo_hooke.cpp
///  @author  Daiki Watanabe
///  @date  June 4, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "neo_hooke.hpp"
#include <Eigen/LU>
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_total.hpp>
#include <basis/bmatrix_updated.hpp>
#include <iostream>
#include <misc/assembling.hpp>
#include <unsupported/Eigen/KroneckerProduct>
namespace icarat
{
using namespace std;
using namespace Eigen;

void NeoHookeTotal::totalKirchhoff(FEM &fem, VectorXd &sigmaVector, MatrixXd &F,
                                   MatrixXd &sigmaMatrix)
{
  // lame constant for initial & current position
  double mu = 0.5 * young_ / (1 + poisson_);
  double lambda =
      (poisson_ * young_) / ((1.0 + poisson_) * (1.0 - 2.0 * poisson_));

  MatrixXd I = MatrixXd::Identity(fem.ndim, fem.ndim);
  MatrixXd C = MatrixXd::Zero(fem.ndim, fem.ndim);
  MatrixXd E = MatrixXd::Zero(fem.ndim, fem.ndim);

  // right C-G tensor
  C = F.transpose() * F;
  // Green strain
  E = 0.5 * (C - I);

  if(fem.ndim == 2)
  {
    // tensol -> vector */
    VectorXd strainvec = VectorXd::Zero(fem.voigt);
    strainvec[0] = E(0, 0);
    strainvec[1] = E(1, 1);
    strainvec[2] = 2.0 * E(0, 1);

    // cal invers of C
    double jac = F.determinant();
    double jacC = pow(jac, 2.0);
    MatrixXd C_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
    C_inv(0, 0) = C(1, 1) / jacC;
    C_inv(0, 1) = -C(0, 1) / jacC;
    C_inv(1, 0) = -C(1, 0) / jacC;
    C_inv(1, 1) = C(0, 0) / jacC;

    // cal stress
    sigmaVector(0) =
        mu * (1.0 - C_inv(0, 0)) + lambda * jac * (jac - 1.0) * C_inv(0, 0);
    sigmaVector(1) =
        mu * (1.0 - C_inv(1, 1)) + lambda * jac * (jac - 1.0) * C_inv(1, 1);
    sigmaVector(2) =
        mu * (-C_inv(0, 1)) + lambda * jac * (jac - 1.0) * C_inv(0, 1);

    // cal 2nd PK stress S
    MatrixXd sigma = MatrixXd::Zero(fem.ndim, fem.ndim);
    sigma(0, 0) = sigmaVector(0);
    sigma(0, 1) = sigmaVector(2);
    sigma(1, 0) = sigmaVector(2);
    sigma(1, 1) = sigmaVector(1);
    sigmaMatrix = kroneckerProduct(sigma, I);

    // cauchy stress(bonet 5.45b)
    // MatrixXd S1 = F * sigma * F.transpose() / jac;

    // De matrix
    // make IJKL data
    MatrixXi IJdata = MatrixXi::Zero(fem.ndim, fem.ndim);
    IJdata(0, 0) = 0;
    IJdata(0, 1) = 2;
    IJdata(1, 0) = 2;
    IJdata(1, 1) = 1;

    double A = lambda * jac * (2.0 * jac - 1.0);
    double B = 2.0 * (mu - lambda * jac * (jac - 1.0));

    for(int i = 0; i < 2; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        for(int k = 0; k < 2; k++)
        {
          for(int l = 0; l < 2; l++)
          {
            De_(IJdata(i, j), IJdata(k, l)) =
                A * C_inv(i, j) * C_inv(k, l) +
                B * 0.5 *
                    (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k));
          }
        }
      }
    }
  }
  else
  {
    // tensol -> vector */
    VectorXd strainvec = VectorXd::Zero(fem.voigt);
    strainvec[0] = E(0, 0);
    strainvec[1] = E(1, 1);
    strainvec[2] = E(2, 2);
    strainvec[3] = 2.0 * E(0, 1);
    strainvec[4] = 2.0 * E(1, 2);
    strainvec[5] = 2.0 * E(2, 0);
    // cal stress S */
    double jac = F.determinant();
    double jacC = pow(jac, 2.0);
    MatrixXd C_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
    C_inv(0, 0) = (C(1, 1) * C(2, 2) - C(1, 2) * C(2, 1)) / jacC;
    C_inv(0, 1) = (C(0, 2) * C(2, 1) - C(0, 1) * C(2, 2)) / jacC;
    C_inv(0, 2) = (C(0, 1) * C(1, 2) - C(0, 2) * C(1, 1)) / jacC;

    C_inv(1, 0) = (C(1, 2) * C(2, 0) - C(1, 0) * C(2, 2)) / jacC;
    C_inv(1, 1) = (C(0, 0) * C(2, 2) - C(0, 2) * C(2, 0)) / jacC;
    C_inv(1, 2) = (C(0, 2) * C(1, 0) - C(0, 0) * C(1, 2)) / jacC;

    C_inv(2, 0) = (C(1, 0) * C(2, 1) - C(1, 1) * C(2, 0)) / jacC;
    C_inv(2, 1) = (C(0, 1) * C(2, 0) - C(0, 0) * C(2, 1)) / jacC;
    C_inv(2, 2) = (C(0, 0) * C(1, 1) - C(0, 1) * C(1, 0)) / jacC;

    /*--------------------------------------------- cal stress */
    sigmaVector(0) =
        mu * (1.0 - C_inv(0, 0)) + lambda * jac * (jac - 1.0) * C_inv(0, 0);
    sigmaVector(1) =
        mu * (1.0 - C_inv(1, 1)) + lambda * jac * (jac - 1.0) * C_inv(1, 1);
    sigmaVector(2) =
        mu * (1.0 - C_inv(2, 2)) + lambda * jac * (jac - 1.0) * C_inv(2, 2);
    sigmaVector(3) =
        mu * (-C_inv(0, 1)) + lambda * jac * (jac - 1.0) * C_inv(0, 1);
    sigmaVector(4) =
        mu * (-C_inv(1, 2)) + lambda * jac * (jac - 1.0) * C_inv(1, 2);
    sigmaVector(5) =
        mu * (-C_inv(2, 0)) + lambda * jac * (jac - 1.0) * C_inv(2, 0);

    MatrixXd sigma = MatrixXd::Zero(fem.ndim, fem.ndim);
    sigma(0, 0) = sigmaVector[0];
    sigma(1, 1) = sigmaVector[1];
    sigma(2, 2) = sigmaVector[2];
    sigma(0, 1) = sigmaVector[3];
    sigma(1, 2) = sigmaVector[4];
    sigma(2, 0) = sigmaVector[5];
    sigma(1, 0) = sigma(0, 1);
    sigma(2, 1) = sigma(1, 2);
    sigma(0, 2) = sigma(2, 0);
    sigmaMatrix.setZero();
    sigmaMatrix.block(0, 0, 3, 3) = sigma;
    sigmaMatrix.block(3, 3, 3, 3) = sigma;
    sigmaMatrix.block(6, 6, 3, 3) = sigma;

    // De matrix
    // make IJKL data
    MatrixXi IJdata = MatrixXi::Zero(fem.ndim, fem.ndim);
    IJdata(0, 0) = 0;
    IJdata(0, 1) = 3;
    IJdata(0, 2) = 5;

    IJdata(1, 0) = 3;
    IJdata(1, 1) = 1;
    IJdata(1, 2) = 4;

    IJdata(2, 0) = 5;
    IJdata(2, 1) = 4;
    IJdata(2, 2) = 2;

    double A = lambda * jac * (2.0 * jac - 1.0);
    double B = 2.0 * (mu - lambda * jac * (jac - 1.0));

    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        for(int k = 0; k < 3; k++)
        {
          for(int l = 0; l < 3; l++)
          {
            De_(IJdata(i, j), IJdata(k, l)) =
                A * C_inv(i, j) * C_inv(k, l) +
                B * 0.5 *
                    (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k));
          }
        }
      }
    }
  }
}

void FiniteNeoHookeTotal::makeKeFint(vector<Triplet<double>> &tripletsK,
                                     vector<Triplet<double>> &tripletsF,
                                     vector<Triplet<double>> &tripletsR,
                                     FEM &fem, vector<Node> &nodes)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = gamma_ * nodes[nodeID[i]].val[j];
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].dof[j];
    }
  }

  volume = 0.0;
  Ke_.setZero();
  fe_.setZero();

  // linear init
  l_.init();
  Bmatrix bLinear;

  // nonlinear init
  MatrixXd nKe = Ke_;
  VectorXd nfe = fe_;
  BmatrixTotal bTotal;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    young_ = SIMP(young1_, young2_, design_s_, pp_);

    // linear part
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
    De(l_.De, young_, poisson_);
    bLinear.make(Be, Ne, jac, eType, ipmax, X, ip);
    l_.strain = Be * disp;
    l_.stress = l_.De * l_.strain;
    l_.fe += Be.transpose() * l_.stress * jac;
    l_.Ke += Be.transpose() * l_.De * Be * jac;
    volume += jac;

    // nonlinear part
    jac = 1.0;
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, fem.ndim, ne, ipmax, X, U, ip);
    assert(jac > 0.0);

    totalKirchhoff(fem, sigmaVector, F, Sbar);

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * De_ * BeL * jac +
           BeNL.transpose() * Sbar * BeNL * jac;
  }

  // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
  Ke_ = (1.0 - gamma_ * gamma_) * l_.Ke + gamma_ * gamma_ * nKe;

  // inner stress(Bonet's book (9.15b))
  fe_ = (1.0 - gamma_ * gamma_) * l_.fe + gamma_ * nfe;

  lfe_ = l_.fe;
  nfe_ = nfe;

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

void mHomoNeoHookeTotal::makeKeHomo(vector<Triplet<double>> &tripletsK,
                                    vector<Triplet<double>> &tripletsF,
                                    vector<Triplet<double>> &tripletsR,
                                    FEM &fem, vector<HomoNode> &nodes, int dir)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = gamma_ * nodes[nodeID[i]].nmtval[dir][j];
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].nmtval[dir][j];
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
    }
  }

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;

  // linear init
  l_.init();
  Bmatrix bLinear;

  // nonlinear init
  MatrixXd nKe = Ke_;
  VectorXd nfe = fe_;
  BmatrixTotal bTotal;
  for(int ip = 0; ip < ipmax; ip++)
  {
    young_ = SIMP(young1_, young2_, design_s_, pp_);

    // linear part
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
    De(l_.De, young_, poisson_);
    bLinear.make(Be, Ne, jac, eType, ipmax, X, ip);
    l_.strain = Be * disp;
    l_.stress = l_.De * l_.strain;
    l_.fe += Be.transpose() * l_.stress * jac;
    l_.Ke += Be.transpose() * l_.De * Be * jac;
    volume += jac;

    // nonlinear part
    jac = 1.0;
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix
    bTotal.make(BeL, BeNL, Ne, jac, F, fem.ndim, ne, ipmax, X, U, ip);
    assert(jac > 0.0);

    totalKirchhoff(fem, sigmaVector, F, Sbar);

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * De_ * BeL * jac +
           BeNL.transpose() * Sbar * BeNL * jac;
  }

  // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
  Ke_ = (1.0 - gamma_ * gamma_) * l_.Ke + gamma_ * gamma_ * nKe;

  // inner stress(Bonet's book (9.15b))
  fe_ = (1.0 - gamma_ * gamma_) * l_.fe + gamma_ * nfe;

  assemblingHomo(ne, numdof, idof, Ke_, fe_, tripletsK, tripletsF, fem, nodeID,
                 nodes);
  assemblingReactHomo(ne, numdof, idof, fe_, tripletsR, fem, nodeID, nodes);
}

void mHomoNeoHookeTotal::makeMicroValues(FEM &fem, vector<HomoNode> &nodes)
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