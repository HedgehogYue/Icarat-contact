///  @file  neo_hooke.hpp
///  @author  Daiki Watanabe
///  @date  June 4, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "elasticity.hpp"

namespace icarat
{
/// compressible Neo-Hookean model of total lagrange
class NeoHookeTotal : public Elastic
{
public:
  NeoHookeTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_)
  {
  }

  void totalKirchhoff(FEM &fem, Eigen::VectorXd &sigmaVector,
                      Eigen::MatrixXd &F, Eigen::MatrixXd &Sbar);
};

/// stabilizer class for numerical instability in the finite deformation's TO
///  @sa Wang G. et al.: Interpolation scheme for fictitious domain techniques
///  and topology optimization of finite strain elastic problems (2014).
struct LinearStabilizer
{
  LinearStabilizer(int voigt_, int numdof_)
      : strain(voigt_), stress(voigt_), fe(numdof_), De(voigt_, voigt_),
        Ke(numdof_, numdof_)
  {
  }

  void init()
  {
    strain.setZero();
    stress.setZero();
    fe.setZero();
    De.setZero();
    Ke.setZero();
  }

  /// linear parts
  Eigen::VectorXd strain;
  Eigen::VectorXd stress;
  Eigen::VectorXd fe;
  Eigen::MatrixXd De;
  Eigen::MatrixXd Ke;
};

/// compressible Neo-Hookean model for finite deformation
class FiniteNeoHookeTotal : public NeoHookeTotal
{
public:
  FiniteNeoHookeTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : NeoHookeTotal(voigt_, ne_, ipmax_, numdof_), gamma_(0.0),
        l_(voigt_, numdof_), lfe_(numdof_), nfe_(numdof_)
  {
  }

  void setGamma(double gamma) { gamma_ = gamma; }

  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);

  const double &gamma() const { return gamma_; }
  const Eigen::VectorXd &lfe() const { return lfe_; }
  const Eigen::VectorXd &nfe() const { return nfe_; }

protected:
  double gamma_;        ///< threshold parameter
  LinearStabilizer l_;  ///< fictitious linear class
  Eigen::VectorXd lfe_; ///< linear inner force
  Eigen::VectorXd nfe_; ///< nonlinear inner force
};

/// micro scale elastic model class based on homogenization theory
///@sa  寺田賢二郎ら，「数値材料試験」，丸善出版，2021.
class mHomoNeoHookeTotal : public NeoHookeTotal
{
public:
  mHomoNeoHookeTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : NeoHookeTotal(voigt_, ne_, ipmax_, numdof_), gamma_(0.0),
        l_(voigt_, numdof_), sEnergy_(voigt_, voigt_)
  {
  }

  void setGamma(double gamma) { gamma_ = gamma; }

  /// used for homogenization in micro analysis
  void makeKeHomo(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<HomoNode> &nodes, int dir);

  /// generate the homogenized elasticity tensor in micro analysis
  void makeMicroValues(FEM &fem, std::vector<HomoNode> &nodes);

  /// getter
  const double &gamma() const { return gamma_; }
  const Eigen::MatrixXd &sEnergy() const { return sEnergy_; }

protected:
  double gamma_;            ///< threshold parameter
  LinearStabilizer l_;      ///< fictitious linear class
  Eigen::MatrixXd sEnergy_; ///< strain energy
};

} // namespace icarat