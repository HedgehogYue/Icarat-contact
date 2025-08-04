///  @file  newton_raphson.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include "analysis_utility.hpp"
#include "solver/solvers.hpp"
#include <problem/base.hpp>

namespace icarat
{
/// FEM procedure for linear elasticity
template <class E, class N> class ArcLength
{
public:
  ArcLength(FEM &fem, std::vector<E> &element, std::vector<N> &node,
            Force &force, Value &value);

  ///  main function (Incremental load or displacement method)
  void solve(int thisStep, int totalStep, std::string solvName, bool log);

  /// initializaiton
  void initialization();

protected:
  void divideValIntoNode(int dofnp, Eigen::VectorXd &val, Eigen::VectorXd &du,
                         int thisStep, int totalStep);

  /// assembling stiffness matrix
  void assemblingK(bool log);

  double getArcLength(Eigen::VectorXd &incre, Eigen::VectorXd &incre_lin);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  Force &force;
  Value &value;
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::VectorXd reaction;            ///< reaction force
  Eigen::VectorXd disp1;               ///< disp increment in a previous step
  Eigen::VectorXd disp2;               ///< disp increment in a previous step
  int nitar;                           ///< iteration in previous step
  double ll;                           ///< length in previous step
  double l_max;
  double length;
  double lambda;
  double lambda1;
  double lambda2;
};

template <class E, class N>
ArcLength<E, N>::ArcLength(FEM &fem, std::vector<E> &element,
                           std::vector<N> &node, Force &force, Value &value)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq), disp1(fem.numeq),
      disp2(fem.numeq), nitar(0), ll(0.0), l_max(0.0), length(0.0), lambda(0.0),
      lambda1(0.0), lambda2(0.0)
{
}

template <class E, class N>
void ArcLength<E, N>::solve(int thisStep, int totalStep, std::string solvName,
                            bool log)
{

  Eigen::VectorXd incre_lin = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd incre = Eigen::VectorXd::Zero(fem.numeq);
  disp1.setZero();
  disp2.setZero();
  Solvers solver0(fem, globalK, force.fext, incre_lin);
  Solvers solver1(fem, globalK, force.fres, incre);

  /// set dirichlet boundary condition
  // for(auto &n : node)
  // {
  //   for(int j = 0; j < 3; j++)
  //   {
  //     n.val[j] = thisStep * n.dirich->val[j] / totalStep;
  //     n.du[j] = n.dirich->val[j] / totalStep;
  //   }
  // }

  // assembling
  assemblingK(log);

  // cal residual force
  solver0.solve(solvName, log);

  double aa = incre_lin.dot(incre_lin);

  double ll, factor = 8.0, dlambda0 = 1.0e-3, pp = 5.0;
  if(thisStep == 1)
  {
    ll = sqrt(pow(dlambda0, 2.0) * aa);
    l_max = pp * ll;
  }
  else
  {
    ll = length * factor / nitar;
    if(ll > l_max)
      ll = l_max;
  }

  double sign = incre.dot(incre_lin);
  if(sign >= 0.0)
    sign = 1.0;
  else
    sign = -1.0;

  if(aa <= 0.0)
  {
    std::cerr << "error in arc_length.hpp" << std::endl;
    exit(1);
  }

  double dlambda = sign * fabs(ll / sqrt(aa));
  lambda1 = lambda1 + dlambda;
  lambda2 = lambda + lambda1;

  disp1 = dlambda * incre_lin;
  disp2 = value.val + disp1;

  // make result to node
  divideValIntoNode(fem.dofnp, disp2, disp1, thisStep, totalStep);

  double fnorm = lambda2 * sqrt(force.fext.dot(force.fext));

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter "
            << " first res "
            << "      res " << std::endl;

  int iterstep = 0;
  double residual = 1.0;
  while(residual > 1.0e-4)
  {
    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = lambda2 * force.fext - force.fint;

    solver0.solve(solvName, log);
    solver1.solve(solvName, log);

    dlambda = getArcLength(incre, incre_lin);

    // update length parameters
    lambda1 = lambda1 + dlambda;
    lambda2 = lambda + lambda1;

    incre += dlambda * incre_lin;

    disp1 += incre;
    disp2 = value.val + disp1;
    value.val = disp2;

    // make result to node
    divideValIntoNode(fem.dofnp, disp2, disp1, thisStep, totalStep);
    divideReactionIntoNode(fem.dofnp, reaction, node);

    /// convergence check
    double resnorm = sqrt(force.fres.dot(force.fres));
    residual = resnorm / fnorm;

    std::cout << "   " << iterstep << "     " << fnorm << "       " << residual
              << std::endl;

    if(iterstep == 100)
    {
      std::cerr << "itaration step > 100 in newton_raphson.hpp" << std::endl;
      exit(1);
      return;
    }
  } ///< NR loop end

  // store value
  length = ll;
  nitar = iterstep;

  value.val = disp2;
  // make result to node
  divideValIntoNode(fem.dofnp, disp2, disp1, thisStep, totalStep);
  divideReactionIntoNode(fem.dofnp, reaction, node);
}

template <class E, class N> void ArcLength<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
    {
      n.val[j] = 0.0;
      n.du[j] = 0.0;
    }

  /// ValueNR class
  value.val.setZero();
}

template <class E, class N>
void ArcLength<E, N>::divideValIntoNode(int dofnp, Eigen::VectorXd &val,
                                        Eigen::VectorXd &du, int thisStep,
                                        int totalStep)
{
  int dofnptmp = std::min(dofnp, 3);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
  {
    for(int j = 0; j < dofnptmp; j++)
    {
      if(node[i].dof[j] < 0)
      {
        node[i].val[j] = thisStep * node[i].dirich->val[j] / totalStep;
        node[i].du[j] = node[i].dirich->val[j] / totalStep;
      }
      else
      {
        node[i].val[j] = val.coeff(node[i].dof[j]);
        node[i].du[j] = du.coeff(node[i].dof[j]);
      }
    }
  }
}

template <class E, class N> void ArcLength<E, N>::assemblingK(bool log)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  reaction.setZero();
  force.fint.setZero();
  globalK.setZero();
  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> reactMat(fem.neq, fem.neq);
  std::vector<Eigen::Triplet<double>> tripletsK, tripletsF, tripletsR;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripK_pri, tripF_pri, tripR_pri;

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].makeKeFint(tripK_pri, tripF_pri, tripR_pri, fem, node);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      /// connect
      tripletsK.insert(tripletsK.end(), tripK_pri.begin(), tripK_pri.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR.insert(tripletsR.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalK.setFromTriplets(tripletsK.begin(), tripletsK.end());
  fintMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  reactMat.setFromTriplets(tripletsR.begin(), tripletsR.end());

  reaction = reactMat.diagonal();
  force.fint = fintMat.diagonal();
  globalK.makeCompressed();

  if(log == 1)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

template <class E, class N>
double ArcLength<E, N>::getArcLength(Eigen::VectorXd &incre,
                                     Eigen::VectorXd &incre_lin)
{
  double dlambda;
  double tol = 1.0e-8;

  double aa = 0.0, bb = 0.0, cc = 0.0;
  for(int i = 0; i < this->fem.numeq; i++)
  {
    aa += incre_lin.coeff(i) * incre_lin.coeff(i);
    bb += 2.0 * (disp1.coeff(i) + incre.coeff(i)) * incre_lin.coeff(i);
    cc += pow(disp1.coeff(i) + incre.coeff(i), 2.0) - pow(ll, 2.0);
  }

  double val0 = pow(bb, 2.0) - 4.0 * aa * cc;

  /// cal dlambda
  std::vector<double> dl(2);
  if(val0 > tol)
  {
    if(bb >= 0.0)
    {
      dl[0] = (-bb - sqrt(val0)) / (2.0 * aa);
      dl[1] = (cc / aa) / dl[0];
    }
    else
    {
      dl[0] = (-bb + sqrt(val0)) / (2.0 * aa);
      dl[1] = (cc / aa) / dl[0];
    }

    double val1 = (disp1 + incre + dl[0] * incre_lin).dot(disp1);
    double val2 = (disp1 + incre + dl[1] * incre_lin).dot(disp1);

    if(val1 > val2)
      dlambda = dl[0];
    else
      dlambda = dl[1];
  }
  else if(fabs(val0) <= tol)
  {
    dlambda = -bb / (2.0 * aa);
  }
  else // if(val0 < tol)
  {
    std::cerr << "arc length parameter is negative in newton_raphson.hpp"
              << std::endl;
    exit(1);
  }

  return dlambda;
}
} // namespace icarat
