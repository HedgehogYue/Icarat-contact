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
template <class E, class N> class NewtonRaphson
{
public:
  NewtonRaphson(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                Force &force, Value &value);

  ///  main function (Incremental load or displacement method)
  ///@return Newton-Raphson iteration
  int solve(int thisTime, int totalTime, double tolNR, int maxIter,
            std::string solvName, int log);

  /// main function (arc length method)
  void solveArcLength(int thisStep, int totalStep, std::string solvName,
                      int log);

  /// initializaiton
  void initialization();

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }

private:
  void assemblingK(int log);

  void divideValIntoNode(int dofnp, Eigen::VectorXd &val, Eigen::VectorXd &du,
                         int thisStep, int totalStep);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  Force &force;
  Value &value;
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::VectorXd reaction;            ///< reaction force
};

template <class E, class N>
NewtonRaphson<E, N>::NewtonRaphson(FEM &fem, std::vector<E> &element,
                                   std::vector<N> &node, Force &force,
                                   Value &value)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq)
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;
}

template <class E, class N>
int NewtonRaphson<E, N>::solve(int thisTime, int totalTime, double tolNR,
                               int maxIter, std::string solvName, int log)
{

  int iterstep = 0;
  double firstNorm, residual = 1;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);

  // set dirichlet boundary condition
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
      for(int j = 0; j < 3; j++)
      {
        if(n.dirich->flag[j] == 1)
        {
          n.val[j] = thisTime * n.dirich->val[j] / totalTime;
          n.du[j] = n.dirich->val[j] / totalTime;
        }
        else
          n.du[j] = 0.0;
      }
    else
      for(int j = 0; j < 3; j++)
      {
        n.du[j] = 0.0;
      }
  }

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter "
            << " first norm "
            << "      res " << std::endl;

  while(residual > tolNR)
  {
    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    double norm = sqrt(force.fres.dot(force.fres));
    if(iterstep == 1)
      firstNorm = norm;

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;

    if(iterstep == maxIter)
    {
      std::cerr << "itaration step reach " << maxIter
                << " in newton_raphson.hpp" << std::endl;
      return iterstep;
    }

    sol.setZero();
    Solvers solv(fem, globalK, force.fres, sol);
    solv.solve(solvName, log);

    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
    divideReactionIntoNode(fem.dofnp, reaction, node);
  } ///< NR loop end

  // make result to node
  divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  return iterstep;
}

template <class E, class N> void NewtonRaphson<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
    {
      n.val[j] = 0.0;
      n.du[j] = 0.0;
      n.reaction[j] = 0.0;
    }

  /// ValueNR class
  value.val.setZero();
}

template <class E, class N> void NewtonRaphson<E, N>::assemblingK(int log)
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
void NewtonRaphson<E, N>::divideValIntoNode(int dofnp, Eigen::VectorXd &val,
                                            Eigen::VectorXd &du, int thisStep,
                                            int totalStep)
{
  int dofnptmp = std::min(3, dofnp);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
  {
    for(int j = 0; j < dofnptmp; j++)
    {
      if(node[i].dof[j] >= fem.numeq)
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

} // namespace icarat
