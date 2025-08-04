///  @file  multibody_linear.hpp
///  @author  Takeshi Chang
///  @date  May 26, 2023.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include "analysis_utility.hpp"
#include "solver/solvers.hpp"
#include <problem/base.hpp>

namespace icarat
{

// using namespace std;
// using namespace Eigen;
template <class E, class N> class LinearAnalysisMB
{
private:
  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  Force &force;
  Value &value;

  Eigen::SparseMatrix<double>
      globalKnon; /// stiffness matrix of the non-design domain
  Eigen::SparseMatrix<double>
      globalKdesign;                   /// stiffness matrix of the design domain
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::VectorXd reaction;            ///< reaction force
public:
  /// main function
  void solve(std::string solvName, bool log);

  ///   initializaiton
  void initialization();

  /// assembling stiffness matrix & inner force
  void assemblingK(bool log);
  void assemblingKnon(bool log);
  void assemblingKdesign(bool log);

  /// getters
  const Eigen::SparseMatrix<double> &Knon() const { return globalKnon; }
  const Eigen::SparseMatrix<double> &Kdesign() const { return globalKdesign; }
  const Eigen::SparseMatrix<double> &K() const { return globalK; }

  LinearAnalysisMB(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                   Force &force, Value &value);
};

template <class E, class N>
LinearAnalysisMB<E, N>::LinearAnalysisMB(FEM &fem, std::vector<E> &element,
                                         std::vector<N> &node, Force &force,
                                         Value &value)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq),
      globalKnon(fem.numeq, fem.numeq), globalKdesign(fem.numeq, fem.numeq)
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;
}

template <class E, class N>
void LinearAnalysisMB<E, N>::solve(std::string solvName, bool log)
{
  // init
  initialization();

  std::cout << "--------------- linear analysis ---------------" << std::endl;

  // set dirichlet boundary condition
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
      for(int j = 0; j < 3; j++)
      {
        if(n.dirich->flag[j] == 1)
          n.val[j] = n.dirich->val[j];
      }
  }

  // cal external force at 1 step
  force.fext = force.fext_org;

  // assembling
  assemblingK(log);
  assemblingKnon(log);
  assemblingKdesign(log);

  // cal residual force
  force.fres = force.fext - force.fint;
  double fnorm = sqrt(force.fres.dot(force.fres));

  // solve KU=F with
  Solvers solv(fem, globalK, force.fres, value.val);
  solv.solve(solvName, log);

  divideValIntoNode(fem.dofnp, fem.numeq, value.val, node);

  // compute reaction force
  assemblingK(log);
  assemblingKnon(log);
  assemblingKdesign(log);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  // check convergence of residual force
  if(log)
  {
    force.fres = force.fext - force.fint;
    double norm = sqrt(force.fres.dot(force.fres));
    double residual = norm / fnorm;

    std::cout << "first norm "
              << "      res " << std::endl;
    std::cout << fnorm << "       " << residual << std::endl;
  }
}

template <class E, class N> void LinearAnalysisMB<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
      n.val[j] = 0.0;

  /// Value class
  value.val.setZero();
}

template <class E, class N> void LinearAnalysisMB<E, N>::assemblingK(bool log)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  force.fint.setZero();
  globalK.setZero();
  reaction.setZero();

  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> globalR(fem.neq, fem.neq);
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
  globalR.setFromTriplets(tripletsR.begin(), tripletsR.end());

  force.fint = fintMat.diagonal();
  reaction = globalR.diagonal();
  globalK.makeCompressed();

  if(log)
    std::cout << "done" << std::endl;
  fflush(stdout);
}
template <class E, class N>
void LinearAnalysisMB<E, N>::assemblingKnon(bool log)
{
  if(log == 1)
    std::cout << "assembling Knon...";
  fflush(stdout);

  // force.fint.setZero();
  // reaction.setZero();
  globalKnon.setZero();

  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> globalR(fem.neq, fem.neq);
  std::vector<Eigen::Triplet<double>> tripletsKnon, tripletsF, tripletsR;

  // #ifdef _OPENMP
  // #pragma omp parallel
  // #endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripKnon_pri, tripF_pri, tripR_pri;

    // #ifdef _OPENMP
    // #pragma omp for nowait
    // #endif
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      if(element[nel].isDesignable == 0)
        element[nel].makeKeFint(tripKnon_pri, tripF_pri, tripR_pri, fem, node);
      else
        continue;
    }

    // #ifdef _OPENMP
    // #pragma omp critical
    // #endif
    {
      /// connect
      tripletsKnon.insert(tripletsKnon.end(), tripKnon_pri.begin(),
                          tripKnon_pri.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR.insert(tripletsR.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalKnon.setFromTriplets(tripletsKnon.begin(), tripletsKnon.end());
  fintMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  globalR.setFromTriplets(tripletsR.begin(), tripletsR.end());

  // force.fint = fintMat.diagonal();
  // reaction = globalR.diagonal();
  globalKnon.makeCompressed();

  if(log)
    std::cout << "done" << std::endl;
  fflush(stdout);
}
template <class E, class N>
void LinearAnalysisMB<E, N>::assemblingKdesign(bool log)
{
  if(log == 1)
    std::cout << "assembling Kdesign...";
  fflush(stdout);

  // force.fint.setZero();
  // reaction.setZero();
  globalKdesign.setZero();

  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> globalR(fem.neq, fem.neq);
  std::vector<Eigen::Triplet<double>> tripletsKdesign, tripletsF, tripletsR;

  // #ifdef _OPENMP
  // #pragma omp parallel
  // #endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripKdesign_pri, tripF_pri, tripR_pri;

    // #ifdef _OPENMP
    // #pragma omp for nowait
    // #endif
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      if(element[nel].isDesignable == 1)
        element[nel].makeKeFint(tripKdesign_pri, tripF_pri, tripR_pri, fem,
                                node);
      else
        continue;
    }

    // #ifdef _OPENMP
    // #pragma omp critical
    // #endif
    {
      /// connect
      tripletsKdesign.insert(tripletsKdesign.end(), tripKdesign_pri.begin(),
                             tripKdesign_pri.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR.insert(tripletsR.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalKdesign.setFromTriplets(tripletsKdesign.begin(), tripletsKdesign.end());
  fintMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  globalR.setFromTriplets(tripletsR.begin(), tripletsR.end());

  // force.fint = fintMat.diagonal();
  // reaction = globalR.diagonal();
  globalKdesign.makeCompressed();

  if(log)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

} // namespace icarat