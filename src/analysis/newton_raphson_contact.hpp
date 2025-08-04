///  @file  newton_raphson_contact.hpp
///  @author  Takeshi Chang
///  @date  March 13, 2025.
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
template <class E, class N> class NewtonRaphsonContact
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
  Eigen::VectorXd reaction;            ///< reaction forcepublic:

public:
  NewtonRaphsonContact(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                       Force &force, Value &value, int numpair);

  ///  main function (Incremental load or displacement method)
  ///@return Newton-Raphson iteration
  int solve(int thisTime, int totalTime, double tolNR, int maxIter,
            std::string solvName, int log, int numpair,
            std::vector<std::vector<int>> node_pair, double epsilon);

  /// main function (arc length method)
  void solveArcLength(int thisStep, int totalStep, std::string solvName,
                      int log);

  /// initializaiton
  void initialization();
  /// calculate gap
  void get_gap(Eigen::VectorXd &gap, int numpair,
               std::vector<std::vector<int>> node_pair);
  /// assemble K
  void assemblingK(bool log);
  void assemblingKnon(bool log);
  void assemblingKdesign(bool log);
  /// assemble Cn_bar
  void assembling_Cn_bar(Eigen::MatrixXd &Cn_bar, Eigen::VectorXd gap,
                         int numpair, std::vector<std::vector<int>> node_pair);

  void divideValIntoNode(int dofnp, Eigen::VectorXd &val, Eigen::VectorXd &du,
                         int thisStep, int totalStep);

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }
  const Eigen::SparseMatrix<double> &Knon() const { return globalKnon; }
  const Eigen::SparseMatrix<double> &Kdesign() const { return globalKdesign; }
};

template <class E, class N>
NewtonRaphsonContact<E, N>::NewtonRaphsonContact(FEM &fem,
                                                 std::vector<E> &element,
                                                 std::vector<N> &node,
                                                 Force &force, Value &value,
                                                 int numpair)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq),
      globalKnon(fem.numeq, fem.numeq), globalKdesign(fem.numeq, fem.numeq)
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;
}

template <class E, class N>
void NewtonRaphsonContact<E, N>::get_gap(
    Eigen::VectorXd &gap, int numpair, std::vector<std::vector<int>> node_pair)
{
  double pos1;
  double pos2;
  double temp_gap;
  for(int i = 0; i < numpair; i++)
  {
    pos1 = node[node_pair[0][i]].x[1] + node[node_pair[0][i]].val[1];
    pos2 = node[node_pair[1][i]].x[1] + node[node_pair[1][i]].val[1];
    temp_gap = pos2 - pos1;
    gap(i) = temp_gap;
    /// can only be used in y direction problem
  }
}

template <class E, class N>
void NewtonRaphsonContact<E, N>::assembling_Cn_bar(
    Eigen::MatrixXd &Cn_bar, Eigen::VectorXd gap, int numpair,
    std::vector<std::vector<int>> node_pair)
{
  for(int i = 0; i < numpair; i++) /// i is pair number counter
  {
    int counter = 0; /// counter is the position in the matrix counter
    int j = 0;       /// j is node counter
    /// searching for the first node in the pair
    do
    {
      counter = counter + abs(node[j].nwd);
      j++;
    } while(!(j >= node_pair[0][i]));
    /// found the first node in the pair and give it Cn
    if(j == node_pair[0][i])
    {
      if(counter < 0 || counter >= Cn_bar.cols())
      {
        std::cout << "Error: counter = " << counter
                  << " exceeds Cn_bar.cols() = " << Cn_bar.cols() << std::endl;
        exit(1);
      }
      if(gap(i) >= 0 || node[node_pair[0][i]].node_density <
                            0.01) /// one of gap or density, out
      {
        counter = counter + abs(node[j].nwd);
        j++;
      }
      else
      {
        node[j].contact_active = 1;
        if(node[j].nwd == 2)
        {
          Cn_bar(i, counter) = node[j].Cn[0];
          counter++;
          Cn_bar(i, counter) = node[j].Cn[1];
          counter++;
          j++;
        }
        else if(abs(node[j].nwd) == 1)
        {
          if(node[j].nwd == -1)
            Cn_bar(i, counter) = node[j].Cn[0];

          else if(node[j].nwd == 1)
            Cn_bar(i, counter) = node[j].Cn[1];

          counter++;
          j++;
        }
        else if(node[j].nwd == 0)
          j++;
      }
    }
    else
    {
      std::cout << std::endl << "counting pairs wrong!!" << std::endl;
      exit(0);
    }
    /// searching for the second node in the pair, and skip all the nodes
    /// betweeen
    do
    {
      counter = counter + abs(node[j].nwd);
      j++;
    } while(!(j >= node_pair[1][i]));
    /// found the second node in the pair and give it Cn
    if(j == node_pair[1][i])
    {
      if(counter < 0 || counter >= Cn_bar.cols())
      {
        std::cout << "Error: counter = " << counter
                  << " exceeds Cn_bar.cols() = " << Cn_bar.cols() << std::endl;
        exit(1);
      }
      if(gap(i) >= 0 || node[node_pair[0][i]].node_density < 0.01)
      {
        counter = counter + abs(node[j].nwd);
        j++;
      }
      else
      {
        node[j].contact_active = 1;
        if(node[j].nwd == 2)
        {
          Cn_bar(i, counter) = node[j].Cn[0];
          counter++;
          Cn_bar(i, counter) = node[j].Cn[1];
          counter++;
          j++;
        }
        else if(abs(node[j].nwd) == 1)
        {
          if(node[j].nwd == -1)
          {
            Cn_bar(i, counter) = node[j].Cn[0];
            counter++;
          }
          else if(node[j].nwd == 1)
          {
            Cn_bar(i, counter) = node[j].Cn[1];
            counter++;
          }
          j++;
        }
        else if(node[j].nwd == 0)
          j++;
      }
    }
    else
    {
      std::cout << std::endl << "counting pairs wrong!!" << std::endl;
      exit(0);
    }
  }
}

template <class E, class N>
int NewtonRaphsonContact<E, N>::solve(int thisTime, int totalTime, double tolNR,
                                      int maxIter, std::string solvName,
                                      int log, int numpair,
                                      std::vector<std::vector<int>> node_pair,
                                      double epsilon)
{

  int iterstep = 0;
  double firstNorm, residual = 1;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::MatrixXd Cn_bar = Eigen::MatrixXd::Zero(numpair, fem.numeq);
  Eigen::VectorXd gap = Eigen::VectorXd::Zero(numpair);
  Eigen::MatrixXd Kn = Eigen::MatrixXd::Zero(fem.numeq, fem.numeq);
  Eigen::MatrixXd Ktan = Eigen::MatrixXd::Zero(fem.numeq, fem.numeq);
  Eigen::VectorXd Fn = Eigen::VectorXd::Zero(fem.numeq);
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

    get_gap(gap, numpair, node_pair);
    // assembling
    assemblingK(log);
    assemblingKnon(log);
    assemblingKdesign(log);

    assembling_Cn_bar(Cn_bar, gap, numpair, node_pair);

    /// calculate Kn and Ktan
    Kn = epsilon * Cn_bar.transpose() * Cn_bar;
    Ktan = Kn + globalK;
    Eigen::SparseMatrix<double> Ktan_sparse = Ktan.sparseView();
    /// calculate Fn
    Fn = Kn * value.val;
    // cal residual force
    force.fres = force.fext - force.fint - Fn;

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
    Solvers solv(fem, Ktan_sparse, force.fres, sol);
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

template <class E, class N> void NewtonRaphsonContact<E, N>::initialization()
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

template <class E, class N>
void NewtonRaphsonContact<E, N>::assemblingK(bool log)
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

  if(log)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

template <class E, class N>
void NewtonRaphsonContact<E, N>::assemblingKnon(bool log)
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
void NewtonRaphsonContact<E, N>::assemblingKdesign(bool log)
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

template <class E, class N>
void NewtonRaphsonContact<E, N>::divideValIntoNode(int dofnp,
                                                   Eigen::VectorXd &val,
                                                   Eigen::VectorXd &du,
                                                   int thisStep, int totalStep)
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
