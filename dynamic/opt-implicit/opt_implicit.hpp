#pragma once
#include <cmath>
#include <experimental/filesystem>
#include <iostream>
#include <vector>

#include <analysis/implicit_dyna.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <optimization/filter.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat
{
namespace opt_implicit
{
struct FEMDyna : FEM
{
  int timestep;
  double deltaT;
  double gamma;
  double beta;
  double aa;
  double bb;
};

class Sensitivity
{
public:
  Sensitivity(FEMDyna &fem, std::vector<DynamicElastic> &element,
              std::vector<Node> &node, Force &force, ValueDyna &dis,
              Optimize &optim, int timestep)
      : fem_(fem), element_(element), node_(node), force_(force), dis_(dis),
        optim_(optim), dMdsu_(timestep), dCdsu_(timestep), dKdsu_(timestep)
  {
    for(int i = 0; i < timestep; i++)
    {
      dMdsu_[i].resize(fem_.nelx);
      dCdsu_[i].resize(fem_.nelx);
      dKdsu_[i].resize(fem_.nelx);
    }
  }

  void getDerivatives(int time);

  void adjointAndAssembling(double V0);

private:
  void assembleDerivatives(int time);

  FEMDyna &fem_;
  std::vector<DynamicElastic> &element_;
  std::vector<Node> &node_;
  Force &force_;
  ValueDyna &dis_;
  Optimize &optim_;
  std::vector<std::vector<Eigen::VectorXd>> dMdsu_; //[timestep][nelx][numdof]
  std::vector<std::vector<Eigen::VectorXd>> dCdsu_;
  std::vector<std::vector<Eigen::VectorXd>> dKdsu_;
};

void fdm(FEMDyna &fem, std::vector<DynamicElastic> &element,
         std::vector<Node> &node, Force &force, ValueDyna &dis, Optimize &optim,
         Filter<DynamicElastic, Node> &filter);

} // namespace implicit
} // namespace icarat