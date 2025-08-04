///  @file  opt_multibody.hpp
///  @author  Takeshi Chang
///  @date  Nov 7, 2022.
/// Copyright (c) 2022 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "multibody_presetting.hpp"
#include <analysis/linear.hpp>
#include <analysis/multibody_linear.hpp>
#include <analysis/solver/solvers.hpp>
#include <basis/bmatrix.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <io/mesh_mdpa.hpp>
#include <io/mesh_vtk.hpp>
#include <misc/assembling.hpp>
#include <optimization/filter_multibody.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt-solver/OCSolver.hpp>
#include <optimization/opt_base.hpp>
using namespace std;
using namespace icarat;
using namespace icarat::multibody;
using namespace Eigen;
namespace icarat
{
namespace multibody
{
void fdm(FEM &fem, std::vector<ElementMB> &element, std::vector<NodeMB> &node,
         Force &force, Value &value, Optimize &optim,
         FilterMultibody<ElementMB, NodeMB> &filter, double &T);
void fdm_constraint(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
                    Force &force, Value &dis, Optimize &optim,
                    FilterMultibody<ElementMB, NodeMB> &filter, double &T,
                    Eigen::VectorXd wvector, double &allconst);
void fdm_compliance_constraint(FEM &fem, std::vector<ElementMB> &element,
                               std::vector<NodeMB> &node, Force &force,
                               Value &value, Optimize &optim,
                               FilterMultibody<ElementMB, NodeMB> &filter,
                               double &T, double &C);
void sensitivity(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
                 Optimize &optim, Value &dis,
                 LinearAnalysisMB<ElementMB, NodeMB> &compt,
                 Eigen::VectorXd &Evector);
void calculate_dKds(FEM &fem, vector<ElementMB> &element,
                    Eigen::VectorXd &Evector);
void objective_sensitivity(FEM &fem, vector<ElementMB> &element,
                           vector<NodeMB> &node, Optimize &optim, Value &dis,
                           LinearAnalysisMB<ElementMB, NodeMB> &compt,
                           Eigen::VectorXd &Evector);
void volume_constraint_sensitivity(FEM &fem, vector<ElementMB> &element,
                                   Optimize &optim, double &V0,
                                   vector<double> &dg1ds);
void displacement_constraint_sensitivity(
    FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node, Optimize &optim,
    double &allconst, Value &dis, LinearAnalysisMB<ElementMB, NodeMB> &compt,
    Eigen::VectorXd &wvector, vector<double> &dg2ds, Eigen::VectorXd &Evector);
void compliance_constraint_sensitivity(
    FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node, Optimize &optim,
    Value &dis, LinearAnalysisMB<ElementMB, NodeMB> &compt,
    vector<double> &dg2ds, Eigen::VectorXd &Evector, double &C);
void volume_constraint(FEM &fem, Optimize &optim, vector<ElementMB> &element,
                       double &V0);
void displacement_constraint(FEM &fem, Optimize &optim, Value &dis,
                             Eigen::VectorXd wvector, double &allconst);
void compliance_constraint(FEM &fem, Optimize &optim, Value &dis,
                           vector<ElementMB> &element,
                           LinearAnalysisMB<ElementMB, NodeMB> compt,
                           const toml::value &config, double &C);
void selector(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
              const toml::value &config, double &numconstnode,
              double &constvalue, int &numDesign);
void makelocal(FEM &fem, vector<NodeMB> &node, Eigen::MatrixXd &W,
               Eigen::VectorXd &wvector);
} // namespace multibody
} // namespace icarat