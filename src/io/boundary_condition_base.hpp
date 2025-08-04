///  @file  boundary_condition_base.hpp
///  @author  Daiki Watanabe
///  @date  February 9, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <Eigen/Dense>
#include <io/mesh_icarat.hpp>
#include <numeric>
#include <problem/base.hpp>

namespace icarat
{

/// set Dirichlet object & value to node
template <class N>
void setDirichValue(int dofnp, N &node, std::vector<int> &flag,
                    std::vector<double> &val)
{
  if(node.dirich == nullptr)
  {
    std::shared_ptr<Dirichlet> dirich(new Dirichlet);
    node.dirich = dirich;
  }
  for(int k = 0; k < dofnp; k++)
  {
    if(node.dirich->flag[k] != 1)
      node.dirich->flag[k] = flag[k];
    node.dirich->val[k] = val[k];
  }
}

template <class N> int countDOF(FEM &fem, N &node)
{
  int counter = 0;
  for(auto &n : node)
  {
    // dirichlet is nothing
    if(n.dirich == nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        n.dof[j] = counter;
        counter++;
      }
    }
    else
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        // dirichlet is not active
        if(n.dirich->flag[j] == 0)
        {
          n.dof[j] = counter;
          counter++;
        }
        // dirichlet is active... through
        else
          continue;
      }
    }
  }

  // constrainted dof
  int constCounter = counter;
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        // dirichlet is active
        if(n.dirich->flag[j] != 0)
        {
          n.dof[j] = constCounter;
          constCounter++;
        }
      }
    }
  }

  assert(constCounter == fem.neq);

  std::cout << "numeq : " << counter << std::endl;
  return counter;
}

/// base class for load input
struct InputLoad
{
  int sideID;              ///< side ID of the element
  std::vector<double> val; ///< load value
};

/// Class for creating distributed loads.
class SurfaceLoad
{
public:
  SurfaceLoad(FEM &fem, Element &actele) : fem_(fem), actele_(actele) {}

  /// Determine which line (surface) of element is subject to external loads.
  int judgeLine(int counter, std::vector<int> hitpoints);

  int setNumNode();

  int surfaceIpmax();

  ///  compute the equivalent nodal force at an element to Force.fext_org
  void generate(Force &force, Eigen::VectorXi &idof, Eigen::MatrixXd &X,
                InputLoad &inp_f);

protected:
  /// Set and resize the parameters for the element's surface
  std::tuple<std::string, int, int>
  resizeSurfaceParameter(Eigen::VectorXi &idof, Eigen::MatrixXd &X,
                         std::string eType, int sideID);

  FEM &fem_;
  Element &actele_;

  /// tria3,6
  /// 2
  /// |  |
  /// 5     4          y
  /// |        |      ↑
  /// 0----3-----1     →x
  ///
  /// side = 1...0-1
  /// side = 2...1-2
  /// side = 3...2-0

  /// tetra4,10
  /// 3
  /// |  9
  /// 7  8 2            z  y
  /// | 6    5      ↑↗
  /// 0---4----1     →x
  ///
  /// side = 1...0-1-2
  /// side = 2...0-1-3
  /// side = 3...1-2-3
  /// side = 4...2-0-3

  /// quad4,8
  /// 3-----6------2
  /// |            |
  /// 7            5       y
  /// |            |      ↑
  /// 0------4-----1       →x
  ///
  /// side = 1...0-1
  /// side = 2...1-2
  /// side = 3...2-3
  /// side = 4...3-0

  /// hexa8,20
  ///    7-----14------6
  ///  15           13 |
  /// 4-----12-----5  18
  /// | 19         |  |
  /// 16  3    10  17 2     z   y
  /// | 11         | 9      ↑↗
  /// 0------8-----1          →x
  ///
  /// side = 1...0-1-2-3
  /// side = 2...3-0-4-7
  /// side = 3...0-1-5-4
  /// side = 4...1-2-6-5
  /// side = 5...2-3-7-6
  /// side = 6...4-5-6-7
};

} // namespace icarat
