///  @file  boundary_condition_base.cpp
///  @author  Daiki Watanabe
///  @date  February 9, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "boundary_condition_base.hpp"
#include <basis/nmatrix_surface.hpp>

namespace icarat
{
using namespace std;
using namespace Eigen;

int SurfaceLoad::judgeLine(int counter, std::vector<int> hitpoints)
{
  int side;
  if(actele_.eType == "tria3" || actele_.eType == "tria6")
  {
    if(counter == 1)
      side = 1; // Ymin
    else
    {
      if(hitpoints[1] - hitpoints[0] == 1)
        side = 2; // Ymax
      else
        side = 3; // Xmin
    }
  }
  else if(actele_.eType == "quad4" || actele_.eType == "quad8")
  {
    if(counter == 1)
      side = 1; // Ymin
    else if(counter == 2)
      side = 2; // Xmax
    else
    {
      if(hitpoints[1] - hitpoints[0] == 1)
        side = 3; // Ymax
      else
        side = 4; // Xmin
    }
  }
  else if(actele_.eType == "tetra4" || actele_.eType == "tria10")
  {
    if(counter == 2)
      side = 1;
    else
    {
      if(hitpoints[2] + hitpoints[1] + hitpoints[0] == 4)
        side = 2;
      else if(hitpoints[2] + hitpoints[1] + hitpoints[0] == 6)
        side = 3;
      else
        side = 4;
    }
  }
  else if(actele_.eType == "hexa8" || actele_.eType == "hexa20")
  {
    if(counter == 3)
      side = 1; // Zmin
    else if(counter == 5)
      side = 3; // Ymin
    else if(counter == 6)
      side = 4; // Xmax
    else
    {
      if(hitpoints[2] - hitpoints[1] == 1)
        side = 2; // Xmin
      else if(hitpoints[1] - hitpoints[0] == 1)
        side = 5; // Ymax
      else
        side = 6; // Zmax
    }
  }
  else
  {
    cerr << "No load can be applied to this element type in "
            "boundary_condition_base.cpp "
         << endl;
    exit(1);
  }

  return side;
}

int SurfaceLoad::setNumNode()
{

  int numnode;
  if(actele_.eType == "tria3" || actele_.eType == "tria6")
    numnode = 2;
  else if(actele_.eType == "quad4" || actele_.eType == "quad8")
    numnode = 2;
  else if(actele_.eType == "tetra4" || actele_.eType == "tetra10")
    numnode = 3;
  else if(actele_.eType == "hexa8" || actele_.eType == "hexa20")
    numnode = 4;
  else
  {
    cerr << "No load can be applied to this element type in "
            "boundary_condition_base.cpp "
         << endl;
    exit(1);
  }

  return numnode;
}

int SurfaceLoad::surfaceIpmax()
{
  // set number of gauss point for appried surface
  int ipmax;
  // line2
  if(actele_.eType == "tria3" || actele_.eType == "quad4" ||
     actele_.eType == "mitc4")
    ipmax = 2;
  // line3
  else if(actele_.eType == "quad8" || actele_.eType == "tria6")
    ipmax = 3;
  // tria3
  else if(actele_.eType == "tetra4")
    ipmax = 1;
  // tria6
  else if(actele_.eType == "tetra10")
    ipmax = 3;
  // quad4
  else if(actele_.eType == "hexa8")
    ipmax = 4;
  // quad8
  else if(actele_.eType == "hexa20")
    ipmax = 9;
  else
  {
    std::cerr << "Error in boundary_condition_base.hpp" << std::endl;
    exit(1);
  }

  return ipmax;
}

void SurfaceLoad::generate(Force &force, Eigen::VectorXi &idof,
                           Eigen::MatrixXd &X, InputLoad &inp_f)
{
  // get surface parameters
  string eTypeSur;
  int neSur, ipmaxSur;
  tie(eTypeSur, neSur, ipmaxSur) =
      resizeSurfaceParameter(idof, X, actele_.eType, inp_f.sideID);

  int min = std::min(3, fem_.dofnp);
  Eigen::VectorXd val = Eigen::VectorXd::Zero(min);
  for(int i = 0; i < min; i++)
    val[i] = inp_f.val[i];

  Eigen::MatrixXd Ne = Eigen::MatrixXd::Zero(min, idof.size());
  Eigen::VectorXd fe = Eigen::VectorXd::Zero(idof.size());

  SurfaceNmatrix nmatrix;
  for(int ip = 0; ip < ipmaxSur; ip++)
  {
    double jac = 0.0;

    // make shape function for surface
    nmatrix.make(Ne, jac, eTypeSur, ipmaxSur, X, ip);

    fe += Ne.transpose() * val * jac;
  }

  // assembling
  for(int i = 0; i < idof.size(); i++)
  {
    // check boundary condition
    if(idof.coeff(i) >= fem_.numeq)
      continue;

    force.fext_org(idof[i]) += fe[i];
  }
}

std::tuple<std::string, int, int>
SurfaceLoad::resizeSurfaceParameter(Eigen::VectorXi &idof, Eigen::MatrixXd &X,
                                    std::string eType, int sideID)
{
  int ne, ipmax;
  string eTypeSur; ///< element type of surface
  VectorXi neID;

  if(eType == "tria3")
  {
    eTypeSur = "line2";
    ne = 2;
    ipmax = 2;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1;
    else if(sideID == 2)
      neID << 1, 2;
    else if(sideID == 3)
      neID << 2, 0;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "tria6")
  {
    eTypeSur = "line3";
    ne = 3;
    ipmax = 3;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 3;
    else if(sideID == 2)
      neID << 1, 2, 4;
    else if(sideID == 3)
      neID << 2, 0, 5;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "quad4")
  {
    eTypeSur = "line2";
    ne = 2;
    ipmax = 2;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1;
    else if(sideID == 2)
      neID << 1, 2;
    else if(sideID == 3)
      neID << 2, 3;
    else if(sideID == 4)
      neID << 3, 0;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "quad8")
  {
    eTypeSur = "line3";
    ne = 3;
    ipmax = 3;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 4;
    else if(sideID == 2)
      neID << 1, 2, 5;
    else if(sideID == 3)
      neID << 2, 3, 6;
    else if(sideID == 4)
      neID << 3, 0, 7;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "tetra4")
  {
    eTypeSur = "tria3";
    ne = 3;
    ipmax = 1;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2;
    else if(sideID == 2)
      neID << 0, 1, 3;
    else if(sideID == 3)
      neID << 1, 2, 3;
    else if(sideID == 4)
      neID << 2, 0, 3;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "tetra10")
  {
    eTypeSur = "tria6";
    ne = 6;
    ipmax = 3;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2, 4, 5, 6;
    else if(sideID == 2)
      neID << 0, 1, 3, 4, 8, 7;
    else if(sideID == 3)
      neID << 1, 2, 3, 5, 9, 8;
    else if(sideID == 4)
      neID << 2, 0, 3, 6, 7, 9;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "hexa8")
  {
    eTypeSur = "quad4";
    ne = 4;
    ipmax = 4;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2, 3;
    else if(sideID == 2)
      neID << 3, 0, 4, 7;
    else if(sideID == 3)
      neID << 0, 1, 5, 4;
    else if(sideID == 4)
      neID << 1, 2, 6, 5;
    else if(sideID == 5)
      neID << 2, 3, 7, 6;
    else if(sideID == 6)
      neID << 4, 5, 6, 7;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "hexa20")
  {
    eTypeSur = "quad8";
    ne = 8;
    ipmax = 9;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2, 3, 8, 9, 10, 11;
    else if(sideID == 2)
      neID << 3, 0, 4, 7, 11, 16, 15, 19;
    else if(sideID == 3)
      neID << 0, 1, 5, 4, 8, 17, 12, 16;
    else if(sideID == 4)
      neID << 1, 2, 6, 5, 9, 18, 13, 17;
    else if(sideID == 5)
      neID << 2, 3, 7, 6, 10, 19, 14, 18;
    else if(sideID == 6)
      neID << 4, 5, 6, 7, 12, 13, 14, 15;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else
  {
    std::cerr << "Distributed loads cannot be applied to this element type in "
                 "boundary_condition_base.cpp."
              << endl;
    exit(1);
  }

  // convert to surface parameter
  MatrixXd Xtmp = X;
  VectorXi idoftmp = idof;
  X.resize(ne, fem_.ndim);
  idof.resize(fem_.dofnp * ne);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem_.dofnp; j++)
      idof.coeffRef(fem_.dofnp * i + j) = idoftmp(fem_.dofnp * neID(i) + j);
    for(int j = 0; j < fem_.ndim; j++)
      X.coeffRef(i, j) = Xtmp(neID(i), j);
  }

  return forward_as_tuple(eTypeSur, ne, ipmax);
}

} // namespace icarat