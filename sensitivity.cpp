///  @file    sensitivity.cpp
///  @author	Takeshi Chang
///  @date		May 9, 2023.
#include "opt_multibody.hpp"
using namespace std;
using namespace Eigen;
namespace icarat
{
namespace multibody
{
void sensitivity(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
                 Optimize &optim, Value &dis,
                 LinearAnalysisMB<ElementMB, NodeMB> &compt,
                 Eigen::VectorXd &Evector)
{
  std::cout << "---- start objective sensitivity ----" << endl;
  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd disappear = 2.0 * compt.Kdesign() * dis.val;
  VectorXd zeta = VectorXd::Zero(fem.numeq);
  Solvers zetasolver(fem, compt.K(), disappear, zeta);
  zetasolver.solve("Eigen_CG", 1);
  vector<VectorXd> zetaEle =
      numeqToElement<ElementMB, NodeMB>(fem, element, node, zeta);
  std::cout << "------- finish to compute -------" << endl;
  int counter = 0;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    if(element[nel].isDesignable == 1)
    {
      Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
        for(int j = 0; j < fem.ndim; j++)
        {
          disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
        }
      double part1 = Evector[nel] * disp.transpose() * element[nel].Ke() * disp;
      double part2 =
          -Evector[nel] * zetaEle[nel].transpose() * element[nel].Ke() * disp;
      optim.dfds[counter] = part1 + part2;
      counter++;
    }
  }
  // output sensivity graph
  if(optim.optstep == 1)
    exportGraph(optim.dfds, "sensitivity objective.csv");
  std::cout << "sensitivity calculation finished" << endl;
}
void calculate_dKds(FEM &fem, vector<ElementMB> &element,
                    Eigen::VectorXd &Evector)
{
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();
    double young, dyoung;
    if(young2 >= young1)
    {
      young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
      dyoung = p * pow(dens, p - 1.0) * (young2 - young1);
    }
    else if(young2 < young1)
    {
      young = pow(1.0 - dens, p) * young1 + (1.0 - pow(1.0 - dens, p)) * young2;
      dyoung = (young2 - young1) * p * pow(1.0 - dens, p - 1.0);
    }
    double E = dyoung / young;
    Evector[nel] = E;
  }
}
void objective_sensitivity(FEM &fem, vector<ElementMB> &element,
                           vector<NodeMB> &node, Optimize &optim, Value &dis,
                           LinearAnalysisMB<ElementMB, NodeMB> &compt,
                           Eigen::VectorXd &Evector)
{
  std::cout << "---- start objective sensitivity ----" << endl;
  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd disappearpart = 2.0 * compt.Knon() * dis.val;
  VectorXd lambda = VectorXd::Zero(fem.numeq);
  Solvers lambdasolver(fem, compt.K(), disappearpart, lambda);
  lambdasolver.solve("Eigen_CG", 1);
  vector<VectorXd> lambdaEle =
      numeqToElement<ElementMB, NodeMB>(fem, element, node, lambda);
  std::cout << "------- finish to compute -------" << endl;
  int counter = 0;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    if(element[nel].isDesignable == 1)
    {
      Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
        for(int j = 0; j < fem.ndim; j++)
        {
          disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
        }
      optim.dfds[counter] =
          -Evector[nel] * lambdaEle[nel].transpose() * element[nel].Ke() * disp;
      counter++;
    }
  }
  if(optim.optstep == 1)
    exportGraph(optim.dfds, "sensitivity objective.csv");
  std::cout << "------- objective sensitivity finish -------" << endl;
}
void volume_constraint_sensitivity(FEM &fem, vector<ElementMB> &element,
                                   Optimize &optim, double &V0,
                                   vector<double> &dg1ds)
{
  std::cout << "---- start volume constraint sensitivity ----" << endl;
  int counter = 0;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    if(element[nel].isDesignable == 1)
    {
      dg1ds[counter] = element[nel].volume / V0;
      counter++;
    }
  }
  if(optim.optstep == 1)
    exportGraph(dg1ds, "sensitivity volume constraint.csv");
  std::cout << "------- volume constraint sensitivity finish -------" << endl;
}
void displacement_constraint_sensitivity(
    FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node, Optimize &optim,
    double &allconst, Value &dis, LinearAnalysisMB<ElementMB, NodeMB> &compt,
    Eigen::VectorXd &wvector, vector<double> &dg2ds, Eigen::VectorXd &Evector)
{
  std::cout << "---- start displacement constraint sensitivity ----" << endl;
  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd all = wvector / allconst;
  VectorXd kappa = VectorXd::Zero(fem.numeq);
  VectorXd dis_sens = VectorXd::Zero(fem.numeq);
  for(int i = 0; i < fem.numeq; i++)
  {
    dis_sens[i] = 2 * wvector[i] * dis.val[i];
  }
  Solvers kappasolver(fem, compt.K(), dis_sens, kappa);
  kappasolver.solve("Eigen_CG", 1);
  vector<VectorXd> kappaEle =
      numeqToElement<ElementMB, NodeMB>(fem, element, node, kappa);
  std::cout << "------- finish to compute -------" << endl;
  int counter = 0;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    if(element[nel].isDesignable == 1)
    {
      Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
        for(int j = 0; j < fem.ndim; j++)
        {
          disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
        }
      dg2ds[counter] =
          -Evector[nel] * kappaEle[nel].transpose() * element[nel].Ke() * disp;
      counter++;
    }
  }
  if(optim.optstep == 1)
    exportGraph(dg2ds, "sensitivity displacement constraint.csv");
  std::cout << "------- displacement constraint sensitivity finish -------"
            << endl;
}
void compliance_constraint_sensitivity(
    FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node, Optimize &optim,
    Value &dis, LinearAnalysisMB<ElementMB, NodeMB> &compt,
    vector<double> &dg2ds, Eigen::VectorXd &Evector, double &C)
{
  std::cout << "---- start compliance constraint sensitivity ----" << endl;
  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd adjointpart = (-2.0 / C) * compt.Kdesign() * dis.val;
  VectorXd miu = VectorXd::Zero(fem.numeq);
  Solvers miusolver(fem, compt.K(), adjointpart, miu);
  miusolver.solve("Eigen_CG", 1);
  vector<VectorXd> miuEle =
      numeqToElement<ElementMB, NodeMB>(fem, element, node, miu);
  std::cout << "------- finish to compute -------" << endl;
  int counter = 0;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    if(element[nel].isDesignable == 1)
    {
      Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
        for(int j = 0; j < fem.ndim; j++)
        {
          disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
        }
      double part1 = (1.0 / C) * Evector[nel] * disp.transpose() *
                     element[nel].Ke() * disp;
      double part2 =
          Evector[nel] * miuEle[nel].transpose() * element[nel].Ke() * disp;
      dg2ds[counter] = part1 + part2;
      counter++;
    }
  }
  if(optim.optstep == 1)
    exportGraph(dg2ds, "sensitivity compliance constraint.csv");
  std::cout << "------- compliance constraint sensitivity finish -------"
            << endl;
}
} // namespace multibody
} // namespace icarat