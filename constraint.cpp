///  @file    constraint.cpp
///  @author	Takeshi Chang
///  @date		June 12, 2023.
#include "opt_multibody.hpp"
using namespace std;
using namespace Eigen;
namespace icarat
{
namespace multibody
{
void volume_constraint(FEM &fem, Optimize &optim, vector<ElementMB> &element,
                       double &V0)
{
  if(optim.optstep == 1)
  {
    V0 = 0.0;
    for(int i = 0; i < fem.nelx; i++)
    {
      if(element[i].isDesignable == 1)
      {
        V0 += element[i].design_s() * element[i].volume;
      }
    }
    if(V0 == 0.0)
      cerr << "V0 = 0. in optimization" << endl;
    // optim.const_h[0] = 0.0;
  }
  else
  {
    double Vratio = 0.0;
    for(int i = 0; i < fem.nelx; i++)
    {
      if(element[i].isDesignable == 1)
      {
        Vratio += element[i].design_s() * element[i].volume / V0;
      }
    }
    optim.const_h[0] = Vratio - 1.0;
    // optim.const_h[0] = 1 - Vratio;
    cout << "const value1: " << optim.const_h[0] << endl;
  }
}
void displacement_constraint(FEM &fem, Optimize &optim, Value &dis,
                             Eigen::VectorXd wvector, double &allconst)
{
  Eigen::VectorXd dis_sqr = Eigen::VectorXd::Zero(fem.numeq);
  for(int i = 0; i < fem.numeq; i++)
  {
    dis_sqr[i] = dis.val[i] * dis.val[i];
  }
  double actualdis = wvector.transpose() * dis_sqr;
  optim.const_h[1] = actualdis - allconst;
  // optim.const_h[1] = (actualdis / allconst) - 1;
  cout << "const value2: " << optim.const_h[1] << endl;
}
void compliance_constraint(FEM &fem, Optimize &optim, Value &dis,
                           vector<ElementMB> &element,
                           LinearAnalysisMB<ElementMB, NodeMB> compt,
                           const toml::value &config, double &C)
{
  double C0;

  if(optim.optstep == 1)
  {
    double alpha =
        toml::find<double>(config, "compliance_constraint_ratio", "ratio");
    C0 = dis.val.transpose() * compt.Kdesign() * dis.val;
    C = C0 * alpha;
    optim.const_h[1] = (C0 / C) - 1;
  }
  else
  {
    double design_domain_compliance =
        dis.val.transpose() * compt.Kdesign() * dis.val;
    optim.const_h[1] = (design_domain_compliance / C) - 1;
    cout << "const value2: " << optim.const_h[1] << endl;
  }
}
} // namespace multibody
} // namespace icarat