///  @file    fdm.cpp
///  @author	Takeshi Chang
///  @date		March 1, 2023.

#include "opt_multibody.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace multibody
{
void fdm(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
         Force &force, Value &dis, Optimize &optim,
         FilterMultibody<ElementMB, NodeMB> &filter, double &T)
{

  double delta = 1.0e-5;
  double func_new;
  vector<double> object(fem.nelx);
  double func_old = optim.object_f;
  auto design_s_old = optim.design_s;
  int counter = 0;
  LinearAnalysisMB<ElementMB, NodeMB> compt(fem, element, node, force, dis);
  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    if(element[nel].isDesignable == 0)
    {
      object[nel] = 0;
      cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
    }
    else if(element[nel].isDesignable == 1)
    {
      // set delta
      optim.design_s[counter] += delta;

      // set threshold design variable
      auto design_filtered = filter.densityFilter(optim.design_s, false);

      double beta = 1.0;
      auto design_threshold = filter.threshold(design_filtered, beta, T);
      int numcount = 0;
      for(int i = 0; i < fem.nelx; i++)
      {
        if(element[i].isDesignable == 1)
        {
          // element[i].setDesignS(design_threshold[numcount]);
          // before doing FDM, please change this part into design_s
          element[i].setDesignS(optim.design_s[numcount]);
          numcount++;
        }
      }

      // cal f(s + delta)
      compt.solve("Eigen_CG", 0);

      func_new = dis.val.transpose() * compt.Knon() * dis.val;

      //   func_new = dis.val.dot(force.fext);
      //   func_new = dis.val.transpose() * W * compt.K() * W * dis.val;

      // objective function initialize
      // func_new = 0.0;
      // for(int n = 0; n < fem.nelx; n++)
      // {
      //   if(element[n].isDesignable == 0)
      //   {
      //     Eigen::VectorXd disp_element =
      //         Eigen::VectorXd::Zero(element[n].numdof);
      //     for(int i = 0; i < element[n].ne; i++)
      //       for(int j = 0; j < fem.ndim; j++)
      //       {
      //         disp_element.coeffRef(fem.ndim * i + j) =
      //             node[element[n].nodeID[i]].val[j];
      //       }
      //     func_new += disp_element.transpose() * element[n].Ke() *
      //     disp_element;
      //   }
      //   else
      //     continue;
      // }
      // cout << func_new << endl << func_old;

      /// dfds
      object[nel] = (func_new - func_old) / delta;
      cout << func_new << " " << func_old << endl;
      // reset design variables
      optim.design_s = design_s_old;

      cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
      counter++;
    }
  }

  // output sensivity graph
  exportGraph(object, "FDM_sensitivity_objective.csv");

  cout << "------------------------------------------------" << endl;
  cout << "         Finish computation with FDM            " << endl;
  cout << "------------------------------------------------" << endl;
  // exit(0);
}

void fdm_constraint(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
                    Force &force, Value &dis, Optimize &optim,
                    FilterMultibody<ElementMB, NodeMB> &filter, double &T,
                    Eigen::VectorXd wvector, double &allconst)
{

  double delta = 1.0e-5;
  double func_new;
  vector<double> object(fem.nelx);
  double func_old = optim.const_h[1];
  auto design_s_old = optim.design_s;
  int counter = 0;
  LinearAnalysisMB<ElementMB, NodeMB> compt(fem, element, node, force, dis);
  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int nel = 0; nel < fem.nelx; nel++)
  {

    if(element[nel].isDesignable == 0)
    {
      object[nel] = 0;
      cout << "Element ID : " << nel << " dgds macro = " << object[nel] << endl;
    }
    else if(element[nel].isDesignable == 1)
    {
      // set delta
      optim.design_s[counter] += delta;

      // set threshold design variable
      auto design_filtered = filter.densityFilter(optim.design_s, false);

      double beta = 1.0;
      auto design_threshold = filter.threshold(design_filtered, beta, T);
      int numcount = 0;
      for(int i = 0; i < fem.nelx; i++)
      {
        if(element[i].isDesignable == 1)
        {
          // element[i].setDesignS(design_threshold[numcount]);
          // before doing FDM, please change this part into design_s
          element[i].setDesignS(optim.design_s[numcount]);
          numcount++;
        }
      }

      // cal f(s + delta)
      compt.solve("Eigen_CG", 0);
      Eigen::VectorXd dis_sqr = Eigen::VectorXd::Zero(fem.numeq);
      for(int i = 0; i < fem.numeq; i++)
      {
        dis_sqr[i] = dis.val[i] * dis.val[i];
      }
      double actualdis = wvector.transpose() * dis_sqr;
      // double actualdis = wvector.transpose() * dis.val;
      func_new = actualdis - allconst;
      // func_new = (actualdis / allconst) - 1;

      //   func_new = dis.val.dot(force.fext);
      //   func_new = dis.val.transpose() * W * compt.K() * W * dis.val;

      /// dfds
      object[nel] = (func_new - func_old) / delta;
      cout << func_new << " " << func_old << endl;
      // reset design variables
      optim.design_s = design_s_old;

      cout << "Element ID : " << nel << " dgds macro = " << object[nel] << endl;
      counter++;
    }
  }

  // output sensivity graph
  exportGraph(object, "FDM_constraint_sensitivity_objective.csv");

  cout << "------------------------------------------------" << endl;
  cout << "         Finish computation with FDM            " << endl;
  cout << "------------------------------------------------" << endl;
  exit(0);
}
void fdm_compliance_constraint(FEM &fem, std::vector<ElementMB> &element,
                               std::vector<NodeMB> &node, Force &force,
                               Value &dis, Optimize &optim,
                               FilterMultibody<ElementMB, NodeMB> &filter,
                               double &T, double &C)
{
  double delta = 1.0e-5;
  double func_new;
  vector<double> object(fem.nelx);
  double func_old = optim.const_h[1];
  auto design_s_old = optim.design_s;
  LinearAnalysisMB<ElementMB, NodeMB> compt(fem, element, node, force, dis);
  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    if(element[nel].isDesignable == 0)
    {
      object[nel] = 0;
      cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
    }
    else if(element[nel].isDesignable == 1)
    {
      // set delta
      optim.design_s[nel] += delta;

      // set threshold design variable
      auto design_filtered = filter.densityFilter(optim.design_s, false);

      double beta = 1.0;
      auto design_threshold = filter.threshold(design_filtered, beta, T);
      int numcount = 0;
      for(int i = 0; i < fem.nelx; i++)
      {
        if(element[i].isDesignable == 1)
        {
          // element[i].setDesignS(design_threshold[numcount]);
          // before doing FDM, please change this part into design_s
          element[i].setDesignS(optim.design_s[numcount]);
          numcount++;
        }
      }

      // cal f(s + delta)
      compt.solve("Eigen_CG", 0);
      double design_domain_compliance =
          dis.val.transpose() * compt.Kdesign() * dis.val;
      func_new = (design_domain_compliance / C) - 1;

      //   func_new = dis.val.dot(force.fext);
      //   func_new = dis.val.transpose() * W * compt.K() * W * dis.val;

      /// dfds
      object[nel] = (func_new - func_old) / delta;
      cout << func_new << " " << func_old << endl;
      // reset design variables
      optim.design_s = design_s_old;

      cout << "Element ID : " << nel << " dgds macro = " << object[nel] << endl;
    }
  }
  // output sensivity graph
  exportGraph(object, "FDM_compliance_constraint_sensitivity.csv");

  cout << "------------------------------------------------" << endl;
  cout << "         Finish computation with FDM            " << endl;
  cout << "------------------------------------------------" << endl;
  if(optim.optstep == 1)
    exit(0);
}
} // namespace multibody
} // namespace icarat