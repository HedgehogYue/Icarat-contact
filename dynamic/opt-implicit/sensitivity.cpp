#include "opt_implicit.hpp"

namespace icarat
{
namespace opt_implicit
{
using namespace std;
using namespace Eigen;

void Sensitivity::getDerivatives(int time)
{
  assert(element_[0].young2() > element_[0].young1());

  for(int nel = 0; nel < fem_.nelx; nel++)
  {
    double design = element_[nel].design_s();

    double young1 = element_[nel].young1();
    double young2 = element_[nel].young2();
    double p = element_[nel].pp();
    double young = (1.0 - pow(design, p)) * young1 + pow(design, p) * young2;
    double dyoung = p * pow(design, p - 1.0) * (young2 - young1);
    double E = dyoung / young;

    double density1 = 0.0; // void density
    double density2 = element_[nel].density();
    double density = design * density2;
    double ddensity = density2;
    double D = ddensity / density;

    // get element value
    Eigen::VectorXd disp = Eigen::VectorXd::Zero(element_[nel].numdof);
    Eigen::VectorXd velo = Eigen::VectorXd::Zero(element_[nel].numdof);
    Eigen::VectorXd acce = Eigen::VectorXd::Zero(element_[nel].numdof);
    for(int i = 0; i < element_[nel].ne; i++)
    {
      for(int j = 0; j < fem_.ndim; j++)
      {
        int dof = node_[element_[nel].nodeID[i]].dof[j];
        // if dirichlet is applied
        if(dof >= fem_.numeq)
        {
          double val = node_[element_[nel].nodeID[i]].dirich->val[j];
          disp.coeffRef(fem_.ndim * i + j) = val;
          velo.coeffRef(fem_.ndim * i + j) = 0.0;
          acce.coeffRef(fem_.ndim * i + j) = 0.0;
        }
        else
        {
          disp.coeffRef(fem_.ndim * i + j) = dis_.val.coeff(dof);
          velo.coeffRef(fem_.ndim * i + j) = dis_.velocity.coeff(dof);
          acce.coeffRef(fem_.ndim * i + j) = dis_.accel.coeff(dof);
        }
      }
    }

    dMdsu_[time - 1][nel] = D * element_[nel].Me() * acce;

    dKdsu_[time - 1][nel] = E * element_[nel].Ke() * disp;

    dCdsu_[time - 1][nel] =
        (fem_.aa * D * element_[nel].Me() + fem_.bb * E * element_[nel].Ke()) *
        velo;
  }
}

void Sensitivity::adjointAndAssembling(double V0)
{
  cout << "---------- computing adjoint vectors ----------" << endl;

  // init
  for(int i = 0; i < fem_.nelx; i++)
  {
    optim_.dfds[i] = 0.0;
    optim_.dgds[i] = 0.0;
  }
  // compute adjoint vector
  ImplicitDyna<DynamicElastic, Node> compt(fem_, element_, node_, force_, dis_);
  compt.initialization();
  for(int time = 1; time < fem_.timestep + 1; time++)
  {
    int invTime = fem_.timestep + 1 - time;
    force_.fext = -invTime * 2 * force_.fext_org / fem_.timestep;

    cout << "---------- timestep no. " << time << " ----------" << endl;

    compt.calRayleigh(time, fem_.timestep, fem_.aa, fem_.bb, fem_.deltaT,
                      fem_.beta, fem_.gamma, "Eigen_CG", 0);

    assembleDerivatives(time);
  }
}

void Sensitivity::assembleDerivatives(int time)
{
  for(int nel = 0; nel < fem_.nelx; nel++)
  {
    // get element value
    Eigen::VectorXd disp = Eigen::VectorXd::Zero(element_[nel].numdof);
    for(int i = 0; i < element_[nel].ne; i++)
    {
      for(int j = 0; j < fem_.ndim; j++)
      {
        int dof = node_[element_[nel].nodeID[i]].dof[j];
        // if dirichlet is applied
        if(dof > fem_.numeq)
        {
          double val = node_[element_[nel].nodeID[i]].dirich->val[j];
          disp.coeffRef(fem_.ndim * i + j) = val;
        }
        else
          disp.coeffRef(fem_.ndim * i + j) = dis_.val.coeff(dof);
      }
    }

    optim_.dfds[nel] += disp.dot(dMdsu_[time - 1][nel] + dCdsu_[time - 1][nel] +
                                 dKdsu_[time - 1][nel]);

    optim_.dgds[nel] = element_[nel].volume;
  }
}

} // namespace implicit
} // namespace opt_icarat