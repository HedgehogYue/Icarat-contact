#pragma once
#include <Eigen/Core>
#include <type_traits>
#include <vector>

namespace icarat
{
template <class N>
void divideValIntoNode(int dofnp, int numeq, Eigen::VectorXd &value,
                       std::vector<N> &node)
{
  int dofnptmp = std::min(3, dofnp);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
  {
    for(int j = 0; j < dofnptmp; j++)
    {
      if(node[i].dof[j] >= numeq)
        node[i].val[j] = node[i].dirich->val[j];
      else
        node[i].val[j] = value.coeff(node[i].dof[j]);

      // Too small value is deleted for paraview
      if(abs(node[i].val[j]) < 1.0e-15)
        node[i].val[j] = 0.0;
    }
  }
}

template <class N>
void divideReactionIntoNode(int dofnp, Eigen::VectorXd &value,
                            std::vector<N> &node)
{
  int dofnptmp = std::min(3, dofnp);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
    for(int j = 0; j < dofnptmp; j++)
    {
      node[i].reaction[j] = value.coeff(node[i].dof[j]);

      // Too small value is deleted for paraview
      if(abs(node[i].reaction[j]) < 1.0e-15)
        node[i].reaction[j] = 0.0;
    }
}
} // namespace icarat