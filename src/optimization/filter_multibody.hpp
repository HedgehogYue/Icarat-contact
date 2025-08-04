///  @file  filter_multibody.hpp
///  @author  Takeshi Chang
///  @date  June 26, 2023.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "opt_base.hpp"
#include <problem/base.hpp>
#include <vector>

namespace icarat
{
/// Density or Sensitivity filter class in TO
/// This class can be only applied for uniform cube or square mesh
template <class E, class N> class FilterMultibody
{
public:
  /// @param [in] radius_ :filter radius
  FilterMultibody(FEM &fem_, std::vector<E> &element_, std::vector<N> &node_,
                  double radius_, int numDesign);

  ///  apply sensitivity Filter
  std::vector<double> sensitivityFilter(const std::vector<double> &design_s,
                                        const std::vector<double> &dfds);

  /// apply density filter
  /// @returns filtered design_s and sensitivity derivative term of design_s
  /// @param  [in] isSens true...sensitivity derivative term,
  /// false...density filter procedure
  std::vector<double> densityFilter(const std::vector<double> &design_s,
                                    bool isSens);

  /// update design variable using threshold function
  std::vector<double> threshold(const std::vector<double> &design_s,
                                double beta, double T);

  /// differential derivative of threshold function
  std::vector<double> thresholdDerivative(const std::vector<double> &design_s,
                                          double beta, double T);

private:
  /// make numele,ID_,dis_,dfilters_ds
  void preNeighbour();

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  double radius_;           ///< filter radius
  std::vector<int> numele_; //< number of the adjacent element for an element i
  std::vector<std::vector<int>>
      ID_; ///< ID[i][j] of the adjacent element j for an element i
  std::vector<std::vector<double>>
      dis_; ///< distance[i][j] from an element i to the adjacent element j
  std::vector<double> dfilters_ds; ///< it is used in "densityFilterSensitivity"
};

template <class E, class N>
FilterMultibody<E, N>::FilterMultibody(FEM &fem_, std::vector<E> &element_,
                                       std::vector<N> &node_, double radius,
                                       int numDesign)
    : fem(fem_), element(element_), node(node_), radius_(radius),
      numele_(numDesign), ID_(numDesign), dis_(numDesign),
      dfilters_ds(numDesign)
{
  preNeighbour();
}

///////////////////////////////////////////////////////////
////////////////////////////public////////////////////////
//////////////////////////////////////////////////////////

template <class E, class N>
std::vector<double>
FilterMultibody<E, N>::sensitivityFilter(const std::vector<double> &design_s,
                                         const std::vector<double> &dfds)
{
  std::vector<double> dfds_new(dfds.size()); ///< non-filtered sensitivity

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < fem.nelx; i++)
  {
    double aa = 0.0;
    double bb = 0.0;

    for(int j = 0; j < numele_[i]; j++)
    {
      double weight = radius_ - dis_[i][j];
      aa += weight;
      bb += (weight * design_s[ID_[i][j]] * dfds[ID_[i][j]]) /
            (design_s[ID_[i][j]] * element[ID_[i][j]].volume);
    }

    double cc = (design_s[i] * element[i].volume) / (design_s[i] * aa);
    dfds_new[i] = cc * bb;
  }

  return dfds_new;
}

template <class E, class N>
std::vector<double>
FilterMultibody<E, N>::densityFilter(const std::vector<double> &design_s,
                                     bool isSens)
{
  std::vector<double> filtered_s(design_s.size());

  // density filter
  if(isSens == false)
  {
    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    int counter = 0;
    for(int i = 0; i < fem.nelx; i++)
    {
      if(element[i].isDesignable == 1)
      {
        double aa = 0.0;
        double bb = 0.0;

        for(int j = 0; j < numele_[counter]; j++)
        {
          double weight = radius_ - dis_[counter][j];
          aa += weight;
          bb += weight * design_s[ID_[counter][j]];
        }

        if(aa < 1.0e-10)
        {
          std::cerr << "coefficient aa is too small in filter.hpp."
                    << std::endl;
          exit(1);
        }
        filtered_s[counter] = bb / aa;
        dfilters_ds[counter] = aa;
        counter++;
      }
    }
  }
  // density filter's sensitivity
  else
  {
    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    for(int i = 0; i < (int)design_s.size(); i++)
    {
      for(int j = 0; j < numele_[i]; j++)
      {
        double weight = radius_ - dis_[i][j];
        filtered_s[i] += weight * design_s[ID_[i][j]] / dfilters_ds[ID_[i][j]];
      }
    }
  }
  return filtered_s;
}

template <class E, class N>
std::vector<double>
FilterMultibody<E, N>::threshold(const std::vector<double> &d_filtered,
                                 double beta, double T)
{
  std::vector<double> hat((int)d_filtered.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)hat.size(); i++)
    hat[i] = 0.5 * (tanh(T * beta) + tanh(beta * (d_filtered[i] - T))) /
             tanh(T * beta);

  return hat;
}

template <class E, class N>
std::vector<double>
FilterMultibody<E, N>::thresholdDerivative(const std::vector<double> &d_fil,
                                           double beta, double T)
{
  std::vector<double> dshat((int)d_fil.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)dshat.size(); i++)
    dshat[i] = (0.5 * (beta * (1.0 - pow(tanh(beta * (d_fil[i] - T)), 2.0))) /
                tanh(T * beta));

  return dshat;
}

///////////////////////////////////////////////////////////
////////////////////////////private////////////////////////
///////////////////////////////////////////////////////////

template <class E, class N> void FilterMultibody<E, N>::preNeighbour()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  int counteri = 0;
  for(int i = 0; i < fem.nelx; i++)
  {
    if(element[i].isDesignable == 0)
      continue;
    if(element[i].isDesignable == 1)
    {
      int nei = element[i].ne;
      // coordinate of element's center
      Eigen::VectorXd centerI(fem.ndim);
      centerI.setZero();
      for(int dim = 0; dim < fem.ndim; dim++)
        for(int e = 0; e < nei; e++)
          centerI[dim] += node[element[i].nodeID[e]].x[dim] / nei;

      std::vector<int> IDi;
      std::vector<double> disi;
      int numj = 0;
      int counterj = 0;
      for(int j = 0; j < fem.nelx; j++)
      {
        if(element[j].isDesignable == 1)
        {
          int nej = element[j].ne;
          Eigen::VectorXd centerJ(fem.ndim);
          centerJ.setZero();

          for(int dim = 0; dim < fem.ndim; dim++)
            for(int e = 0; e < nej; e++)
              centerJ[dim] += node[element[j].nodeID[e]].x[dim] / nej;

          double disTmp = (centerI - centerJ).norm();

          if(disTmp <= radius_ && element[j].isDesignable == 1)
          {
            IDi.push_back(counterj);
            disi.push_back(disTmp);
            numj++;
          }
          counterj++;
        }
      } // j loop end

      // stock filter data of element i
      numele_[counteri] = numj;
      ID_[counteri] = IDi;
      dis_[counteri] = disi;
      counteri++;
    }
  }
}

} // namespace icarat