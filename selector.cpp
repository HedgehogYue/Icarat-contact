///  @file    selector.cpp
///  @author	Takeshi Chang
///  @date		May 9, 2023.
#include "opt_multibody.hpp"
namespace icarat
{
namespace multibody
{
/// check every node's cooordinate and give each node a w property 1/0
void selector(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
              const toml::value &config, double &numconstnode,
              double &constvalue, int &numDesign)
{
  vector<int> designIDs, nondesignIDs;
  vector<int> outputw;
  /// non-design domain range
  double xmin = 0;
  double xmax = 0;
  double ymin = 0;
  double ymax = 0;
  double zmin = 0;
  double zmax = 0;
  /// constraint range
  double conxmin = 0;
  double conxmax = 0;
  double conymin = 0;
  double conymax = 0;
  double conzmin = 0;
  double conzmax = 0;
  double constvaluex = 0;
  double constvaluey = 0;
  double constvaluez = 0;
  /// get the non-design domain data from toml
  const auto &non = toml::find(config, "non_design");
  int nonnum = toml::find<int>(non, "non_design_num");
  const auto &range = toml::find<toml::array>(config, "non_design", "range");
  /// get the constraint data from toml
  const auto &con = toml::find(config, "constraint");
  int connum = toml::find<int>(con, "constraint_num");
  const auto &const_range =
      toml::find<toml::array>(config, "constraint", "constrange");
  const auto &const_max =
      toml::find<toml::array>(config, "constraint", "const_max");

  /// 3D problem
  if(fem.ndim == 3)
  {
    /// select the non-design domain range
    for(int nonID = 0; nonID < nonnum; nonID++)
    {
      /// get non-design domain range
      xmin = toml::get<double>(range[nonID][0][0]);
      xmax = toml::get<double>(range[nonID][0][1]);
      ymin = toml::get<double>(range[nonID][1][0]);
      ymax = toml::get<double>(range[nonID][1][1]);
      zmin = toml::get<double>(range[nonID][2][0]);
      zmax = toml::get<double>(range[nonID][2][1]);
      ///  select the nodes
      for(int i = 0; i < fem.numnp; i++)
      {
        if(xmin <= node[i].x[0] && node[i].x[0] <= xmax && //
           ymin <= node[i].x[1] && node[i].x[1] <= ymax && //
           zmin <= node[i].x[2] && node[i].x[2] <= zmax)
        {
          node[i].w = 1;
          outputw.push_back(1);
        }
        else
        {
          node[i].w = 0;
          outputw.push_back(0);
        }
      }
      /// select the elements
      for(int i = 0; i < fem.nelx; i++)
      {
        // get coordinates of cell center
        double centerX = 0.0, centerY = 0.0, centerZ = 0.0;
        for(int j = 0; j < element[i].ne; j++)
        {
          centerX += node[element[i].nodeID[j]].x[0] / element[i].ne;
          centerY += node[element[i].nodeID[j]].x[1] / element[i].ne;
          centerZ += node[element[i].nodeID[j]].x[2] / element[i].ne;
        }
        if(xmin < centerX && centerX < xmax && //
           ymin < centerY && centerY < ymax && //
           zmin < centerZ && centerZ < zmax)
        {
          nondesignIDs.push_back(i);
          element[i].isDesignable = 0;
          element[i].nonID = nonID;
        }
        else
        {
          designIDs.push_back(i);
          element[i].isDesignable = 1;
          numDesign++;
        }
      }
    } /// select the constraint range
    for(int conID = 0; conID < connum; conID++)
    {
      /// get constraint range
      conxmin = toml::get<double>(const_range[conID][0][0]);
      conxmax = toml::get<double>(const_range[conID][0][1]);
      conymin = toml::get<double>(const_range[conID][1][0]);
      conymax = toml::get<double>(const_range[conID][1][1]);
      conzmin = toml::get<double>(const_range[conID][2][0]);
      conzmax = toml::get<double>(const_range[conID][2][1]);
      constvaluex = toml::get<double>(const_max[conID][0]);
      constvaluey = toml::get<double>(const_max[conID][1]);
      constvaluez = toml::get<double>(const_max[conID][2]);
      constvalue = constvaluey;
      for(int i = 0; i < fem.numnp; i++)
      {
        if(conxmin <= node[i].x[0] && node[i].x[0] <= conxmax && //
           conymin <= node[i].x[1] && node[i].x[1] <= conymax && //
           conzmin <= node[i].x[2] && node[i].x[2] <= conzmax)
        {
          node[i].isConstraint = 1;
          node[i].constvalue = constvalue;
          numconstnode++;
        }

        else
          node[i].isConstraint = 0;
      }
    }
  }
  /// 2D problem
  else if(fem.ndim == 2)
  {
    /// select the non-design domain range
    for(int nonID = 0; nonID < nonnum; nonID++)
    {
      /// get non-design domain range
      xmin = toml::get<double>(range[nonID][0][0]);
      xmax = toml::get<double>(range[nonID][0][1]);
      ymin = toml::get<double>(range[nonID][1][0]);
      ymax = toml::get<double>(range[nonID][1][1]);
      ///  select the nodes
      for(int i = 0; i < fem.numnp; i++)
      {
        if(xmin <= node[i].x[0] && node[i].x[0] <= xmax && //
           ymin <= node[i].x[1] && node[i].x[1] <= ymax)
        {
          node[i].w = 1;
          outputw.push_back(1);
        }
        else
        {
          node[i].w = 0;
          outputw.push_back(0);
        }
      }
      /// select the elements
      for(int i = 0; i < fem.nelx; i++)
      {
        // get coordinates of cell center
        double centerX = 0.0, centerY = 0.0, centerZ = 0.0;
        for(int j = 0; j < element[i].ne; j++)
        {
          centerX += node[element[i].nodeID[j]].x[0] / element[i].ne;
          centerY += node[element[i].nodeID[j]].x[1] / element[i].ne;
        }

        if(xmin <= centerX && centerX <= xmax && ymin <= centerY &&
           centerY <= ymax)
        {
          nondesignIDs.push_back(i);
          element[i].isDesignable = 0;
          element[i].nonID = nonID;
        }
        else
        {
          designIDs.push_back(i);
          element[i].isDesignable = 1;
          numDesign++;
        }
      }
    } /// select the constraint range
    for(int conID = 0; conID < connum; conID++)
    {
      /// get constraint range
      conxmin = toml::get<double>(const_range[conID][0][0]);
      conxmax = toml::get<double>(const_range[conID][0][1]);
      conymin = toml::get<double>(const_range[conID][1][0]);
      conymax = toml::get<double>(const_range[conID][1][1]);
      constvaluex = toml::get<double>(const_max[conID][0]);
      constvaluey = toml::get<double>(const_max[conID][1]);
      constvalue = constvaluey;
      for(int i = 0; i < fem.numnp; i++)
      {
        if(conxmin <= node[i].x[0] && node[i].x[0] <= conxmax && //
           conymin <= node[i].x[1] && node[i].x[1] <= conymax)
        {
          node[i].isConstraint = 1;
          node[i].constvalue = constvalue;
          numconstnode++;
        }
        else
          node[i].isConstraint = 0;
      }
    }
  }
}

void makelocal(FEM &fem, vector<NodeMB> &node, Eigen::MatrixXd &W,
               Eigen::VectorXd &wvector)
{
  for(int i = 0; i < fem.numnp; i++)
  {
    if(node[i].dirich == nullptr)
    {
      node[i].nwd = 2;
    }
    else
    {
      node[i].nwd = 0;
      if(node[i].dirich->flag[0] == 1)
        node[i].nwd += 1;
      if(node[i].dirich->flag[1] == 1)
        node[i].nwd -= 1;
      // for(int j = 0; j < 2; j++)
      // {
      //   if(node[i].dirich->flag[j] == 0)
      //   {
      //     node[i].nwd++;
      //   }
      // }
    }
    // cout << "dirichnwd" << i << "=" << node[i].nwd << endl;
  }
  int counter;
  counter = 0;
  for(int i = 0; i < fem.numnp; i++)
  {
    if(node[i].isConstraint == 0)
    {
      for(int j = 0; j < abs(node[i].nwd); j++)
      {

        wvector(counter) = node[i].isConstraint;
        counter++;
      }
    }
    else if(node[i].isConstraint == 1)
    {
      double counterb = 0;
      for(int j = 0; j < abs(node[i].nwd); j++)
      {
        if(node[i].nwd == -1)
          wvector(counter) = 0;
        else if(node[i].nwd == 1)
          wvector(counter) = node[i].isConstraint;
        else if(node[i].nwd == 2)
        {
          wvector(counter) = counterb;
          counterb++;
        }
        counter++;
      }
    }
  }
  // cout << "w vector is" << endl << wvector << endl;
  for(int i = 0; i < fem.numeq; i++)
  {
    for(int j = 0; j < fem.numeq; j++)
    {
      if(i == j)
      {
        W(i, j) = wvector(i);
      }
      else
      {
        W(i, j) = 0;
      }
    }
  }
  // cout << "w matrix is" << endl << W << endl;
}
} // namespace multibody
} // namespace icarat