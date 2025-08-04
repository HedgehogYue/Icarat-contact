///  @file    node_pair_selector.cpp
///  @author	Takeshi Chang
///  @date		Nov 14, 2024.
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
}
} // namespace multibody
} // namespace icarat