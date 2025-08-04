#include "opt_multibody.hpp"

using namespace std;
using namespace icarat;
using namespace icarat::multibody;
using namespace Eigen;

void output(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
            Optimize &optim);
// void selectelements(FEM &fem, vector<ElementDD> &element, vector<NodeDD>
// &node,
//                     ReadCFG &cfg);
void checker(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file> <.mdpa file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  double V0 = 0.0;     // first volume
  double T = 0.5;
  double xmin = 1.0e-4;
  double xmax = 1.0;
  double C = 0.0;

  // set input file
  const auto config = toml::parse(argv[1]);
  // vector<string> inputMDPA = read_file(argv[2]);
  string vtk = "output2400.vtk";
  // set filter radius
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  // set number of constraint function
  int numConst = toml::find<int>(config, "this_problem", "numConst");
  // set the type of extra constraint function
  int typeConst = toml::find<int>(config, "this_problem", "typeConst");

  // set problem
  FEM fem;
  // MeshIcarat<ElementDD, NodeDD> mesh(fem, cfg);
  MeshVTK<ElementMB, NodeMB> mesh(fem, vtk, config);
  mesh.setParameter("structure");

  // set node
  vector<NodeMB> node(fem.numnp);
  // set element(unstructured mesh version)
  vector<ElementMB> element;
  // set element
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      ElementMB actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof);
      // actele.setParameter(cfg);
      element.push_back(actele);
    }
  }

  // input mesh to element & node
  mesh.generate(element, node);

  // get pair of design area's ID & non-design area's ID
  double numconstnode = 0;
  double constvalue = 0;
  int numDesign = 0;
  selector(fem, element, node, config, numconstnode, constvalue, numDesign);
  double allconst = constvalue * constvalue * numconstnode;
  cout << "all const=" << allconst << endl;
  cout << "num const=" << numConst << endl;
  cout << "type const=" << typeConst << endl;

  // for(int i = 0; i < fem.nelx; i++)
  //   cout << element[i].isDesignable << endl;
  // exit(1);

  for(int i = 0; i < fem.nelx; i++)
  {
    element[i].setParameter(config);
  }

  // boundary condition class
  BCIcarat<ElementMB, NodeMB> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // localizer W matrix
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(fem.numeq, fem.numeq);
  Eigen::VectorXd wvector = Eigen::VectorXd::Zero(fem.numeq);
  // get the W matrix and vector
  makelocal(fem, node, W, wvector);

  // checker(fem, element, node);
  // set optimization's classes
  // Optimize optim(fem.nelx, numConst, xmin, xmax);
  // MMASolver mma(fem.nelx, numConst);
  Optimize optim(numDesign, numConst, xmin, xmax);
  MMASolver mma(numDesign, numConst);
  FilterMultibody<ElementMB, NodeMB> filter(fem, element, node, filterR,
                                            numDesign);

  int numcount = 0;
  // initialization
  optim.optstep = 0;
  // for(int i = 0; i < fem.nelx; i++)
  //   optim.design_s[i] = element[i].design_s();
  for(int i = 0; i < fem.nelx; i++)
  {
    if(element[i].isDesignable == 1)
    {
      optim.design_s[numcount] = element[i].design_s();
      numcount++;
    }
  }
  if(numcount != numDesign)
  {
    cout << "something wrong in counting the non design number" << endl;
    exit(1);
  }

  LinearAnalysisMB<ElementMB, NodeMB> compt(fem, element, node, force, dis);

  //  convergence = 1;
  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    auto design_filtered = filter.densityFilter(optim.design_s, false);

    /// set threshold function
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    cout << "beta: " << beta << endl;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDeriv = filter.thresholdDerivative(design_filtered, beta, T);

    // set threshold design variable
    numcount = 0;
    for(int i = 0; i < fem.nelx; i++)
    {
      if(element[i].isDesignable == 1)
      {
        // this is two filter version
        // element[i].setDesignS(design_threshold[numcount]);
        // this is only density filter version
        element[i].setDesignS(design_filtered[numcount]);
        // before doing FDM, please change this part into design_s
        // element[i].setDesignS(optim.design_s[numcount]);
        numcount++;
      }
    }

    // structural analysis
    compt.solve("Eigen_CG", 1);

    // Eigen::VectorXd disp = Eigen::VectorXd::Zero(fem.numeq);

    /// get  stress
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get objective function
    // here！！！！
    optim.object_f = dis.val.transpose() * compt.Knon() * dis.val;
    // optim.object_f = dis.val.transpose() * compt.Kdesign() * dis.val;
    cout << "function value = " << optim.object_f << endl;

    // get const_func
    if(numConst == 1)
    {
      volume_constraint(fem, optim, element, V0);
    }
    else if(numConst == 2)
    {
      volume_constraint(fem, optim, element, V0);
      if(typeConst == 1)
      {
        displacement_constraint(fem, optim, dis, wvector, allconst);
      }
      else if(typeConst == 2)
      {
        compliance_constraint(fem, optim, dis, element, compt, config, C);
      }
      else
      {
        cout << "there is a wrong constraint function type" << endl;
        exit(1);
      }
    }
    else
    {
      cout << "there is a wrong constraint function number" << endl;
      exit(1);
    }

    vector<double> dg1ds(numDesign), dg2ds(numDesign);

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    addPlot(optim.optstep, optim.const_h[0], "constraint1_value.csv",
            optim.optstep);
    addPlot(optim.optstep, optim.const_h[1], "constraint2_value.csv",
            optim.optstep);
    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-7, 50);

    // analytical
    Eigen::VectorXd Evector = Eigen::VectorXd::Zero(fem.nelx);
    calculate_dKds(fem, element, Evector);
    // sensitivity(fem, element, node, optim, dis, compt, Evector);
    objective_sensitivity(fem, element, node, optim, dis, compt, Evector);
    if(numConst == 1)
    {
      volume_constraint_sensitivity(fem, element, optim, V0, dg1ds);
    }
    else if(numConst == 2)
    {
      volume_constraint_sensitivity(fem, element, optim, V0, dg1ds);
      if(typeConst == 1)
      {
        displacement_constraint_sensitivity(fem, element, node, optim, allconst,
                                            dis, compt, wvector, dg2ds,
                                            Evector);
      }
      else if(typeConst == 2)
      {
        compliance_constraint_sensitivity(fem, element, node, optim, dis, compt,
                                          dg2ds, Evector, C);
      }
      else
      {
        cout << "there is a wrong constraint function type" << endl;
        exit(1);
      }
    }
    else
    {
      cout << "there is a wrong constraint function number" << endl;
      exit(1);
    }
    // binary filter
    // if need to cancel the threshold filter, change here
    // for(int i = 0; i < numDesign; i++)
    // {
    //   optim.dfds[i] *= threDeriv[i];
    //   dg1ds[i] *= threDeriv[i];
    //   dg2ds[i] *= threDeriv[i];
    // }

    // filtering sensitivity
    optim.dfds = filter.densityFilter(optim.dfds, true);
    dg1ds = filter.densityFilter(dg1ds, true);
    dg2ds = filter.densityFilter(dg2ds, true);

    // FDM
    if(toml::find<bool>(config, "this_problem", "fdm"))
    {
      fdm(fem, element, node, force, dis, optim, filter, T);
      if(numConst == 2)
      {
        if(typeConst == 1)
        {
          fdm_constraint(fem, element, node, force, dis, optim, filter, T,
                         wvector, allconst);
        }
        else if(typeConst == 2)
        {
          fdm_compliance_constraint(fem, element, node, force, dis, optim,
                                    filter, T, C);
        }
      }
    }

    // for(int i = 0; i < fem.nelx; i++)
    for(int i = 0; i < numDesign; i++)
    {
      if(numConst == 1)
      {
        optim.dgds[i] = dg1ds[i];
      }
      else if(numConst == 2)
      {
        optim.dgds[numConst * i] = dg1ds[i];
        optim.dgds[numConst * i + 1] = dg2ds[i];
      }
    }

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 0;
    }
  }
  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node,
            Optimize &optim)
{
  ofstream fout("res" + to_string(optim.optstep) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](NodeMB &n) { return n.val; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](ElementMB &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](ElementMB &e) { return e.mises(); });

  addElementScalarVTK("young", element, fout, false,
                      [&](ElementMB &e) { return e.young(); });

  addElementScalarVTK("isdesignable", element, fout, false,
                      [&](ElementMB &e) { return e.isDesignable; });

  addElementScalarVTK("strainenergy", element, fout, false,
                      [&](ElementMB &e) { return e.strainenergy(); });

  fout.close();
}

/// make the w matrix from the nwd for the global calculation
/// the size of w matrix is the same as global K

// checker to check the parameter of elements or nodes
void checker(FEM &fem, vector<ElementMB> &element, vector<NodeMB> &node)
{
  for(int i = 0; i < fem.nelx; i++)
  {
    cout << "element No." << element[i].ID << endl;
    cout << "element isDesign" << element[i].isDesignable << endl;
    cout << "element young" << element[i].young2() << endl;
    cout << "element poisson" << element[i].poisson() << endl;
    cout << "element nonid" << element[i].nonID << endl;
  }
}
