// =====================================
// Example code: how to compute the chi2
// =====================================

#include "GlobalFunc.h"

using namespace cosmobl;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;

int main () {

  string input_data = "../input/gaussian.dat";
  double xmin = -5, xmax = 5;

  auto data = Data::Create(input_data,xmin,xmax);

  vector<shared_ptr<Parameter>> parameters(3);
  vector<double> pars = {0, 1, 10000};

  parameters[0] = make_shared<Parameter>(Parameter(0., 0, "Mean"));
  parameters[1] = make_shared<Parameter>(Parameter(1, 0, "Sigma"));
  parameters[2] = make_shared<Parameter>(Parameter(10000, 0, "Norm"));

  parameters[0]->set_prior({-1, 1});
  parameters[1]->set_prior({0.01, 2});
  parameters[2]->set_prior({9000, 11000});

  auto model = make_shared<Model1D>(Model1D(parameters, NULL, &gaussian<double>));

  Chi2 chi2(data, model);
  chi2.minimize(pars, "error", 1, 1000);

  cout << parameters[0]->name() << " = " << parameters[0]->value() << endl;
  cout << parameters[1]->name() << " = " << parameters[1]->value() << endl;
  cout << parameters[2]->name() << " = " << parameters[2]->value() << endl;

  return 0;
} 
