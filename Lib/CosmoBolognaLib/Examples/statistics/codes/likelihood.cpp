// ===========================================
// Example code: how to construct a likelihood
// ===========================================

#include "GlobalFunc.h"

using namespace cosmobl;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;

double glikelihood(vector<double> parameters, shared_ptr<void> fixed_parameters)
{
  return exp(-0.5*chi2_1d_error_npar(parameters,fixed_parameters));
}

int main () {

  string input_data = "../input/gaussian.dat";
  double xmin = -5, xmax = 5;

  auto data = Data::Create(input_data,xmin,xmax);

  vector<shared_ptr<Parameter>> parameters(3);
  vector<double> pars = {0, 1, 10000};

  parameters[0] = make_shared<Parameter>(Parameter(0., 1, "Mean"));
  parameters[1] = make_shared<Parameter>(Parameter(1, 0, "Sigma"));
  parameters[2] = make_shared<Parameter>(Parameter(10000, 0, "Norm"));

  parameters[0]->set_prior({-1, 1});
  parameters[1]->set_prior({0.01, 2});
  parameters[2]->set_prior({9000, 11000});

  auto model = make_shared<Model1D>(Model1D(parameters, NULL, &gaussian<double>));

  Likelihood likelihood(data, model, glikelihood);
  
  int nchains = 1000;
  int chain_size = 1000;
  likelihood.sample(nchains, chain_size);

  string output_file = "../output/chains_npar";
  double start = 0.5;
  double end = 1;
  int thin = 10;
  likelihood.write_chain(output_file, start, end, thin);

  return 0;
} 
