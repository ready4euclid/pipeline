/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/** @file Statistics/Likelihood.cpp
 *
 *  @brief Methods of the class Likelihood 
 *
 *  This file contains the implementation of the methods of the class
 *  Likelihood, used for bayesian analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Likelihood.h"
using namespace cosmobl;


// ============================================================================================


void cosmobl::Likelihood::minimize(double parameter, unsigned int max_iter, double min, double max){

  unsigned int npar = m_model->npar_eff();
  if(npar!=1)
    ErrorMsg("Error in minimize of Chi2, wrong number of parameters");

  if(!m_model->parameter(0)->prior()->isIncluded(parameter)){ 
    string msg = "Warning, starting value for parameter '"+m_model->parameter(0)->name()+"' is out of range, it will be changed";
    WarningMsg(msg);
    parameter=0.5*(m_model->parameter(0)->prior()->xmax()+m_model->parameter(0)->prior()->xmin());
  }

  min = (min>-1.e29) ? min : m_model->parameter(0)->prior()->xmin();
  max = (max>1.e29) ? max : m_model->parameter(0)->prior()->xmax();

  auto fixed_parameters = make_shared<STR_params>(STR_params(m_data,m_model));

  auto ff = [&](double x, shared_ptr<void> p){ return -2*Log(m_likelihood_1par(x,p));};

  GSLfunction_1D_1 func(ff,fixed_parameters);
  func.minimize(parameter,max_iter,min,max);
  m_model->update_parameters(parameter);
  cout << "Done" << endl;
  
}


// ============================================================================================


void cosmobl::Likelihood::minimize(vector<double> parameters, unsigned int max_iter, double tol){

  unsigned int npar = m_model->npar_eff();
  if(npar==0)
    ErrorMsg("Error in minimize of Chi2, there is no parameter to vary");

  auto fixed_parameters = make_shared<STR_params>(STR_params(m_data,m_model));

  vector<double> pars ;
  if(parameters.size() == npar){
    pars = parameters;
  }
  else if(parameters.size() == m_model->npar()){
    for(unsigned int i=0;i<m_model->npar();i++){
      if(!m_model->parameter(i)->freeze())
        pars.push_back(parameters[i]);
    }
  }
  else{ErrorMsg("Error in minimize of Chi2, unrecognized number of parameters");}

  vector<double> step_size(npar,1);
  int nn=0;
  for(unsigned int i=0;i<m_model->npar();i++){
    if (!m_model->parameter(i)->freeze()) 
    { 
      if(!m_model->parameter(i)->prior()->isIncluded(pars[nn])){ 
        string msg = "Warning, starting value for parameter '"+m_model->parameter(i)->name()+"' is out of range, it will be changed";
        WarningMsg(msg);
        pars[nn]=0.5*(m_model->parameter(i)->prior()->xmax()+m_model->parameter(i)->prior()->xmin());
      }
    
      if ( m_model->parameter(i)->interval_size()<1.e29){
        step_size[nn]=1.e-1*m_model->parameter(i)->interval_size();
      }
      nn++;
    }
  }
  auto ff = [&](vector<double> x, shared_ptr<void> p){ return -2*Log(m_likelihood_npar(x,p));};

  GSLfunction_nD_1 func(npar,ff,fixed_parameters);
  func.minimize(pars,step_size,max_iter,tol);
  m_model->update_parameters(pars);
  cout << "Done" << endl;
}


// ============================================================================================


double cosmobl::Likelihood::sample(int nchains, int chain_size){
  
  vector< shared_ptr<Parameter> > pp = m_model->parameters();
  m_nchains=nchains;
  m_chain_size=chain_size;

  vector<double> val_likelihood(m_nchains,0);
  vector<double> accept(m_nchains,0);

  auto fixed_parameters = make_shared<STR_params>(STR_params(m_data,m_model));

  for(size_t i=0;i<pp.size();i++)
    pp[i]->set_chains(m_nchains, m_chain_size);

  default_random_engine generator;
  generator.seed(time(NULL));
  uniform_real_distribution<double> distribution(0.,1.);
  uniform_int_distribution<int> idist(0.,m_nchains-1);

  auto random_number = bind(distribution,generator);
  auto kchain = bind(idist,generator);

  double gzpar=2;
  double zmin = 1./gzpar;
  double zmax = gzpar;

  uniform_real_distribution<double> gz(zmin,zmax);
  vector<double> gzz;
  vector<double> gen_z;
  vector<double> prior_probs(m_nchains,1);
  
  vector<double> zz = linear_bin_vector(1000,zmin,zmax);
  for(auto &&zzz : zz)
    gzz.push_back(1./sqrt(zzz));

  // Initialite chains
  for (int i=0;i<m_nchains;i++){
    vector<double> starting_parameters(m_model->npar());
    for(size_t j=0;j<pp.size();j++){
      starting_parameters[j] = (pp[j]->freeze()) ? pp[j]->value() : pp[j]->random_value(random_number());
      pp[j]->set_value(starting_parameters[j]);
      pp[j]->chain(i)->set_chain_value(0,starting_parameters[j]);
      prior_probs[i] *= pp[j]->PriorProbability(starting_parameters[j]);
//      cout << pp[j]->value() << " " << pp[j]->freeze() << " " ;
    }
//    cout << endl;
    val_likelihood[i] = m_likelihood_npar(starting_parameters,fixed_parameters);
  }

  //Populate chains
  for(int n=1;n<m_chain_size;n++){
    vector<double> generated_z;
    fill_distr (m_nchains, zz, gzz, generated_z, zmin, zmax, n);

    for(int i=0;i<m_nchains;i++){

      double prob=1;
      vector<double> proposed(pp.size(),0);
      int kk=i;
      while(kk==i)
	kk=kchain();

      for(size_t j=0;j<pp.size();j++){
	if(!pp[j]->freeze()){
	  proposed[j] =  pp[j]->chain(kk)->chain_value(n-1)+ generated_z[i]*(pp[j]->chain(i)->chain_value(n-1)-pp[j]->chain(kk)->chain_value(n-1));
	}
	else {proposed[j] = pp[j]->chain(i)->chain_value(n-1);}
	prob*=pp[j]->PriorProbability(proposed[j]);
	if(prob>1.e29){prob=0;}
      }
      prob = prob/prior_probs[i];


      double proposed_likelihood = m_likelihood_npar(proposed,fixed_parameters);
      prob *= min(1. , pow(generated_z[i],m_model->npar()-1)*proposed_likelihood/val_likelihood[i] );

      if (random_number() < prob){
	prior_probs[i]=1;
	for(size_t j=0;j<pp.size();j++){
	  pp[j]->chain(i)->set_chain_value(n,proposed[j]);
	  prior_probs[i]*=pp[j]->PriorProbability(proposed[j]);
	}
	val_likelihood[i] = proposed_likelihood;
	accept[i]+=1;
      }
      else {
	for(size_t j=0;j<pp.size();j++)
	  pp[j]->chain(i)->set_chain_value(n,pp[j]->chain(i)->chain_value(n-1));
      }

    }
  }

  cout << "mean acceptance ratio =" << Average(accept)/m_chain_size << endl;

  for(size_t i = 0; i< pp.size();i++){
    if(pp[i]->freeze()) { cout << "Parameter " << i << " freezed" << endl;}
    else{
      auto chain = pp[i]->merge_chains(m_chain_size,m_chain_size*0.5);
      chain->Statistics();
      cout << chain->chain_size() << " " << chain->mean() << " " << chain->std() << " " << chain->median() << endl;
      pp[i]->set_value(chain->mean());

      vector<double> mean,var;
      for(auto &&cc : pp[i]->chains()){
	cc->Statistics(m_chain_size,0.5*m_chain_size);
	mean.push_back(cc->mean());
	var.push_back(cc->std()*cc->std());
      }
      double W  = Average(var);
      double B = m_chain_size*Sigma(mean)*Sigma(mean);
      double R = sqrt((double(m_chain_size-1)/m_chain_size*W+B/m_chain_size)/W);
      cout <<"sqrt(R) = " << R << endl;
    }
  }

  return Average(accept)/m_chain_size;
  
}


// ============================================================================================


double cosmobl::Likelihood::sample(int nchains, int chain_size, double shift, int nsubstep){
  
  shared_ptr<Parameter> pp = m_model->parameter(0);
  m_nchains=nchains;
  m_chain_size=chain_size;
  double par_shift = 0.01*shift*pp->interval_size();
  vector<double> val_likelihood(m_nchains,0);
  vector<double> accept(m_nchains,0);

  auto fixed_parameters = make_shared<STR_params>(STR_params(m_data,m_model));
  pp->set_chains(m_nchains, m_chain_size);

  default_random_engine generator(time(NULL));
  uniform_real_distribution<double> distribution(0.,1.);
  auto random_number = bind(distribution,generator);
  normal_distribution<double> ndistribution(0,par_shift);

  for(int i=0;i<m_nchains;i++){
    double starting_point = pp->random_value(random_number());
    pp->set_value(starting_point);
    pp->chain(i)->set_chain_value(0,starting_point);
    val_likelihood[i]=m_likelihood_1par(starting_point,fixed_parameters);
    vector<double> acc(m_chain_size,0);

    for(int n=1;n<m_chain_size;n++){
      double aa=0;
      for(int nsub = 0;nsub <nsubstep;nsub++){
	// propose //
	double proposed = pp->value()+ndistribution(generator);
	double vlikelihood = m_likelihood_1par(proposed,fixed_parameters);
	double prob = pp->eval_proposed(proposed);
	if(prob>1.e29){prob=0;}
	prob=min(1.,prob*vlikelihood/val_likelihood[i]);

	if(random_number()<prob){
	  pp->confirm_proposed_value();
	  pp->chain(i)->set_chain_value(n,proposed);
	  val_likelihood[i] = vlikelihood;
	  aa += 1.;
	}
      }
      if(aa==0) {pp->chain(i)->set_chain_value(n,pp->value());};
      acc[n]=aa/nsubstep;
    }
    accept[i]=Average(acc);
  }
  print(accept);

  cout << "mean acceptance ratio =" << Average(accept) << endl;

  auto chain = pp->merge_chains(m_chain_size,m_chain_size*0.5);
  chain->Statistics();
  cout << chain->chain_size() << " " << chain->mean() << " " << chain->std() << " " << chain->median() << endl;
  pp->set_value(chain->mean());
  vector<double> mean,var;
  for(auto &&cc : pp->chains()){
    cc->Statistics(m_chain_size,m_chain_size*0.5);
    mean.push_back(cc->mean());
    var.push_back(cc->std()*cc->std());
  }
  double w  = Average(var);
  double b = m_chain_size*Sigma(mean)*Sigma(mean);
  double r = sqrt((double(m_chain_size-1)/m_chain_size*w+b/m_chain_size)/w);
  cout <<"sqrt(r) = " << r << endl;
  return Average(accept);
}

void cosmobl::Likelihood::write_chain(string output_file, double start, double stop, int thin){

  auto fixed_parameters = make_shared<STR_params>(STR_params(m_data,m_model));

  ofstream fout(output_file.c_str());
  start = int(m_nchains*start);
  stop = int(m_nchains*stop);
  int n=0;

  if(m_npar){
    vector< shared_ptr<Parameter> > pp = m_model->parameters();
    for(int j=start;j<stop;j+=thin){
      for(int i=0;i<m_nchains;i++){
	fout << n << " ";
	vector<double> pars;
	for(size_t k = 0;k<pp.size();k++){
	  pars.push_back(pp[k]->chain(i)->chain_value(j));
	  fout << " " << pp[k]->chain(i)->chain_value(j) << " ";
	}
	fout << -2*Log(m_likelihood_npar(pars,fixed_parameters)) << endl;
	n++;
      }
    }
  }
  else{
    shared_ptr<Parameter> pp = m_model->parameter(0);
    for(int j=start;j<stop;j+=thin){
      for(int i=0;i<m_nchains;i++){
	fout << i << " " << pp->chain(i)->chain_value(j) << " " <<  -2*Log(m_likelihood_1par(pp->chain(i)->chain_value(j),fixed_parameters)) <<endl;
      }
    }
  }

  fout.clear(); fout.close();
}
