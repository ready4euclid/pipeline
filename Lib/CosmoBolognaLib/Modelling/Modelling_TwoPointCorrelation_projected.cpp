/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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

/**
 *  @file Modelling/Modelling_TwoPointCorrelation_projected.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation_projected, 
 *  used for modelling projected 2pcf
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_projected
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "Modelling_TwoPointCorrelation_projected.h"

using namespace cosmobl;


// ============================================================================================


cosmobl::modelling::Modelling_TwoPointCorrelation_projected::Modelling_TwoPointCorrelation_projected(const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop, const double redshift, const Cosmology cosmology)
{

  m_data = twop->dataset();
  m_cosmology = make_shared<Cosmology>(cosmology);
  m_redshift = redshift;
  m_twoPType = twopt::TwoPType::_1D_projected_;

}


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation_projected::fit_bias (const string LikelihoodType, const vector<double> xlimits, const double bias_value, const statistics::Prior bias_prior, const int nChains, const int chain_size, const string dir_output, const double start, const double stop, const int thin)
{

  m_data->set_limits(xlimits[0],xlimits[1]);

  ModelBias model(bias_value,bias_prior);
  double pimax = 0; //WORK IN PROGRESS
  model.set_wp_parameters(m_data->xx(),m_cosmology,m_redshift, pimax);
  m_model = make_shared<ModelBias>(model);
  
  statistics::likelihood_npar func;
  if(LikelihoodType == "model")
    func = &statistics::likelihood_gaussian_1D_model_npar;
  else if(LikelihoodType == "error")
    func = &statistics::likelihood_gaussian_1D_error_npar;
  else if(LikelihoodType == "covariance")
    func = &statistics::likelihood_gaussian_1D_covariance_npar;
  else
    ErrorMsg("Error in fit_bias of Modelling_TwoPointCorrelation1D_projected, no such type of likelihood");

  statistics::Likelihood lik(m_data,m_model,func);
  lik.sample(nChains, chain_size); 
  string dir_output_file=dir_output+"projected_bias_xmin="+conv(m_data->x_down(),par::fDP1)+"_xmax="+conv(m_data->x_up(),par::fDP1)+"_"+LikelihoodType;
  lik.write_chain(dir_output_file, start, stop, thin);

}


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation_projected::fit_bias_cosmology(const string LikelihoodType, const vector<double> xlimits, const double bias_value, const statistics::Prior bias_prior, const vector<CosmoPar> CosmoPars, const vector<statistics::Prior> prior_CosmoPars, const int nChains, const int chain_size, const string dir_output, const double start, const double stop, const int thin){}
