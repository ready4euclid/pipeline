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

/**
 *  @file Statistics/Parameter.cpp
 *
 *  @brief Methods of the  class Parameter 
 *
 *  This file contains the implementation of the methods of the class
 *  Parameter, used to manage model parameters in chi2 fitting or monte carlo
 *  analysis
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Parameter.h"
using namespace cosmobl;


// ============================================================================================


cosmobl::Parameter::Parameter(double value, bool freeze, string name)
{
  m_value=value;
  m_freeze = freeze;
  m_name = name;
  m_prior= make_shared<Prior>();
}


// ============================================================================================


cosmobl::Parameter::Parameter(double value, vector<double> limits, string name){
  m_value=value;
  m_freeze = 0;
  m_name = name;
  m_prior= make_shared<Prior>(limits);
}


// ============================================================================================


cosmobl::Parameter::Parameter(double value, shared_ptr<Prior> prior , string name){
  m_value=value;
  m_freeze = 0;
  m_name = name;
  m_prior = prior; 
}


// ============================================================================================


void cosmobl::Parameter::set_chains(int nchains, int chain_size){
  m_nchains = nchains;
  m_chain_size = chain_size;
  m_chains.resize(nchains);

  for (int i=0;i<m_nchains;i++){
    auto chain = make_shared<Chain>(Chain(m_chain_size));
    m_chains[i] = chain;
  }
  

}


// ============================================================================================


shared_ptr<Chain> cosmobl::Parameter::merge_chains(int max, int min, int thin){

  vector<double> values;

  for(auto &&cc : m_chains){
    int cmin = (min<=0) ? 0 : min;
    int cmax = (max<=0) ? cc->chain_size() : max;
    for(int i=cmin;i<cmax;i+=thin)
      values.push_back(cc->chain_value(i));
  }
  
  shared_ptr<Chain> chain=make_shared<Chain>(Chain(values.size()));

  for(size_t i=0;i<values.size();i++)
    chain->set_chain_value(i,values[i]);

  return chain;
}



// ============================================================================================


double cosmobl::Parameter::eval_proposed(double proposed_value){
  m_proposed_value = proposed_value;
  return eval_proposed();
}


// ============================================================================================


double cosmobl::Parameter::eval_proposed(){
  return PriorProbability(m_proposed_value)/PriorProbability();
}


// ============================================================================================


void cosmobl::Parameter::confirm_proposed_value(){
  m_value = m_proposed_value;
}


// ============================================================================================


double cosmobl::Parameter::random_value(double random)
{
  m_value = (m_prior->xmax()-m_prior->xmin())*random+m_prior->xmin();
  return m_value;
}
