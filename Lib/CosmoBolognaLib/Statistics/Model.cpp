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
 *  @file Statistics/Model.cpp
 *
 *  @brief Methods of the class Model
 *
 *  This file contains the implementation of the methods of the class
 *  Model
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Model.h"
using namespace cosmobl;


// ======================================================================================


cosmobl::Model::Model (const vector<shared_ptr<Parameter> > parameters, const shared_ptr<void> model_parameters)
  : m_parameters(parameters), m_model_parameters(model_parameters)
{
  m_npar = m_parameters.size();
  m_npar_eff = npar_eff();
}


// ======================================================================================


double cosmobl::Model::update_parameters (const double new_parameter)
{
  double parameter;

  if (!m_parameters[0]->freeze()) {
    parameter = (m_parameters[0]->prior()->isIncluded(new_parameter)) ? new_parameter : closest(new_parameter,m_parameters[0]->prior()->xmin(),m_parameters[0]->prior()->xmax());
    m_parameters[0]->set_value(parameter);
  }
  
  else
    parameter = m_parameters[0]->value();
  
  m_parameters[0]->set_value(parameter);

  return parameter;
}


// ======================================================================================


vector<double> cosmobl::Model::update_parameters (const vector<double> new_parameters)
{  
  vector<double> parameters;
  int nn = 0;
  
  for (unsigned int i=0; i<m_npar; i++) {
    
    if (!m_parameters[i]->freeze()) {
      double value = (m_parameters[i]->prior()->isIncluded(new_parameters[nn])) ? new_parameters[nn] : closest(new_parameters[nn], m_parameters[i]->prior()->xmin(), m_parameters[i]->prior()->xmax());
      m_parameters[i]->set_value(value);
      nn ++;
    }

    parameters.push_back(m_parameters[i]->value());
  }
  
  return parameters;
}

