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
 *  @file Statistics/Prior.cpp
 *
 *  @brief Methods of the class Prior
 *
 *  This file contains the implementation of the methods of the class
 *  Prior
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Prior.h"
using namespace cosmobl;


// ======================================================================================


cosmobl::Prior::Prior () 
{
  m_xmin = -1.e30; 
  m_xmax = 1.e30;
  m_func = &identity<double>;  
}


// ======================================================================================


cosmobl::Prior::Prior (vector<double> Limits) 
{
  m_xmin = Limits[0]; 
  m_xmax = Limits[1];
  m_func = &identity<double>;  
}


// =====================================================================================


cosmobl::Prior::Prior (vector<double> Limits, double mean, double sigma) 
{
  m_xmin = Limits[0]; 
  m_xmax = Limits[1];
  m_prior_func_pars.push_back(mean);
  m_prior_func_pars.push_back(sigma);
  m_func = &gaussian<double>;
}


// =====================================================================================


cosmobl::Prior::Prior (double mean, double sigma) 
{
  m_xmin = -1.e30; 
  m_xmax = 1.e30;
  m_prior_func_pars.push_back(mean);
  m_prior_func_pars.push_back(sigma);
  m_func = &gaussian<double>;
}


// =====================================================================================


cosmobl::Prior::Prior (vector<double> Limits, prior_func func, vector<double> prior_func_pars) 
: m_func(func), m_prior_func_pars(prior_func_pars) 
{
   m_xmin = Limits[0]; 
   m_xmax = Limits[1];
}


// =====================================================================================


void cosmobl::Prior::set_limits (vector<double> Limits) {
  m_xmin = Limits[0]; 
  m_xmax = Limits[1];
}


// =====================================================================================


void cosmobl::Prior::set_gaussian_parameters (double mean, double sigma) 
{
  m_func = &gaussian<double>;
  m_prior_func_pars.erase(m_prior_func_pars.begin(), m_prior_func_pars.end());
  m_prior_func_pars.push_back(mean);
  m_prior_func_pars.push_back(sigma);
}


// =====================================================================================


void cosmobl::Prior::set_func_parameters (prior_func _func, vector<double> pars) 
{
   m_func = _func;
   m_prior_func_pars.erase(m_prior_func_pars.begin(),m_prior_func_pars.end());

   for (unsigned int i=0; i<pars.size(); i++)
      m_prior_func_pars.push_back(pars[i]);
}


// =====================================================================================


bool cosmobl::Prior::isIncluded (double value) {
  if (value > m_xmin && m_xmax > value)
    return 1;
  else 
    return 0;
}


