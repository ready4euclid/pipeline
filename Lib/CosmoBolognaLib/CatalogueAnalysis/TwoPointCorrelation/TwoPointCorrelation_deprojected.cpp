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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_deprojected.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_deprojected used to
 *  measure the projected two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_deprojected used to measure the projected
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_deprojected.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measure (const double piMax_integral, const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  
  // ----------- measure the projected two-point correlation function, w(rp) ----------- 

  TwoPointCorrelation_projected::measure(piMax_integral, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
  
  
  // ----------- estimate the deprojected (real-space) two-point correlation function -----------
  // ------------------ (following Saunders, Rowan-Robinson and Lawrence 1992) ------------------ 

  vector<double> rad, xi, error_xi;

  rad.resize(m_dd->nbins_D1());
  xi.resize(0); xi.resize(m_dd->nbins_D1(), 0.);
  error_xi.resize(0); error_xi.resize(m_dd->nbins_D1(), 0.);

  vector<double> rp = TwoPointCorrelation_projected::m_dataset->xx();
  vector<double> ww = TwoPointCorrelation_projected::m_dataset->fx();
  vector<double> error_ww = TwoPointCorrelation_projected::m_dataset->error_fx();

  for (int i=0; i<m_dd->nbins_D1(); i++) {
    
    double ri = m_dd->scale_D1(i);
    rad[i] = ri;

    for (int j=i; j<m_dd->nbins_D2()-1; j++) {
      double rj = m_dd->scale_D1(j);
      double rj1 = m_dd->scale_D1(j+1);
      double fact = 1./(rj1-rj)*log((rj1+sqrt(max(0., rj1*rj1-ri*ri)))/(rj+sqrt(max(0., rj*rj-ri*ri))));

      xi[i] -= (ww[j+1]-ww[j])*fact;
      error_xi[i] += pow(error_ww[j+1]*fact, 2)+pow(error_ww[j]*fact, 2);
    }
    
  }

  for_each( xi.begin(), xi.end(), [] (double &vv) { vv /= par::pi;} );
  for_each( error_xi.begin(), error_xi.end(), [] (double &vv) { vv = sqrt(vv/par::pi);} );

  m_dataset->set_xx(rad);
  m_dataset->set_fx(xi);
  m_dataset->set_error_fx(error_xi);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::write (const string dir, const string file, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins_D1(), "rad");
  m_dataset->write(dir, file, "rad", "xi_deprojected", rank);
}
