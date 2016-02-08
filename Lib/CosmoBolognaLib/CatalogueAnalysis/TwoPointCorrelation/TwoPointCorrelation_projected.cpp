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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_projected.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_projected used to
 *  measure the projected two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_projected used to measure the projected
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_projected.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measure (const double piMax_integral, const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  if (piMax_integral>m_dd->sMax_D2()) 
    WarningMsg("Attention: the upper limit of the integral is larger than piMax -> the measured two-point correlation function will be extrapolated!");

  
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation2D_cartesian::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
 
 
  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 

  vector<double> rp, ww, error;

  vector<vector<double>> xi2d, error_xi2d;
  xi2d = TwoPointCorrelation2D_cartesian::m_dataset->fxy();
  error_xi2d = TwoPointCorrelation2D_cartesian::m_dataset->error_fxy();

  rp.resize(m_dd->nbins_D1());
  ww.resize(0); ww.resize(m_dd->nbins_D1(), 0.);
  error.resize(0); error.resize(m_dd->nbins_D1(), 0.);
 
  int pim = nint((piMax_integral-m_dd->sMin_D2())*m_dd->binSize_inv_D2()); // to convert from Mpc/h into the vector index

  double binSize = 1./m_dd->binSize_inv_D2();

  for (int i=0; i<m_dd->nbins_D1(); i++) {

    rp[i] = m_dd->scale_D1(i);
     
    for (int j=0; j<pim; j++) {

      ww[i] = 0.;
      error[i] = 0.;
      
      if (xi2d[i][j]>-1.e29) { // check!!!
	ww[i] += 2.*binSize*xi2d[i][j];
	error[i] += pow(2.*binSize*error_xi2d[i][j], 2); // check!!!!
      }
      
    }
  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv);} );

  m_dataset->set_xx(rp);
  m_dataset->set_fx(ww);
  m_dataset->set_error_fx(error);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::read(const string dir, const string file) 
{
  m_dataset->read(dir+file);
}



// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::write (const string dir, const string file, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins_D1(), "rp");
  
  m_dataset->write(dir, file, "rp", "xi_projected", rank);
}
