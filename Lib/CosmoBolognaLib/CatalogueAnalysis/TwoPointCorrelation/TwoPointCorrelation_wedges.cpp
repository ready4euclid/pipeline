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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_wedges.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_multipoles used to
 *  measure the wedges of the two-point correlation
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_wedges used to measure the wedges
 *  of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_wedges.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::measure (const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{ 
  // ----------- measure the 2D two-point correlation function, xi(r,mu) ----------- 
  
  TwoPointCorrelation2D_polar::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  
  // ----------- measure the first three multipoles of the two-point correlation function ----------- 

  vector<double> rad, xiw, error;

  rad.resize(m_dd->nbins_D1()*2);
  xiw.resize(0); xiw.resize(m_dd->nbins_D1()*2, 0.);
  error.resize(0); error.resize(m_dd->nbins_D1()*2, 0.);

  vector<vector<double> > xi2d = TwoPointCorrelation2D_polar::m_dataset->fxy();
  vector<vector<double> > error_xi2d = TwoPointCorrelation2D_polar::m_dataset->error_fxy();

  double binSize = 1./m_dd->binSize_inv_D2();

  int half = index_closest(0.5,m_dd->scale_D2());

  for (int i=0; i<m_dd->nbins_D1(); i++) {

    rad[i] = m_dd->scale_D1(i);
    rad[i+m_dd->nbins_D1()] = m_dd->scale_D1(i);

    for (int j=0; j<half; j++) {
      xiw[i] += 0.5*xi2d[i][j]*binSize;   		// xi_perp
      error[i] += 0.5*pow(error_xi2d[i][j]*binSize, 2); // error[xi_perp]
    }

    for (int j=half; j<m_dd->nbins_D2(); j++) {
      xiw[i+m_dd->nbins_D1()] +=  0.5*xi2d[i][j]*binSize;                // xi_par
      error[i+m_dd->nbins_D1()] += 0.5*pow(error_xi2d[i][j]*binSize, 2); // error[xi_par]
    }  
  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv);} );
  
  m_dataset->set_xx(rad);
  m_dataset->set_fx(xiw);
  m_dataset->set_error_fx(error);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::write (const string dir, const string file, const int rank) const 
{
  vector<double> rad = m_dataset->xx();
  vector<double> xiw = m_dataset->fx();
  vector<double> error = m_dataset->error_fx();

  checkDim(rad, m_dd->nbins_D1()*3, "rad");
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);

  fout << "### rad  xi_perp  error[xi_perp]  xi_par  error[xi_par] ###" << endl;

  for (int i=0; i<m_dd->nbins_D1(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << rad[i] << "  " << setw(8) << xiw[i] << "  " << setw(8) << error[i] << "  " << setw(8) << rad[i+m_dd->nbins_D1()] << "  " << setw(8) << xiw[i+m_dd->nbins_D1()] << "  " << setw(8) << error[i+m_dd->nbins_D1()] << endl;
   
  fout.close(); cout << endl << "I wrote the file: " << file_out << endl << endl;
}
