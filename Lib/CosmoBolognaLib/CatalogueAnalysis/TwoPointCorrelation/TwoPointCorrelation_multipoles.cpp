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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_multipoles.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_multipoles used to
 *  measure the first three multipoles of the two-point correlation
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_multipoles used to measure the first three
 *  multipoles of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_multipoles.h"
#include "Data1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation_multipoles::MultipolesTwoP(const vector<double> rr, const vector<double> mu, const vector<vector<double> > xi, const vector<vector<double> > error_xi)
{
  vector<double> rad, xil, error;

  rad.resize(rr.size()*3);
  xil.resize(0); xil.resize(rr.size()*3, 0.);
  error.resize(0); error.resize(rr.size()*3, 0.);

  double binSize = rr[1]-rr[0];

  for (size_t i=0; i<rr.size(); i++) {

    rad[i] = rr[i];
    rad[i+rr.size()] = rr[i];
    rad[i+2*rr.size()] = rr[i];

    for (size_t j=0; j<mu.size(); j++) {
      double mmu = mu[j];

      xil[i]             += xi[i][j]*binSize;            // xi_0
      xil[i+rr.size()]   += 5.*xi[i][j]*P_2(mmu)*binSize; // xi_2
      xil[i+2*rr.size()] += 9.*xi[i][j]*P_4(mmu)*binSize; // xi_4

      error[i]           += pow(error_xi[i][j]*binSize, 2);            // error[xi_0]
      error[i+rr.size()] += pow(5.*error_xi[i][j]*P_2(mmu)*binSize, 2); // error[xi_2]
      error[i+2*rr.size()] += pow(9.*error_xi[i][j]*P_4(mmu)*binSize, 2); // error[xi_4]
    }

  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv); } );
  return move(unique_ptr<Data1D>(new Data1D(rad,xil,error)));

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_multipoles::measure(const ErrorType errType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int nMocks, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  switch(errType){
    case(ErrorType::_Poisson_):
      measurePoisson(dir_output_pairs,dir_input_pairs,count_dd,count_rr,count_dr,tcount);
      break;
    case(ErrorType::_Jackknife_):
      measureJackknife(dir_output_pairs,dir_input_pairs,dir_output_ResampleXi,count_dd,count_rr,count_dr,tcount);
      break;
    case(ErrorType::_Bootstrap_):
      measureBootstrap(nMocks, dir_output_pairs,dir_input_pairs,dir_output_ResampleXi,count_dd,count_rr,count_dr,tcount);
      break;
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_multipoles::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation2D_polar::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 
  
  m_dataset = MultipolesTwoP(TwoPointCorrelation2D_polar::xx(),TwoPointCorrelation2D_polar::yy(),TwoPointCorrelation2D_polar::xi2D(),TwoPointCorrelation2D_polar::error2D());

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_multipoles::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{

  if (dir_output_ResampleXi != par::defaultString){
    string mkdir = "mkdir -p "+dir_output_ResampleXi;
    if(system(mkdir.c_str())){}
  }

  vector<shared_ptr<Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_polar_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount);

  
  auto data_polar = (count_dr>=0) ? LandySzalayEstimatorTwoP(m_dd,m_rr,m_dr,m_data->weightedN(),m_random->weightedN()) : NaturalEstimatorTwoP(m_dd,m_rr,m_data->weightedN(),m_random->weightedN());

  if (count_dr==1) 
    data = XiJackknife(dd_regions, rr_regions,dr_regions);
  else
    data = XiJackknife(dd_regions, rr_regions);

  vector<vector<double> > ww,covariance;
  for(size_t i=0;i<data.size();i++){
    ww.push_back(data[i]->fx());
    if (dir_output_ResampleXi != par::defaultString){
      string filename = dir_output_ResampleXi+"xi_multipoles_Bootstrap_"+conv(i,par::fINT);
      ofstream fout(filename.c_str());
      fout.clear(); fout.close();

    }
  }

  covariance_matrix(ww,covariance,1);
  m_dataset = MultipolesTwoP(data_polar->xx(),data_polar->yy(),data_polar->fxy(),data_polar->error_fxy());
  m_dataset->set_covariance_fx(covariance);

}

// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_multipoles::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  if (dir_output_ResampleXi != par::defaultString){
    string mkdir = "mkdir -p "+dir_output_ResampleXi;
    if(system(mkdir.c_str())){}
  }
  vector<shared_ptr<Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_polar_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount);

  auto data_polar = (count_dr>=0) ? LandySzalayEstimatorTwoP(m_dd,m_rr,m_dr,m_data->weightedN(),m_random->weightedN()) : NaturalEstimatorTwoP(m_dd,m_rr,m_data->weightedN(),m_random->weightedN());

  if (count_dr==1) 
    data = XiBootstrap(nMocks, dd_regions, rr_regions,dr_regions);
  else
    data = XiBootstrap(nMocks, dd_regions, rr_regions);

  vector<vector<double> > ww,covariance;
  for(size_t i=0;i<data.size();i++){
    ww.push_back(data[i]->fx());
    if (dir_output_ResampleXi != par::defaultString){
      string filename = dir_output_ResampleXi+"xi_multipoles_Bootstrap_"+conv(i,par::fINT);
      ofstream fout(filename.c_str());
      fout.clear(); fout.close();

    }
  }
  covariance_matrix(ww,covariance,0);

  m_dataset = MultipolesTwoP(data_polar->xx(),data_polar->yy(),data_polar->fxy(),data_polar->error_fxy());
  m_dataset->set_covariance_fx(covariance);
}

// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_multipoles::XiJackknife(const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<Data> > data;
  auto data2d = TwoPointCorrelation2D_polar::XiJackknife(dd,rr);
  for (size_t i=0;i<data2d.size();i++){
    data.push_back(move(MultipolesTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  }
  return data;
}

// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_multipoles::XiJackknife(const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<Data> > data;
  auto data2d = TwoPointCorrelation2D_polar::XiJackknife(dd,rr,dr);
  for (size_t i=0;i<data2d.size();i++){
    data.push_back(move(MultipolesTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  }
  return data;
}

// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_multipoles::XiBootstrap(const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<Data> > data;
  auto data2d = TwoPointCorrelation2D_polar::XiBootstrap(nMocks,dd,rr);
  for (size_t i=0;i<data2d.size();i++){
    data.push_back(move(MultipolesTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  }
  return data;
}

// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_multipoles::XiBootstrap(const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<Data> > data;
  auto data2d = TwoPointCorrelation2D_polar::XiBootstrap(nMocks,dd,rr,dr);
  for (size_t i=0;i<data2d.size();i++){
    data.push_back(move(MultipolesTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  }
  return data;
}
/*
void cosmobl::twopt::TwoPointCorrelation_multipoles::measure (const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{ 
  // ----------- measure the 2D two-point correlation function, xi(r,mu) ----------- 
  
  TwoPointCorrelation2D_polar::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  
  // ----------- measure the first three multipoles of the two-point correlation function ----------- 

  vector<double> rad, xil, error;

  rad.resize(m_dd->nbins_D1()*3);
  xil.resize(0); xil.resize(m_dd->nbins_D1()*3, 0.);
  error.resize(0); error.resize(m_dd->nbins_D1()*3, 0.);

  double binSize = 1./m_dd->binSize_inv_D2();
  vector<vector<double> > xi2d = TwoPointCorrelation2D_polar::m_dataset->fxy();
  vector<vector<double> > error_xi2d = TwoPointCorrelation2D_polar::m_dataset->error_fxy();

  for (int i=0; i<m_dd->nbins_D1(); i++) {

    rad[i] = m_dd->scale_D1(i);
    rad[i+m_dd->nbins_D1()] = m_dd->scale_D1(i);
    rad[i+2*m_dd->nbins_D1()] = m_dd->scale_D1(i);
    
    for (int j=0; j<m_dd->nbins_D2(); j++) {
      double mu = m_dd->scale_D2(j);
      
      xil[i]                    += xi2d[i][j]*binSize;            // xi_0
      xil[i+m_dd->nbins_D1()]   += 5.*xi2d[i][j]*P_2(mu)*binSize; // xi_2
      xil[i+2*m_dd->nbins_D1()] += 9.*xi2d[i][j]*P_4(mu)*binSize; // xi_4

      error[i]                    += pow(error_xi2d[i][j]*binSize, 2);            // error[xi_0]
      error[i+m_dd->nbins_D1()]   += pow(5.*error_xi2d[i][j]*P_2(mu)*binSize, 2); // error[xi_2]
      error[i+2*m_dd->nbins_D1()] += pow(9.*error_xi2d[i][j]*P_4(mu)*binSize, 2); // error[xi_4]
    }

  }
  
  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv); } );

  m_dataset->set_xx(rad);
  m_dataset->set_fx(xil);
  m_dataset->set_error_fx(error);
}

*/
// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_multipoles::write (const string dir, const string file, const int rank) const 
{
  vector<double> rad = m_dataset->xx();
  vector<double> xil = m_dataset->fx();
  vector<double> error = m_dataset->error_fx();

  checkDim(rad, m_dd->nbins_D1()*3, "rad");
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);

  fout << "### rad  xi_0  error[xi_0]  xi_2  error[xi_2]  xi_4  error[xi_4] ###" << endl;

  for (int i=0; i<m_dd->nbins_D1(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << rad[i] << "  " << setw(8) << xil[i] << "  " << setw(8) << error[i] << "  " << setw(8) << xil[i+m_dd->nbins_D1()] << "  " << setw(8) << error[i+m_dd->nbins_D1()] << "  " << setw(8) << xil[i+2*m_dd->nbins_D1()] << "  " << setw(8) << error[i+2*m_dd->nbins_D1()] << endl;
   
  fout.close(); cout << endl << "I wrote the file: " << file_out << endl << endl;
}
