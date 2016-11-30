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
 *  @file
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_filtered.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D_filtered used to
 *  measure the filtered monopole of the two-point correlation
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D_filtered used to measure the filtered
 *  monopole of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D_filtered.h"
#include "Data1D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters (const binType binType, const double rMin, const double rMax, const int nbins, const double shift) 
{
  if (binType==_linear_){
    double binSize = ((rMax-rMin)/nbins);

    m_rc.resize(nbins);
    for (int i=0; i<nbins; i++)
      m_rc[i] = binSize*(i+shift)+rMin;
  }
  else if (binType==_logarithmic_){
    if (rMin<1.e-30) ErrorCBL("Error in cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters of TwoPointCorrelation1D_filtered.cpp: Min must be >0!");

    double binSize = ((log10(rMax)-log10(rMin))/nbins);

    m_rc.resize(nbins);
    for (int i=0; i<nbins; i++)
      m_rc[i] = pow(10.,(i+shift)*binSize+log10(rMin));
  }
  else ErrorCBL("Error in cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters of TwoPointCorrelation1D_filtered.cpp: no such type of binning!");

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters (const binType binType, const double rMin, const double rMax, const double binSize, const double shift)
{
  if (binType==_linear_) {
    int nbins = nint((rMax-rMin)/binSize);

    m_rc.resize(nbins);
    for (int i=0; i<nbins; i++)
      m_rc[i] = (i+shift)*binSize+rMin;
  }
  
  else if (binType==_logarithmic_) {
    if (rMin<1.e-30) ErrorCBL("Error in cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters of TwoPointCorrelation1D_filtered.cpp: Min must be >0!");

    int nbins = nint((log10(rMax)-log10(rMin))/binSize);

    m_rc.resize(nbins);
    for (int i=0; i<nbins; i++)
      m_rc[i] = pow(10.,(i+shift)*binSize+log10(rMin));

  }

  else ErrorCBL("Error in cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters of TwoPointCorrelation1D_filtered.cpp: no such type of binning!");
}


// ============================================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation1D_filtered::Filtered (const shared_ptr<data::Data> monopole)
{
  vector<double> rad = (!m_compute_extra_info) ? monopole->xx() : monopole->extra_info()[0];
  vector<double> xi = monopole->fx();
  vector<double> error_xi = monopole->error_fx();
  double binSize = rad[1]-rad[0];

  vector<double> wc(m_rc.size(), 0);
  vector<double> error_wc(m_rc.size(), 0);

  for (size_t i=0; i<m_rc.size(); i++) 
    for (size_t j=0; j<rad.size(); j++) {
      if (rad[j]<m_rc[i]) {
	wc[i] += rad[j]*rad[j]*binSize*xi[j]*Filter(rad[j], m_rc[i]);
	error_wc[i] += pow(rad[j]*rad[j]*binSize*error_xi[j]*Filter(rad[j], m_rc[i]),2);
      }
    }

  for_each( error_wc.begin(), error_wc.end(), [] (double &vv) { vv = sqrt(vv);} );

  //return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(m_rc, wc, error_wc))) : data_with_extra_info(m_rc, wc, error_wc);
  return move(unique_ptr<data::Data1D>(new data::Data1D(m_rc, wc, error_wc)));
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  switch (errorType) {
  case (ErrorType::_Poisson_) :
    measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case (ErrorType::_Jackknife_) :
    measureJackknife(dir_output_pairs, dir_input_pairs, dir_output_ResampleXi, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case (ErrorType::_Bootstrap_) :
    measureBootstrap(nMocks, dir_output_pairs, dir_input_pairs, dir_output_ResampleXi, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  default:
    ErrorCBL("Error in measure() of TwoPointCorrelation1D_filtered.cpp, unknown type of error");
  }
}

// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator) 
{
  // ----------- weigthed number of objects in the real and random catalogues ----------- 
  
  int nData = m_data->weightedN();
  int nRandom = m_random->weightedN();
  
  if (nData==0 || nRandom==0)  
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation.cpp!");

  
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
  
  
  // ----------- compute the monopole of the two-point correlation function ----------- 

  if (estimator==_natural_)
    m_dataset = Filtered(NaturalEstimator(m_dd, m_rr, nData, nRandom));
  else if (estimator==_LandySzalay_)
    m_dataset = Filtered(LandySzalayEstimator(m_dd, m_rr, m_dr, nData, nRandom));
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation1D_filtered.cpp: the chosen estimator is not implemented!");
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();
  
  vector<vector<double> > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
  
  vector<shared_ptr<Data> > data_SS = (estimator==_natural_) ? XiJackknife(dd_regions, rr_regions) : XiJackknife(dd_regions, rr_regions, dr_regions);
  
  for (size_t i=0; i<nRegions; i++) {
    shared_ptr<Data> data_filtered = Filtered(data_SS[i]);
    
    if (dir_output_resample != par::defaultString && dir_output_resample != "") {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      data_filtered->write(dir_output_resample, file, "[1] separation at the bin centre # [2] filtered two-point correlation function # [3] error", 0);
    }

    xi_SubSamples.push_back(data_filtered->fx());
  }

  covariance_matrix(xi_SubSamples, covariance, true);

  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();
  
  if (estimator==_natural_)
    m_dataset = Filtered(NaturalEstimator(m_dd, m_rr, nData, nRandom));
  else if (estimator==_LandySzalay_)
    m_dataset = Filtered(LandySzalayEstimator(m_dd, m_rr, m_dr, nData, nRandom));
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation1D_filtered.cpp: the chosen estimator is not implemented!");
  
  m_dataset->set_covariance(covariance);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (nMocks<=0)
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation1D_filtered.cpp, number of mocks must be >0");

  if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }

  vector<vector<double> > xi_SubSamples,covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
  vector<shared_ptr<Data> > data_SS = (estimator==_natural_) ? XiBootstrap(nMocks, dd_regions, rr_regions) : XiBootstrap(nMocks,dd_regions,rr_regions,dr_regions);

  for (int i=0; i<nMocks; i++) {
    shared_ptr<Data> data_filtered = Filtered(data_SS[i]);

     if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
      string file = "xi_Bootstrap_"+conv(i, par::fINT)+".dat";
      data_filtered->write(dir_output_resample, file, "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error", 0);
    }

    xi_SubSamples.push_back(data_filtered->fx());
  }

  covariance_matrix(xi_SubSamples, covariance, 0);

  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();
  
  if (estimator==_natural_)
    m_dataset = Filtered(NaturalEstimator(m_dd, m_rr, nData, nRandom));
  else if (estimator==_LandySzalay_)
    m_dataset = Filtered(LandySzalayEstimator(m_dd, m_rr, m_dr, nData, nRandom));
  ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation1D_filtered.cpp: the chosen estimator is not implemented!");
  
  m_dataset->set_covariance(covariance);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::write (const string dir, const string file, const int rank) const 
{
  string header = "[1] separation at the bin centre # [2] filtered two-point correlation function # [3] error";
  if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution";
  
  m_dataset->write(dir, file, header, rank);
}

