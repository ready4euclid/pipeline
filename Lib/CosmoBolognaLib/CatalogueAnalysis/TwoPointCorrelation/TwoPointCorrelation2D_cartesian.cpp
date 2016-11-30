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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D_cartesian.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation2D_cartesian used
 *  to measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation2D_cartesian used to measure the monopole of
 *  the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation2D_cartesian.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::set_parameters (const binType binType_rp, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const binType binType_pi, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info) 
{
  if (!compute_extra_info) {
    if (binType_rp==_logarithmic_) {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_loglog_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglog_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglog_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_loglin_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglin_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglin_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
    }
    else {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_linlog_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlog_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlog_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_linlin_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlin_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlin_, _standard_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
    }
  }
  else {
    if (binType_rp==_logarithmic_) {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_loglog_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglog_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglog_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_loglin_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglin_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglin_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
    }
    else {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_linlog_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlog_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlog_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_linlin_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlin_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlin_, _extra_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, angularUnits));
      }
    }
  }    
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::set_parameters (const binType binType_rp, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const binType binType_pi, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info)
{
  if (!compute_extra_info) {
    if (binType_rp==_logarithmic_) {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_loglog_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglog_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglog_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_loglin_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglin_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglin_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
    }
    else {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_linlog_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlog_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlog_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_linlin_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlin_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlin_, _standard_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
    }
  }
  else {
    if (binType_rp==_logarithmic_) {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_loglog_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglog_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglog_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_loglin_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_loglin_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_loglin_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
    }
    else {
      if (binType_pi==_logarithmic_) {
	m_dd = move(Pair::Create(_comovingCartesian_linlog_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlog_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlog_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
      else {
	m_dd = move(Pair::Create(_comovingCartesian_linlin_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight));
	m_rr = move(Pair::Create(_comovingCartesian_linlin_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
	m_dr = move(Pair::Create(_comovingCartesian_linlin_, _extra_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, angularUnits));
      }
    }
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::write (const string dir, const string file, const bool full, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins_D1(), "rp"); checkDim(m_dataset->yy(), m_dd->nbins_D2(), "pi");
  
  string header = "[1] perpendicular separation at the bin centre # [2] parallel separation at the bin centre # [3] 2D two-point correlation function # [4] error";
  if (m_compute_extra_info) header += " # [5] mean perpendicular separation # [6] standard deviation of the distribution of perpendicular separations # [7] mean parallel separation # [8] standard deviation of the distribution of parallel separations # [9] mean redshift # [10] standard deviation of the redshift distribution";
  
  m_dataset->write(dir, file, header, full, rank);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  switch (errorType) {
  case (ErrorType::_Poisson_) :
    measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case (ErrorType::_Jackknife_) :
    measureJackknife(dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case (ErrorType::_Bootstrap_) :
    measureBootstrap(nMocks, dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  default:
    ErrorCBL("Error in measure() of TwoPointCorrelation2D_cartesian.cpp, unknown type of error");
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  // ----------- weigthed number of objects in the real and random catalogues ----------- 
  
  int nData = m_data->weightedN();
  int nRandom = m_random->weightedN();
  
  if (nData==0 || nRandom==0)  
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation2D_cartesian.cpp!");

  
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
  
  
  // ----------- compute the monopole of the two-point correlation function ----------- 

  if (estimator==_natural_)
    m_dataset = NaturalEstimator(m_dd, m_rr, nData, nRandom);
  else if (estimator==_LandySzalay_)
    m_dataset = LandySzalayEstimator(m_dd, m_rr, m_dr, nData, nRandom);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation2D_cartesian.cpp: the chosen estimator is not implemented!");
  
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_JackknifeXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_JackknifeXi!=par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_JackknifeXi;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector< vector< vector<double > > > xi_SubSample;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  vector<shared_ptr<Data> > data_SS = (estimator==_natural_) ? XiJackknife(dd_regions, rr_regions) : XiJackknife(dd_regions, rr_regions, dr_regions);

  for (size_t i=0; i<nRegions; i++) {

    if (dir_output_JackknifeXi !=par::defaultString) {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      data_SS[i]->write(dir_output_JackknifeXi, file, "[1] perpendicular separation at the bin centre # [2] parallel separation at the bin centre # [3] 2D two-point correlation function # [4] error", 0);
    }

    xi_SubSample.push_back(data_SS[i]->fxy());
  }

  vector<vector<double> > error(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(), 0));

  double fact = pow(nRegions-1, 2)/nRegions;
  for (int i=0; i<m_dd->nbins_D1(); i++) {
    for (int j=0; j<m_dd->nbins_D2(); j++) {
      vector<double> temp;
      for (size_t nm=0; nm<nRegions; nm++) 
	temp.push_back(xi_SubSample[nm][i][j]);
      error[i][j] = sqrt(pow(Sigma(temp), 2)*fact);
    }
  }

  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();

  if (estimator==_natural_)
    m_dataset = NaturalEstimator(m_dd, m_rr, nData, nRandom);
  else if (estimator==_LandySzalay_)
    m_dataset = LandySzalayEstimator(m_dd, m_rr, m_dr, nData, nRandom);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation2D_cartesian.cpp: the chosen estimator is not implemented!");

  m_dataset->set_error_fxy(error);
}

// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_cartesian::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_BootstrapXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (nMocks<1)
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation2D_cartesian.cpp, number of mocks must be >0");

  if (dir_output_BootstrapXi!=par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_BootstrapXi;
    if (system(mkdir.c_str())) {}
  }

  vector< vector< vector<double > > > xi_SubSample;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  vector<shared_ptr<Data> > data_SS = (estimator==_natural_) ? XiBootstrap(nMocks, dd_regions, rr_regions) : XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions);

  for (int i=0; i<nMocks; i++) {
    if (dir_output_BootstrapXi!=par::defaultString) {
      string file = "xi_Bootstrap_"+conv(i, par::fINT)+".dat";
      data_SS[i]->write(dir_output_BootstrapXi, file, "[1] perpendicular separation at the bin centre # [2] parallel separation at the bin centre # [3] 2D two-point correlation function # [4] error", 0);
    }
    xi_SubSample.push_back(data_SS[i]->fxy());
  }

  vector<vector<double> > error(m_dd->nbins_D1(), vector<double>(m_dd->nbins_D2(),0));

  for (int i=0; i<m_dd->nbins_D1(); i++) {
    for (int j=0; j<m_dd->nbins_D2(); j++) {
      vector<double> temp;
      for (int nm = 0; nm<nMocks;nm++) 
	temp.push_back(xi_SubSample[nm][i][j]);
      error[i][j] = Sigma(temp);
    }
  }

  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();

  if (estimator==_natural_)
    m_dataset = NaturalEstimator(m_dd, m_rr, nData, nRandom);
  else if (estimator==_LandySzalay_)
    m_dataset = LandySzalayEstimator(m_dd, m_rr, m_dr, nData, nRandom);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation2D_cartesian.cpp: the chosen estimator is not implemented!");

  m_dataset->set_error_fxy(error);
}
