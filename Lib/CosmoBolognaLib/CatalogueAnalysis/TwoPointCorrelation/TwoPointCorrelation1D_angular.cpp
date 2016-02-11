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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_angular.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D_angular used to
 *  measure the angular two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D_angular used to measure the angular two-point
 *  correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D_angular.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::set_parameters (const binType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift) 
{
  m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, thetaMin, thetaMax, nbins, shift))
    : move(Pair::Create(_comoving_lin_, thetaMin, thetaMax, nbins, shift));
  
  m_rr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, thetaMin, thetaMax, nbins, shift))
    : move(Pair::Create(_comoving_lin_, thetaMin, thetaMax, nbins, shift));
  
  m_dr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, thetaMin, thetaMax, nbins, shift))
    : move(Pair::Create(_comoving_lin_, thetaMin, thetaMax, nbins, shift));
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::set_parameters (const binType binType, const double thetaMin, const double thetaMax, const double binSize, const double shift)
{
  m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, thetaMin, thetaMax, binSize, shift))
    : move(Pair::Create(_comoving_lin_, thetaMin, thetaMax, binSize, shift));
  
  m_rr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, thetaMin, thetaMax, binSize, shift))
    : move(Pair::Create(_comoving_lin_, thetaMin, thetaMax, binSize, shift));
  
  m_dr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, thetaMin, thetaMax, binSize, shift))
    : move(Pair::Create(_comoving_lin_, thetaMin, thetaMax, binSize, shift));
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::read(const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::write (const string dir, const string file, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins(), "theta");
  m_dataset->write(dir, file, "theta", "w", rank);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs,  const string dir_output_ResampleXi, int nMocks, int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  switch(errorType){
    case(ErrorType::_Poisson_):
      measurePoisson(dir_output_pairs,dir_input_pairs,count_dd,count_rr,count_dr,tcount);
      break;
    case(ErrorType::_Jackknife_):
      measureJackknife(dir_output_pairs,dir_input_pairs,dir_output_ResampleXi,count_dd,count_rr,count_dr,tcount);
      break;
    case(ErrorType::_Bootstrap_):
      measureBootstrap(nMocks,dir_output_pairs,dir_input_pairs,dir_output_ResampleXi,count_dd,count_rr,count_dr,tcount);
      break;
    default:
      ErrorMsg("Error in measure() of TwoPointCorrelation1D_angular.cpp, unknown type of error");
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  // ----------- weigthed number of objects in the real and random catalogues ----------- 
  
  int nData = m_data->weightedN();
  int nRandom = m_random->weightedN();
  
  if (nData==0 || nRandom==0)  
    ErrorMsg("Error in measurePoisson() of TwoPointCorrelation.cpp!");

  
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  cout << m_twoPType << endl;
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
  
  // ----------- compute the angular of the two-point correlation function ----------- 

  if(count_dr>-1)
    m_dataset = LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, nData, nRandom);
  else
    m_dataset = NaturalEstimatorTwoP(m_dd, m_rr, nData, nRandom);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_JackknifeXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{

  if(dir_output_JackknifeXi!="NULL"){
    string mkdir = "mkdir -p "+dir_output_JackknifeXi;
    if(system(mkdir.c_str())){}
  }
  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector<vector<double> > xi_SubSamples,covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  vector<shared_ptr<Data> > data_SS = (count_dr) ? XiJackknife(dd_regions,rr_regions,dr_regions) : XiJackknife(dd_regions,rr_regions);

  for (int i=0; i<nRegions; i++) {

    if(dir_output_JackknifeXi !="NULL"){
      string file = "xi_Jackknife_"+conv(i, par::fINT);
      data_SS[i]->write(dir_output_JackknifeXi, file, "theta", "w", 0);
    }

    xi_SubSamples.push_back(data_SS[i]->fx());
  }

  covariance_matrix(xi_SubSamples, covariance, 1);

  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();

  if(count_dr>-1)
    m_dataset = LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, nData, nRandom);
  else
    m_dataset = NaturalEstimatorTwoP(m_dd, m_rr, nData, nRandom);

  m_dataset->set_covariance_fx(covariance);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_angular::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_BootstrapXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  if (nMocks <=0)
    ErrorMsg("Error in measureBootstrap() of TwoPointCorrelation1D_angular.cpp, number of mocks must be >0");

  if(dir_output_BootstrapXi!="NULL"){
    string mkdir = "mkdir -p "+dir_output_BootstrapXi;
    if(system(mkdir.c_str())){}
  }

  vector<vector<double> > xi_SubSamples,covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
  vector<shared_ptr<Data> > data_SS = (count_dr) ? XiBootstrap(nMocks,dd_regions,rr_regions,dr_regions) : XiBootstrap(nMocks,dd_regions,rr_regions);

  for (int i=0; i<nMocks; i++) {

     if (dir_output_BootstrapXi!="NULL") {
      string file = "xi_Bootstrap_"+conv(i, par::fINT);
      data_SS[i]->write(dir_output_BootstrapXi, file, "theta", "w", 0);
    }

    xi_SubSamples.push_back(data_SS[i]->fx());
  }

  covariance_matrix(xi_SubSamples, covariance, 0);

  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();

  if (count_dr>-1)
    m_dataset = LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, nData, nRandom);
  else
    m_dataset = NaturalEstimatorTwoP(m_dd, m_rr, nData, nRandom);

  m_dataset->set_covariance_fx(covariance);

}
