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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D_polar.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation2D_polar used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation2D_polar used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation2D_polar.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::set_parameters (const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu) 
{
  if (binType_rad==_logarithmic_) {
    if (binType_mu==_logarithmic_) {
      m_dd = move(Pair::Create(_comovingPolar_loglog_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_loglog_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_loglog_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
    }
    else {
      m_dd = move(Pair::Create(_comovingPolar_loglin_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_loglin_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_loglin_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
    }
  }
  else {
    if (binType_mu==_logarithmic_) {
      m_dd = move(Pair::Create(_comovingPolar_linlog_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_linlog_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_linlog_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
    }
    else {
      m_dd = move(Pair::Create(_comovingPolar_linlin_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_linlin_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_linlin_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu));
    }
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::set_parameters (const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu)
{
  if (binType_rad==_logarithmic_) {
    if (binType_mu==_logarithmic_) {
      m_dd = move(Pair::Create(_comovingPolar_loglog_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_loglog_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_loglog_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
    }
    else {
      m_dd = move(Pair::Create(_comovingPolar_loglin_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_loglin_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_loglin_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
    }
  }
  else {
    if (binType_mu==_logarithmic_) {
      m_dd = move(Pair::Create(_comovingPolar_linlog_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_linlog_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_linlog_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
    }
    else {
      m_dd = move(Pair::Create(_comovingPolar_linlin_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_rr = move(Pair::Create(_comovingPolar_linlin_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
      m_dr = move(Pair::Create(_comovingPolar_linlin_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu));
    }
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::read(const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::write (const string dir, const string file, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins_D1(), "rad"); checkDim(m_dataset->yy(), m_dd->nbins_D2(), "mu");
  m_dataset->write(dir, file, "r", "mu", "xi", rank);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::measure (const string dir_output_pairs, const vector<string> dir_input_pairs, const ErrorType errorType, const string dir_output_ResampleXi, int nMocks, int count_dd, const int count_rr, const int count_dr, const bool tcount)
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
      ErrorMsg("Error in measure of TwoPointCorrelation2D.cpp, unknown type of error");
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  // ----------- weigthed number of objects in the real and random catalogues ----------- 
  
  int nData = m_data->weightedN();
  int nRandom = m_random->weightedN();
  
  if (nData==0 || nRandom==0)  
    ErrorMsg("Error in cosmobl::twopt::TwoPointCorrelation::measurePoisson of TwoPointCorrelation.cpp!");

  
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
  
  // ----------- compute the monopole of the two-point correlation function ----------- 

  if(count_dr>-1)
    m_dataset = LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, nData, nRandom);
  else
    m_dataset = NaturalEstimatorTwoP(m_dd, m_rr, nData, nRandom);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_JackknifeXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{

  if(dir_output_JackknifeXi!="NULL"){
    string mkdir = "mkdir -p "+dir_output_JackknifeXi;
    if(system(mkdir.c_str())){}
  }
  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector< vector< vector<double > > > xi_SubSample;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  for (int i=0; i<nRegions; i++) {
    double nData_SubSample = m_data->weightedN_condition(Var::_REGION_, region_list[i], region_list[i]+1, 1);
    double nRandom_SubSample = m_random->weightedN_condition(Var::_REGION_, region_list[i], region_list[i]+1, 1);

    auto dd_SubSample = Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2());
    auto rr_SubSample = Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2());
    auto dr_SubSample = Pair::Create(m_dr->pairType(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2());

    vector<int> w(nRegions, 1);
    w[i] = 0;

    for (int j=0; j<nRegions; j++) {
      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = w[j]*w[k];
	if (ww>0) {
	  for (int bin1=0;bin1<dd_SubSample->nbins_D1();bin1++) {
	    for (int bin2=0;bin1<dd_SubSample->nbins_D2();bin2++) {
	      dd_SubSample->add_PP2D(bin1 ,bin2, dd_regions[index]->PP2D(bin1,bin2));
	      rr_SubSample->add_PP2D(bin1,bin2 , rr_regions[index]->PP2D(bin1,bin2));
	      if (count_dr>-1)
		dr_SubSample->add_PP2D(bin1, bin2, dr_regions[index]->PP2D(bin1,bin2));
	    }
	  }
	}
      }
    }
    
    shared_ptr<Data> data;
    if(count_dr>-1)
      data = LandySzalayEstimatorTwoP(dd_SubSample, rr_SubSample, dr_SubSample, nData_SubSample, nRandom_SubSample);
    else
      data = NaturalEstimatorTwoP(dd_SubSample, rr_SubSample, nData_SubSample, nRandom_SubSample);

    if(dir_output_JackknifeXi !="NULL"){
      string file = "xi_Jackknife_"+conv(i, par::fINT);
      data->write(dir_output_JackknifeXi, file, "rad", "xi", 0);
    }

    xi_SubSample.push_back(data->fxy());

  }

  vector<vector<double> > error;

  double fact = pow(nRegions-1,2)/nRegions;
  for(int i=0;i<m_dd->nbins_D1();i++){
    for(int j=0;j<m_dd->nbins_D2();j++){
      vector<double> temp;
      for(int nm = 0; nm<nRegions;nm++){
	temp.push_back(xi_SubSample[nm][i][j]);
      }
      error[i][j] = pow(Sigma(temp),2)*fact;
    }
  }


  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();

  if(count_dr>-1)
    m_dataset = LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, nData, nRandom);
  else
    m_dataset = NaturalEstimatorTwoP(m_dd, m_rr, nData, nRandom);

  m_dataset->set_error_fxy(error);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation2D_polar::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_BootstrapXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  if (nMocks <=0)
    ErrorMsg("Error in measureBootstrap() of TwoPointCorrelation2D_polar.cpp, number of mocks must be >0");

  if(dir_output_BootstrapXi!="NULL"){
    string mkdir = "mkdir -p "+dir_output_BootstrapXi;
    if(system(mkdir.c_str())){}
  }

  vector< vector< vector<double > > > xi_SubSample;

  vector<long> region_list = m_data->get_region_list();
  int nRegions = region_list.size();

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  vector<double> nData_Region, nRandom_Region;

  for (int i=0; i<nRegions; i++) {
    nData_Region.push_back(m_data->weightedN_condition(Var::_REGION_, region_list[i], region_list[i]+1, 0));
    nRandom_Region.push_back(m_random->weightedN_condition(Var::_REGION_, region_list[i], region_list[i]+1, 0));
    cout << nData_Region[i] << " " << nRandom_Region[i] << endl;

  }

  uniform_int_distribution<int> uni(0, nRegions-1);
  default_random_engine rng;
  int val=3; //See Norberg et al. 2009

  for (int i=0; i<nMocks; i++) {
    auto dd_SubSample = Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2());
    auto rr_SubSample = Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2());
    auto dr_SubSample = Pair::Create(m_dr->pairType(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2());

    double nData_SubSample=0, nRandom_SubSample=0;

    vector<int> w(nRegions, 0);
    for (int n=0; n<val*nRegions; n++)
      w[uni(rng)] +=1;

    for (int j=0; j<nRegions; j++) {
      nData_SubSample += w[j]*nData_Region[j];
      nRandom_SubSample += w[j]*nRandom_Region[j];

      for (int k=j; k<nRegions; k++) {
	int index = j*nRegions-(j-1)*j/2+k-j;
	double ww = (k==j) ? w[k] : w[j]*w[k];
	if (ww>0) {
	  for (int bin1=0;bin1<dd_SubSample->nbins_D1();bin1++) {
	    for (int bin2=0;bin1<dd_SubSample->nbins_D2();bin2++) {
	      dd_SubSample->add_PP2D(bin1 ,bin2, dd_regions[index]->PP2D(bin1,bin2));
	      rr_SubSample->add_PP2D(bin1, bin2 , rr_regions[index]->PP2D(bin1,bin2));
	      if (count_dr>-1)
		dr_SubSample->add_PP2D(bin1, bin2, dr_regions[index]->PP2D(bin1,bin2));
	    }
	  }
	}
      }
    }
    
    shared_ptr<Data> data;
    if (count_dr>-1)
      data = LandySzalayEstimatorTwoP(dd_SubSample, rr_SubSample, dr_SubSample, nData_SubSample, nRandom_SubSample);
    else
      data = NaturalEstimatorTwoP(dd_SubSample, rr_SubSample, nData_SubSample, nRandom_SubSample);

    if (dir_output_BootstrapXi!="NULL") {
      string file = "xi_Bootstrap_"+conv(i, par::fINT);
      data->write(dir_output_BootstrapXi, file, "rad", "xi", 0);
    }

    xi_SubSample.push_back(data->fxy());

  }

  vector<vector<double> > error;

  for(int i=0;i<m_dd->nbins_D1();i++){
    for(int j=0;j<m_dd->nbins_D2();j++){
      vector<double> temp;
      for(int nm = 0; nm<nRegions;nm++){
	temp.push_back(xi_SubSample[nm][i][j]);
      }
      error[i][j] = pow(Sigma(temp),2);
    }
  }

  double nData = m_data->weightedN();
  double nRandom = m_random->weightedN();

  if (count_dr>-1)
    m_dataset = LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, nData, nRandom);
  else
    m_dataset = NaturalEstimatorTwoP(m_dd, m_rr, nData, nRandom);

  m_dataset->set_error_fxy(error);
}
