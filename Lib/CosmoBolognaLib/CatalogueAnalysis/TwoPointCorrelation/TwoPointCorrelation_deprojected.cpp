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
#include "Data1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation_deprojected::DeProjectedTwoP (const vector<double> rp, const vector<double> ww, const vector<double> error_ww)
{

  vector<double> rad(rp.size(),0),xi(rp.size(),0),error_xi(rp.size(),0);

  for (size_t i=0; i<rp.size(); i++) {
    
    double ri = rp[i];
    rad[i] = ri;

    for (size_t j=i; j<rp.size()-1; j++) {
      double rj = rp[j];
      double rj1 = rp[j+1];
      double fact = 1./(rj1-rj)*log((rj1+sqrt(max(0., rj1*rj1-ri*ri)))/(rj+sqrt(max(0., rj*rj-ri*ri))));

      xi[i] -= (ww[j+1]-ww[j])*fact;
      error_xi[i] += pow(error_ww[j+1]*fact, 2)+pow(error_ww[j]*fact, 2);
    }
    
  }

  for_each( xi.begin(), xi.end(), [] (double &vv) { vv /= par::pi;} );
  for_each( error_xi.begin(), error_xi.end(), [] (double &vv) { vv = sqrt(vv/par::pi);} );

  return move(unique_ptr<Data1D>(new Data1D(rad,xi,error_xi)));

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measure(const double piMax_integral, const ErrorType errType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int nMocks, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  switch(errType){
    case(ErrorType::_Poisson_):
      measurePoisson(piMax_integral,dir_output_pairs,dir_input_pairs,count_dd,count_rr,count_dr,tcount);
      break;
    case(ErrorType::_Jackknife_):
      measureJackknife(piMax_integral,dir_output_pairs,dir_input_pairs,dir_output_ResampleXi,count_dd,count_rr,count_dr,tcount);
      break;
    case(ErrorType::_Bootstrap_):
      measureBootstrap(piMax_integral, nMocks, dir_output_pairs,dir_input_pairs,dir_output_ResampleXi,count_dd,count_rr,count_dr,tcount);
      break;
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measurePoisson (const double piMax_integral, const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{

  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation_projected::measurePoisson(piMax_integral,dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 
  
  m_dataset = DeProjectedTwoP(TwoPointCorrelation_projected::dataset()->xx(),TwoPointCorrelation_projected::dataset()->fx(),TwoPointCorrelation_projected::dataset()->error_fx());

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measureJackknife (const double piMax_integral, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{

  if (dir_output_ResampleXi != par::defaultString){
    string mkdir = "mkdir -p "+dir_output_ResampleXi;
    if(system(mkdir.c_str())){}
  }

  vector<shared_ptr<Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount);

  shared_ptr<Pair> dd_cart = move(Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2()));

  shared_ptr<Pair> rr_cart = move(Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2()));

  shared_ptr<Pair> dr_cart = move(Pair::Create(m_dr->pairType(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2()));

  for (size_t k=0; k<dd_regions.size(); k++)
    for (int i=0; i<dd_regions[k]->nbins_D1(); i++)
      for (int j=0; j<dd_regions[k]->nbins_D2(); j++) {
	dd_cart->add_PP2D(i, j, dd_regions[k]->PP2D(i, j));
	rr_cart->add_PP2D(i, j, rr_regions[k]->PP2D(i, j));
	dr_cart->add_PP2D(i, j, dr_regions[k]->PP2D(i, j));
      }
  
  auto data_cart = (count_dr>=0) ? LandySzalayEstimatorTwoP(dd_cart, rr_cart, dr_cart, m_data->weightedN(), m_random->weightedN()) : NaturalEstimatorTwoP(dd_cart, rr_cart, m_data->weightedN(), m_random->weightedN());

  auto data_proj = TwoPointCorrelation_projected::ProjectedTwoP(piMax_integral, data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());

  if (count_dr==1) 
    data = XiJackknife(piMax_integral, dd_regions, rr_regions,dr_regions);
  else
    data = XiJackknife(piMax_integral, dd_regions, rr_regions);

  vector<vector<double> > ww,covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());

    if (dir_output_ResampleXi != par::defaultString) {
      string filename = "xi_deprojected_Jackkknife_"+conv(i,par::fINT);
      data[i]->write(dir_output_ResampleXi, filename, "rp", "xi_deprojected", 0);
    }
  }

  covariance_matrix(ww, covariance, 1);

  m_dataset = DeProjectedTwoP(data_proj->xx(), data_proj->fx(), data_proj->error_fx());
  m_dataset->set_covariance_fx(covariance);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measureBootstrap (const double piMax_integral, const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  if (dir_output_ResampleXi != par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_ResampleXi;
    if (system(mkdir.c_str())) {}
  }
  
  vector<shared_ptr<Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr,  tcount);

  shared_ptr<Pair> dd_cart = move(Pair::Create(m_dd->pairType(), m_dd->sMin_D1(), m_dd->sMax_D1(), m_dd->nbins_D1(), m_dd->shift_D1(), m_dd->sMin_D2(), m_dd->sMax_D2(), m_dd->nbins_D2(), m_dd->shift_D2()));

  shared_ptr<Pair> rr_cart = move(Pair::Create(m_rr->pairType(), m_rr->sMin_D1(), m_rr->sMax_D1(), m_rr->nbins_D1(), m_rr->shift_D1(), m_rr->sMin_D2(), m_rr->sMax_D2(), m_rr->nbins_D2(), m_rr->shift_D2()));

  shared_ptr<Pair> dr_cart = move(Pair::Create(m_dr->pairType(), m_dr->sMin_D1(), m_dr->sMax_D1(), m_dr->nbins_D1(), m_dr->shift_D1(), m_dr->sMin_D2(), m_dr->sMax_D2(), m_dr->nbins_D2(), m_dr->shift_D2()));
  
  for (size_t k=0; k<dd_regions.size(); k++)
    for (int i=0; i<dd_regions[k]->nbins_D1(); i++)
      for (int j=0; j<dd_regions[k]->nbins_D2(); j++) {
	dd_cart->add_PP2D(i, j, dd_regions[k]->PP2D(i, j));
	rr_cart->add_PP2D(i, j, rr_regions[k]->PP2D(i, j));
	dr_cart->add_PP2D(i, j, dr_regions[k]->PP2D(i, j));
      }

  auto data_cart = (count_dr>=0) ? LandySzalayEstimatorTwoP(dd_cart, rr_cart, dr_cart, m_data->weightedN(), m_random->weightedN()) : NaturalEstimatorTwoP(dd_cart, rr_cart, m_data->weightedN(), m_random->weightedN());

  auto data_proj = TwoPointCorrelation_projected::ProjectedTwoP(piMax_integral, data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());

  if (count_dr==1) 
    data = XiBootstrap(piMax_integral, nMocks, dd_regions, rr_regions, dr_regions);
  else
    data = XiBootstrap(piMax_integral, nMocks, dd_regions, rr_regions);

  vector<vector<double> > ww,covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());
    if (dir_output_ResampleXi != par::defaultString) {
      string filename = "xi_deprojected_Jackkknife_"+conv(i,par::fINT);
      data[i]->write(dir_output_ResampleXi, filename, "rp", "xi_deprojected", 0);
    }
  }
  covariance_matrix(ww, covariance, 0);

  m_dataset = DeProjectedTwoP(data_proj->xx(), data_proj->fx(), data_proj->error_fx());
  m_dataset->set_covariance_fx(covariance);
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_deprojected::XiJackknife(const double piMax_integral, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<Data> > data;
  
  auto data_proj = TwoPointCorrelation_projected::XiJackknife(piMax_integral, dd, rr);

  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(DeProjectedTwoP(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_deprojected::XiJackknife(const double piMax_integral, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<Data> > data;
  
  auto data_proj = TwoPointCorrelation_projected::XiJackknife(piMax_integral, dd, rr, dr);
  
  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(DeProjectedTwoP(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_deprojected::XiBootstrap(const double piMax_integral, const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<Data> > data;

  auto data_proj = TwoPointCorrelation_projected::XiBootstrap(piMax_integral, nMocks, dd, rr);

  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(DeProjectedTwoP(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_deprojected::XiBootstrap(const double piMax_integral, const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<Data> > data;
  
  auto data_proj = TwoPointCorrelation_projected::XiBootstrap(piMax_integral,nMocks,dd,rr,dr);
  
  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(DeProjectedTwoP(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
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
