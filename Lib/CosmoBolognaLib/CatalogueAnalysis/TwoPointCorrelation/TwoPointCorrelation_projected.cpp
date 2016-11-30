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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_projected.cpp
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
#include "Data1D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace pairs;
using namespace twopt;


// ============================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation_projected::data_with_extra_info (const vector<double> rp, const vector<double> ww, const vector<double> error) const
{
  auto dd2D = cosmobl::twopt::TwoPointCorrelation2D_cartesian::m_dd;
  
  vector<double> scale_mean(dd2D->nbins_D1(), 0.), scale_sigma(dd2D->nbins_D1(), 0.), z_mean(dd2D->nbins_D1(), 0.), z_sigma(dd2D->nbins_D1(), 0.);

  for (int i=0; i<dd2D->nbins_D1(); ++i) 
    for (int j=0; j<dd2D->nbins_D2(); ++j) {
      scale_mean[i] += dd2D->scale_D1_mean(i, j);
      scale_sigma[i] += dd2D->scale_D1_sigma(i, j);
      z_mean[i] += dd2D->z_mean(i, j);
      z_sigma[i] += dd2D->z_sigma(i, j);
    }
  
  vector<vector<double>> extra(4);
  
  for (int i=0; i<dd2D->nbins_D1(); ++i) {
    extra[0].push_back(scale_mean[i]/dd2D->nbins_D2());
    extra[1].push_back(scale_sigma[i]/dd2D->nbins_D2());
    extra[2].push_back(z_mean[i]/dd2D->nbins_D2());
    extra[3].push_back(z_sigma[i]/dd2D->nbins_D2());
  }
  
  return move(unique_ptr<data::Data1D_extra>(new data::Data1D_extra(rp, ww, error, extra)));
}


// ============================================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation_projected::Projected (const vector<double> rp, const vector<double> pi, const vector<vector<double>> xi, const vector<vector<double>> error_xi)
{
  vector<double> ww, error;

  ww.resize(0); ww.resize(rp.size(), 0.);
  error.resize(0); error.resize(rp.size(), 0.);
 
  double binSize = pi[1]-pi[0];

  int pim = nint((m_piMax_integral-Min(pi))/binSize); // to convert from Mpc/h into the vector index

  for (size_t i=0; i<rp.size(); i++) {

    ww[i] = 0.;
    error[i] = 0.;
    
    for (int j=0; j<pim; j++) {  
      if (xi[i][j]>-1.e29) { // check!!!
	ww[i] = ww[i]+2.*binSize*xi[i][j];
	error[i] += pow(2.*binSize*error_xi[i][j], 2); // check!!!!
      }
    }
  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv);} );

  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rp, ww, error))) : data_with_extra_info(rp, ww, error);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::write (const string dir, const string file, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins_D1(), "rp");

  string header =  "[1] perpendicular separation at the bin centre # [2] projected two-point correlation function # [3] error";
  if (m_compute_extra_info) header += " # [4] mean perpendicular separation # [5] standard deviation of the distribution of perpendicular separations";

  m_dataset->write(dir, file, header, rank);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
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
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation2D_cartesian::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  
  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 
  
  m_dataset = Projected(TwoPointCorrelation2D_cartesian::xx(), TwoPointCorrelation2D_cartesian::yy(), TwoPointCorrelation2D_cartesian::xi2D(), TwoPointCorrelation2D_cartesian::error2D());
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample != par::defaultString && dir_output_resample != "") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }

  vector<shared_ptr<data::Data>> data;
  vector<shared_ptr<pairs::Pair>> dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  auto data_cart = (estimator==_natural_) ? NaturalEstimator(m_dd, m_rr, m_data->weightedN(), m_random->weightedN()) : LandySzalayEstimator(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN());
  
  if (estimator==_natural_) 
    data = XiJackknife(dd_regions, rr_regions);
  else if (estimator==_LandySzalay_)
    data = XiJackknife(dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation_projected.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double>> ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());
    if (dir_output_resample != par::defaultString && dir_output_resample != "") {
      string file = "xi_projected_Jackkknife_"+conv(i, par::fINT)+".dat";
      data[i]->write(dir_output_resample, file, "[1] perpendicular separation at the bin centre # [2] projected two-point correlation function # [3] error", 0);
    }
  }
  
  covariance_matrix(ww, covariance, true);

  m_dataset = Projected(data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());
  m_dataset->set_covariance(covariance);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample != par::defaultString && dir_output_resample != "") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<shared_ptr<data::Data>> data;
  vector<shared_ptr<pairs::Pair>> dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  auto data_cart = (estimator==_natural_) ? NaturalEstimator(m_dd, m_rr, m_data->weightedN(), m_random->weightedN()) : LandySzalayEstimator(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN());

  if (estimator==_natural_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions);
  else if (estimator==_LandySzalay_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation_projected.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double>> ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());
    if (dir_output_resample != par::defaultString && dir_output_resample != "") {
      string file = "xi_projected_Bootstrap_"+conv(i, par::fINT)+".dat";
      data[i]->write(dir_output_resample, file, "[1] perpendicular separation at the bin centre # [2] projected two-point correlation function # [3] error", 0);
    }
  }
  
  covariance_matrix(ww, covariance, false);

  m_dataset = Projected(data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());
  m_dataset->set_covariance(covariance);
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_projected::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<shared_ptr<data::Data>> data;
  
  auto data2d = TwoPointCorrelation2D_cartesian::XiJackknife(dd, rr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(Projected(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));

  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_projected::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<shared_ptr<data::Data>> data;
 
  auto data2d = TwoPointCorrelation2D_cartesian::XiJackknife(dd, rr, dr);
  
  for (size_t i=0; i<data2d.size(); i++) 
    data.push_back(move(Projected(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));

  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_projected::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<shared_ptr<data::Data>> data;

  auto data2d = TwoPointCorrelation2D_cartesian::XiBootstrap(nMocks, dd, rr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(Projected(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));

  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_projected::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<shared_ptr<data::Data>> data;

  auto data2d = TwoPointCorrelation2D_cartesian::XiBootstrap(nMocks, dd, rr, dr);
  
  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(Projected(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  
  return data;
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::read_covariance_matrix (const string dir, const string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::write_covariance_matrix (const string dir, const string file) const
{
  m_dataset->write_covariance(dir, file, "r");
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::compute_covariance_matrix (const vector<shared_ptr<data::Data>> xi_collection, const bool doJK)
{
  vector<vector<double>> xi;

  for(size_t i=0;i<xi_collection.size();i++)
    xi.push_back(xi_collection[i]->fx());

  vector<vector<double>> cov_mat;
  cosmobl::covariance_matrix(xi, cov_mat, doJK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::compute_covariance_matrix (const vector<string> file_xi, const bool doJK)
{
  vector<double> rad, mean;
  vector<vector<double>> cov_mat;

  cosmobl::covariance_matrix (file_xi, rad, mean, cov_mat, doJK); 
  m_dataset->set_covariance(cov_mat);
}
