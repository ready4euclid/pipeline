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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_wedges.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_multipoles used to
 *  measure the wedges of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_wedges used to measure the wedges of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_wedges.h"
#include "Data1D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace pairs;
using namespace twopt;


// ============================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation_wedges::data_with_extra_info (const vector<double> rad, const vector<double> wedges, const vector<double> error) const
{
  auto dd2D = cosmobl::twopt::TwoPointCorrelation2D_polar::m_dd;
  
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
  
  return move(unique_ptr<data::Data1D_extra>(new data::Data1D_extra(rad, wedges, error, extra)));
}


// ============================================================================================


vector<double> cosmobl::twopt::TwoPointCorrelation_wedges::xx () const
{
  vector<double> rad;

  for (size_t i=0; i<m_dataset->xx().size()/2; i++)
    rad.push_back(m_dataset->xx()[i]);

  return rad;
}


// ============================================================================================


vector<double> cosmobl::twopt::TwoPointCorrelation_wedges::xiPerpendicular () const
{
  size_t sz = m_dataset->fx().size();
  vector<double> xiPerp;
  
  for (size_t i=0; i<sz/2; i++)
    xiPerp.push_back(m_dataset->fx()[i]);

  return xiPerp;
}


// ============================================================================================


vector<double> cosmobl::twopt::TwoPointCorrelation_wedges::errorPerpendicular () const
{
  size_t sz = m_dataset->error_fx().size();
  vector<double> error_xiPerp;
  
  for (size_t i=0; i<sz/2; i++)
    error_xiPerp.push_back(m_dataset->error_fx()[i]);

  return error_xiPerp;
}


// ============================================================================================


vector<double> cosmobl::twopt::TwoPointCorrelation_wedges::xiParallel () const 
{
  size_t sz = m_dataset->fx().size();
  vector<double> xiParallel;

  for (size_t i=sz/2; i<sz; i++)
    xiParallel.push_back(m_dataset->fx()[i]);

  return xiParallel;
}


// ============================================================================================


vector<double> cosmobl::twopt::TwoPointCorrelation_wedges::errorParallel () const 
{
  size_t sz = m_dataset->error_fx().size();
  vector<double> error_xiParallel;

  for (size_t i=sz/2; i<sz; i++)
    error_xiParallel.push_back(m_dataset->error_fx()[i]);

  return error_xiParallel;
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::write (const string dir, const string file, const int rank) const 
{
  (void)rank;
  
  vector<double> rad = m_dataset->xx();
  vector<double> xiw = m_dataset->fx();
  vector<double> error = m_dataset->error_fx();

  checkDim(rad, m_dd->nbins_D1()*2, "rad");
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  string header = "[1] separation at the bin centre # [2] perpendicular wedge # [3] error on the perpendicular wedge # [4] parallel wedge # [5] error on the parallel wedge";
  if (m_compute_extra_info)
    header += " # [4] mean separation # [5] standard deviation of the separation distribution";
  
  fout << "### " << header << " ###" <<endl;

  for (int i=0; i<m_dd->nbins_D1(); i++) {
    fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << rad[i] << "  " << setw(8) << xiw[i] << "  " << setw(8) << error[i] << "  " << setw(8) << xiw[i+m_dd->nbins_D1()] << "  " << setw(8) << error[i+m_dd->nbins_D1()];
    if (m_compute_extra_info)
      for (size_t ex=0; ex<m_dataset->extra_info().size(); ++ex)
	fout << "  " << setw(8) << m_dataset->extra_info(ex, i);
    fout << endl;
  }
  
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}


// ============================================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation_wedges::Wedges (const vector<double> rr, const vector<double> mu, const vector<vector<double> > xi, const vector<vector<double> > error_xi)
{
  double binSize = mu[1]-mu[0];
  int muSize = mu.size();

  int half = max(0, min(int(0.5/binSize), muSize));
  half = (mu[half]<0.5) ? half+1 : half;
  
  vector<double> rad(2*rr.size(),0), wedges(2*rr.size(),0), error(2*rr.size(),0);

  for (size_t i=0; i<rr.size(); i++) {
    rad[i] = rr[i];
    rad[i+rr.size()] = rr[i];

    for (int j=0; j<half; j++) {
      wedges[i] += 2.*xi[i][j]*binSize;   	     // xi_perp
      error[i] += 2.*pow(error_xi[i][j]*binSize, 2); // error[xi_perp]
    }

    for (size_t j=half; j<mu.size(); j++) {
      wedges[i+rr.size()] += 2.*xi[i][j]*binSize;              // xi_par
      error[i+rr.size()] += 2.*pow(error_xi[i][j]*binSize, 2); // error[xi_par]
    }  
  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv);} );

  return (!m_compute_extra_info) ? unique_ptr<data::Data1D>(new data::Data1D(rad, wedges, error)) : data_with_extra_info(rad, wedges, error);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::measure(const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  switch (errorType) {
  case (ErrorType::_Poisson_) :
    measurePoisson(dir_output_pairs,dir_input_pairs,count_dd,count_rr,count_dr,tcount, estimator);
    break;
  case (ErrorType::_Jackknife_) :
    measureJackknife(dir_output_pairs,dir_input_pairs,dir_output_resample,count_dd,count_rr,count_dr,tcount, estimator);
    break;
  case (ErrorType::_Bootstrap_) :
    measureBootstrap(nMocks, dir_output_pairs,dir_input_pairs,dir_output_resample,count_dd,count_rr,count_dr,tcount, estimator);
    break;
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation2D_polar::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  
  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 
  
  m_dataset = Wedges(TwoPointCorrelation2D_polar::xx(), TwoPointCorrelation2D_polar::yy(), TwoPointCorrelation2D_polar::xi2D(), TwoPointCorrelation2D_polar::error2D());
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample != par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }

  vector<shared_ptr<data::Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_polar_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount, estimator);

  auto data_polar = (estimator==_natural_) ? NaturalEstimator(m_dd, m_rr, m_data->weightedN(), m_random->weightedN()) : LandySzalayEstimator(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN());

  if (estimator==_natural_) 
    data = XiJackknife(dd_regions, rr_regions);
  else if (estimator==_LandySzalay_)
    data = XiJackknife(dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation_wedges.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double> > ww, covariance;
  
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());
    
    if (dir_output_resample != par::defaultString) {
      string file_out = dir_output_resample+"xi_wedges_Jackknife_"+conv(i, par::fINT)+".dat";

      vector<double> rad = data[i]->xx();
      vector<double> xiw = data[i]->fx();
      vector<double> error = data[i]->error_fx();

      checkDim(rad, m_dd->nbins_D1()*2, "rad");
  
      ofstream fout(file_out.c_str()); checkIO(fout, file_out);

      fout << "### [1] separation at the bin centre # [2] perpendicular wedge # [3] error on the perpendicular wedge # [4] parallel wedge # [5] error on the parallel wedge ###" << endl;

      for (int i=0; i<m_dd->nbins_D1(); i++) 
	fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << rad[i] << "  " << setw(8) << xiw[i] << "  " << setw(8) << error[i] << "  " << setw(8) << xiw[i+m_dd->nbins_D1()] << "  " << setw(8) << error[i+m_dd->nbins_D1()] << endl;
   
      fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
    }
  }

  covariance_matrix(ww, covariance, true);

  m_dataset = Wedges(data_polar->xx(), data_polar->yy(), data_polar->fxy(), data_polar->error_fxy());
  m_dataset->set_covariance(covariance);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample != par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<shared_ptr<data::Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_polar_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount, estimator);

  auto data_polar = (estimator==_natural_) ? NaturalEstimator(m_dd, m_rr, m_data->weightedN(), m_random->weightedN()) : LandySzalayEstimator(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN());

  if (estimator==_natural_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions);
  else if (estimator==_LandySzalay_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation_wedges.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double> > ww, covariance;
  
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());
    
    if (dir_output_resample != par::defaultString) {
      string file_out = dir_output_resample+"xi_wedges_Bootstrap_"+conv(i, par::fINT)+".dat";

      vector<double> rad = data[i]->xx();
      vector<double> xiw = data[i]->fx();
      vector<double> error = data[i]->error_fx();

      checkDim(rad, m_dd->nbins_D1()*2, "rad");
  
      ofstream fout(file_out.c_str()); checkIO(fout, file_out);

      fout << "### [1] separation at the bin centre # [2] perpendicular wedge # [3] error on the perpendicular wedge # [4] parallel wedge # [5] error on the parallel wedge ###" << endl;

      for (int i=0; i<m_dd->nbins_D1(); i++) 
	fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << rad[i] << "  " << setw(8) << xiw[i] << "  " << setw(8) << error[i] << "  " << setw(8) << xiw[i+m_dd->nbins_D1()] << "  " << setw(8) << error[i+m_dd->nbins_D1()] << endl;
   
      fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
    }
  }
  covariance_matrix(ww, covariance, false);

  m_dataset = Wedges(data_polar->xx(), data_polar->yy(), data_polar->fxy(), data_polar->error_fxy());
  m_dataset->set_covariance(covariance);
}


// ============================================================================================


vector<shared_ptr<data::Data> > cosmobl::twopt::TwoPointCorrelation_wedges::XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<data::Data> > data;

  auto data2d = TwoPointCorrelation2D_polar::XiJackknife(dd, rr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(Wedges(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data> > cosmobl::twopt::TwoPointCorrelation_wedges::XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<data::Data> > data;
  
  auto data2d = TwoPointCorrelation2D_polar::XiJackknife(dd, rr, dr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(Wedges(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data> > cosmobl::twopt::TwoPointCorrelation_wedges::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<data::Data> > data;
  
  auto data2d = TwoPointCorrelation2D_polar::XiBootstrap(nMocks, dd, rr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(Wedges(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data> > cosmobl::twopt::TwoPointCorrelation_wedges::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<data::Data> > data;
  
  auto data2d = TwoPointCorrelation2D_polar::XiBootstrap(nMocks, dd, rr, dr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(Wedges(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  
  return data;
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::read_covariance_matrix (const string dir, const string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::write_covariance_matrix (const string dir, const string file) const
{
  m_dataset->write_covariance(dir, file, "r");
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::compute_covariance_matrix (const vector<shared_ptr<data::Data> > xi_collection, const bool doJK)
{
  vector<vector<double> > xi;

  for(size_t i=0;i<xi_collection.size();i++)
    xi.push_back(xi_collection[i]->fx());

  vector<vector<double> > cov_mat;
  cosmobl::covariance_matrix(xi, cov_mat, doJK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_wedges::compute_covariance_matrix (const vector<string> file_xi, const bool doJK)
{
  vector<double> rad, mean;
  vector<vector<double> > cov_mat;

  cosmobl::covariance_matrix (file_xi, rad, mean, cov_mat, doJK); 
  m_dataset->set_covariance(cov_mat);
}
