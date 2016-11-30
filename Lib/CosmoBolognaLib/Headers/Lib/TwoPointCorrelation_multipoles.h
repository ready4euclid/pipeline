/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Lib/TwoPointCorrelation_multipoles.h
 *
 *  @brief The class TwoPointCorrelation_multipoles
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation_multipoles, used to measure the first three
 *  multipole moments of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTMULT__
#define __TWOPOINTMULT__


#include "TwoPointCorrelation2D_polar.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
  
    /**
     *  @class TwoPointCorrelation_multipoles
     *  TwoPointCorrelation_multipoles.h
     *  "Headers/Lib/TwoPointCorrelation_multipoles.h"
     *
     *  @brief The class TwoPointCorrelation_multipoles
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation_multipoles </EM>. It is used to measure
     *  the first three multipole moments of the two-point correlation
     *  function, i.e. the monopole, the quadrupole and the
     *  exadecupole \f$\xi_l(r) = (2l+1) \int_{\mu_{min}}^{\mu_{max}}
     *  \xi(s,\mu) L_l(\mu) d\mu\f$, with \f$l=0,2,4\f$. If
     *  \f$\mu_{min}\neq0\f$ and \f$\mu_{max}\neq1\f$, the algorithm
     *  provides the truncated multipoles (
     *  http://arxiv.org/abs/1502.05045 ).
     *
     *  This class uses the so-called <EM> integrated </EM> method
     *  (e.g. http://arxiv.org/abs/1105.2037 , appendix E). This
     *  method is unbiased even in presence of wide angle and observer
     *  angle effects, i.e. if the random-random (RR) pairs depend
     *  also on \f$\mu\f$. However it is more affected by numerical
     *  uncertainties, relative to the <EM> direct </EM> method, due
     *  to the numerical integration. Moreover, it requires a large
     *  enough random catalogue, in order to have random pairs in each
     *  angular bin.
     *
     *  The monopole estimated with the <EM> direct </EM> method can be
     *  obtained with the class TwoPointCorrelation1D_monopole.
     */
    class TwoPointCorrelation_multipoles : public TwoPointCorrelation2D_polar {

    protected:
      
      /**
       *  @brief measure the multipoles of the two-point correlation
       *  function
       *  
       *  @param rr absolute separation 
       *
       *  @param mu angular separation
       *
       *  @param xi 2dD cartesian two-point correlation function
       *
       *  @param error_xi errors on the 2d polar two-point correlation
       *  function
       *
       *  @return pointer to an object of type Data
       */
      shared_ptr<data::Data> Multipoles (const vector<double> rr, const vector<double> mu, const vector<vector<double>> xi, const vector<vector<double>> error_xi) override;

      /**
       *  @brief return a data object with extra info
       *  
       *  @param rad vector containing the binned separations
       *
       *  @param xil vector containing the binned multipoles of the
       *  correlation function
       *
       *  @param error vector containing the errors
       *
       *  @return pointer to an object of type Data
       */
      shared_ptr<data::Data> data_with_extra_info (const vector<double> rad, const vector<double> xil, const vector<double> error) const;

      /**
       *  @brief measure the first three multipoles of the two-point
       *  correlation function with Poisson errors
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of
       *  random-random pairs; false &rarr; read the number of
       *  random-random pairs from file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      void measurePoisson (const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      /**
       *  @brief measure the first three multipoles of the two-point
       *  correlation function estimating the covariance with
       *  Jackknife resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_resample output directory used to store
       *  the Jackknife resampling correlation functions, with Poisson
       *  errors
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of
       *  random-random pairs; false &rarr; read the number of
       *  random-random pairs from file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      void measureJackknife (const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      /**
       *  @brief measure the first three multipoles of the two-point
       *  correlation function estimating the covariance with
       *  Bootstrap resampling
       *
       *  @param nMocks number of mocks to be generated with bootstrap
       *  resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_resample output directory used to store
       *  the Bootstrap resampling correlation function, with Poisson
       *  errors
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of
       *  random-random pairs; false &rarr; read the number of
       *  random-random pairs from file
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      void measureBootstrap (const int nMocks, const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      /**
       *  @brief measure the jackknife resampling of the first three
       *  multipoles of the two-point correlation function
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions
       *
       *  @return none
       */
      vector<shared_ptr<data::Data>> XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr) override;

      /**
       *  @brief measure the jackknife resampling of the first three
       *  multipoles of the two-point correlation function
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions
       *
       *  @param dr vector of random-random pairs, divided per regions  
       *
       *  @return a vector of pointers to objects of type Data
       */
      vector<shared_ptr<data::Data>> XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr) override;

      /**
       *  @brief measure the bootstrap resampling of the first three
       *  multipoles of the two-point correlation function
       * 
       *  @param nMocks number of bootstrap resamplings
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions
       *
       *  @return a vector of pointers to objects of type Data
       */
      vector<shared_ptr<data::Data>> XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr) override;

      /**
       *  @brief measure the bootstrap resampling of the two-point correlation
       *  function, &xi;(r)     
       *
       *  @param nMocks number of bootstrap resamplings
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions 
       *
       *  @param dr vector of random-random pairs, divided per regions  
       *
       *  @return a vector of pointers to objects of type Data
       */
      vector<shared_ptr<data::Data>> XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr) override;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation_multipoles
       */
      TwoPointCorrelation_multipoles () { m_twoPType = _1D_multipoles_; }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType_rad binning type in absolute separations
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param nbins_rad number of bins in the absolute
       *  separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param muMin minimum angular used to count the pairs
       *  @param muMax maximum angular used to count the pairs
       *  @param nbins_mu number of bins in the angular
       *  separation
       *  @param shift_mu shift parameter in the angular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation_multipoles (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false)
	: TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, nbins_rad, shift_rad, _linear_, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight, compute_extra_info)
	{ m_twoPType = _1D_multipoles_; }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType_rad binning type in absolute separations
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param binSize_rad bin size in the absolute separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param muMin minimum angular separation used to count
       *  the pairs
       *  @param muMax maximum angular separation used to count
       *  the pairs
       *  @param binSize_mu bin size in the angular separation
       *  @param shift_mu shift parameter in the angular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @param compute_extra_info true &rarr; compute extra
       *  information related to the pairs, such as the mean pair
       *  separation and redshift
       *  @return object of class TwoPointCorrelation_multipoles
       */
      TwoPointCorrelation_multipoles (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false)
	: TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, binSize_rad, shift_rad, _linear_, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight, compute_extra_info)
	{ m_twoPType = _1D_multipoles_; }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation_multipoles () = default;

      ///@}


      /**
       *  @name Member functions to count the number of pairs and measure the two-point correlation function
       */
      ///@{

      /**
       *  @brief get the x coordinates
       *  @return the x coordinates
       */
      vector<double> xx () const  override;

      /**
       *  @brief get the y coordinates
       *  @return the y coordinates
       */
      vector<double> yy () const 
        { cosmobl::ErrorCBL("Error in yy() of TwoPointCorrelation_multipoles.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the the binned correlation function 
       *  @return the binned correlation function 
       */
      vector<double> xi1D () const
        { cosmobl::ErrorCBL("Error in xi1D() of TwoPointCorrelation_multipoles.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the error on the binned correlation function
       *  function
       *  @return the error on the binned correlation function
       *  function
       */
      vector<double> error1D () const
        { cosmobl::ErrorCBL("Error in error1D() of TwoPointCorrelation_multipoles.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the the binned correlation function 
       *  @return the binned correlation function 
       */
      vector<vector<double>> xi2D () const 
        { cosmobl::ErrorCBL("Error in xi2D() of TwoPointCorrelation_multipoles.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the error on the binned correlation function
       *  function
       *  @return the error on the binned correlation function
       *  function
       */
      vector<vector<double>> error2D () const 
        { cosmobl::ErrorCBL("Error in error2D() of TwoPointCorrelation_multipoles.h!"); vector<vector<double>> vv; return vv; }

      /**
       *  @brief get the monopole of the polar xi
       *  @return the xiMonopole
       */
      vector<double> xiMonopole () const override;

      /**
       *  @brief get the error on the monopole of the polar xi
       *  @return the error on the Monopole
       */
      vector<double> errorMonopole () const override;

      /**
       *  @brief get the quadrupole of the polar xi
       *  @return the Quadrupole
       */
      vector<double> xiQuadrupole () const override;

      /**
       *  @brief get the error on the quadrupole of the polar xi
       *  @return the error on the Quadrupole
       */
      vector<double> errorQuadrupole () const override;

      /**
       *  @brief get the octupole of the polar xi
       *  @return the Octupole
       */
      vector<double> xiOctupole () const override;

      /**
       *  @brief get the error on the octupole of the polar xi
       *  @return the error on Octupole
       */
      vector<double> errorOctupole () const override;

      /**
       *  @brief measure the first three multipoles of the two-point
       *  correlation function
       *
       *  @param errorType type of error
       *  
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_resample output directory of the
       *  resampled correlation function
       *
       *  @param nMocks number of resampling used for bootstrap
       *
       *  @param count_dd true &rarr; count the number of data-data
       *  pairs; false &rarr; read the number of data-random pairs from
       *  file
       *
       *  @param count_rr true &rarr; count the number of random-random
       *  pairs; false &rarr; read the number of random-random pairs
       *
       *  @param count_dr true &rarr; count the number of data-random
       *  pairs; false &rarr; read the number of data-random pairs
       *
       *  @param tcount true &rarr; activate the time counter; false
       *  &rarr; no time counter
       *
       *  @param estimator the estimator used to measure the two-point
       *  correlation function
       *
       *  @return none
       */
      void measure (const ErrorType errorType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={},  const string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      ///@}


      /**
       *  @name Input/Output methods
       */
      ///@{

      /**
       *  @brief read the multipoles of the two-point correlation
       *  function
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read (const string dir, const string file) override
      { (void)dir; (void)file; ErrorCBL("Error in TwoPointCorrelation_multipoles::read of TwoPointCorrelation_multipoles.h: work in progress!", glob::ExitCode::_workInProgress_); }

      /**
       *  @brief write the multipoles of the two-point correlation
       *  function
       *  @param dir output directory
       *  @param file output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      void write (const string dir=par::defaultString, const string file=par::defaultString, const int rank=0) const override;

      ///@}

      /**
       *  @name Member functions to compute, read and write covariance matrix
       */
      ///@{ 

      /**
       *  @brief read the measured covariance matrix
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      virtual void read_covariance_matrix (const string dir, const string file) override;

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      virtual void write_covariance_matrix (const string dir, const string file) const override;

      /**
       *  @brief compute the covariance matrix
       *  @param xi_collection vector containing the xi to compute the covariance matrix
       *  @param doJK 1 &rarr; compute jackknife covariance matrix; 0 compute standard covariance matrix
       *  @return none
       */
      virtual void compute_covariance_matrix (vector<shared_ptr<data::Data>> xi_collection, bool doJK) override;

      /**
       *  @brief compute the covariance matrix
       *  @param file_xi vector containing the path to the xi to compute the covariance matrix
       *  @param doJK 1 &rarr; compute jackknife covariance matrix; 0 compute standard covariance matrix
       *  @return none
       */
      virtual void compute_covariance_matrix (vector<string> file_xi, bool doJK) override;

      ///@} 


    };
  }
}

#endif
