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
 *  @file Headers/Lib/TwoPointCorrelation2D_polar.h
 *
 *  @brief The class TwoPointCorrelation2D_polar
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D_polar, used to measure the 2D two-point
 *  correlation function in polar coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2DPOL__
#define __TWOPOINT2DPOL__


#include "TwoPointCorrelation2D.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
    
    /**
     *  @class TwoPointCorrelation2D_polar TwoPointCorrelation2D_polar.h
     *  "Headers/Lib/TwoPointCorrelation2D_polar.h"
     *
     *  @brief The class TwoPointCorrelation2D_polar
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation2D_polar </EM>. It is used to measure the
     *  2D two-point correlation function in polar coordinates,
     *  \f$\xi(r,\mu)\f$, that is as a function of absolute
     *  separation, \f$r=\sqrt{r_p^2+\pi^2}\f$, and the cosine of the
     *  angle between the separation vector and the line of sight,
     *  \f$\mu\equiv\cos\theta=s_\parallel/s\f$.
     */
    class TwoPointCorrelation2D_polar : public TwoPointCorrelation2D {

    protected:

      /**
       *  @brief measure the 2D two-point correlation function in
       *  polar coordinates, with Poisson errors
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
      void measurePoisson (const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      /**
       *  @brief measure the 2D two-point correlation function in
       *  polar coordinates, estimating the covariance with Jackknife
       *  resampling
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
      void measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs={}, const string dir_output_resample = "NULL", const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      /**
       *  @brief measure the 2D two-point correlation function in
       *  polar coordinates, estimating the covariance with Bootstrap
       *  resampling
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
      void measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs={}, const string dir_output_resample="NULL", const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation2D_polar () { m_twoPType = _2D_polar_; }

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
       *  @param binType_mu binning type in angular separations
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
      TwoPointCorrelation2D_polar (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false)
	: TwoPointCorrelation2D(data, random, compute_extra_info) { m_twoPType = _2D_polar_; set_parameters(binType_rad, rMin, rMax, nbins_rad, shift_rad, binType_mu, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight, compute_extra_info); }

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
       *  @param binType_mu binning type in angular separations
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
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation2D_polar (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false)
	: TwoPointCorrelation2D(data, random, compute_extra_info) { m_twoPType = _2D_polar_; set_parameters(binType_rad, rMin, rMax, binSize_rad, shift_rad, binType_mu, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight, compute_extra_info); }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation2D_polar () = default;

      ///@}


      /**
       *  @name Member functions to set the binning parameters
       */
      ///@{

      /**
       *  @brief set the binning parameters
       *  @param binType_rad binning type in absolute separations
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param nbins_rad number of bins in the absolute
       *  separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param binType_mu binning type in angular separations
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
       *  @return none
       */
      void set_parameters (const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

      /**
       *  @brief set the binning parameters
       *  @param binType_rad binning type in absolute separations
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param binSize_rad bin size in the absolute separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param binType_mu binning type in angular separations
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
       *  @return none
       */
      void set_parameters (const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

      ///@}


      /**
       *  @name Methods to set the binning parameters
       */
      ///@{

      /**
       *  @brief measure the 2D two-point correlation function in
       *  polar coordinates
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
      void measure (const ErrorType errorType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      ///@}


      /**
       *  @name Input/Output methods
       */
      ///@{

      /**
       *  @brief read the 2D two-point correlation function
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read (const string dir, const string file) override;     

      /**
       *  @brief write the 2D two-point correlation function
       *  @param dir output directory
       *  @param file output file
       *  @param full false &rarr; simply store the data; true &rarr;
       *  duplicate the data in the other three quadrands (usefull
       *  e.g. when storing the 2D correlation function)
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      void write (const string dir, const string file, const bool full, const int rank=0) const override;

      /**
       *  @brief write the 2D two-point correlation function
       *  @param dir output directory
       *  @param file output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      void write (const string dir=par::defaultString, const string file=par::defaultString, const int rank=0) const override
      { write(dir, file, true, rank); }
      
      ///@}

    };
  }
}

#endif
