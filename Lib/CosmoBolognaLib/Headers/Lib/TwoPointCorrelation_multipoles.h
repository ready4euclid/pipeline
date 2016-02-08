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
 *  multipoles of the two-point correlation function
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
     *  @class TwoPointCorrelation_multipoles TwoPointCorrelation_multipoles.h
     *  "Headers/Lib/TwoPointCorrelation_multipoles.h"
     *
     *  @brief The class TwoPointCorrelation_multipoles
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation_multipoles </EM>. It is used to measure
     *  the first three multipoles of the two-point correlation
     *  function, \f$\xi_l(r) = \frac{2l+1}{2} \int_{-1}^{1}
     *  \xi(s,\mu) L_l(\mu) d\mu\f$, with $l=0,2,4$.
     */
    class TwoPointCorrelation_multipoles : public TwoPointCorrelation2D_polar {

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
       *  @param binType_rad binning type in absolute separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
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
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation_multipoles (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const double muMin, const double muMax, const int nbins_mu, const double shift_mu)
	: TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, nbins_rad, shift_rad, _linear_, muMin, muMax, nbins_mu, shift_mu)
	{ m_twoPType = _1D_multipoles_; }
      
      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType_rad binning type in absolute separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
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
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation_multipoles (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const double muMin, const double muMax, const double binSize_mu, const double shift_mu)
	: TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, binSize_rad, shift_rad, _linear_, muMin, muMax, binSize_mu, shift_mu)
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
       *  @brief measure the first three monopoles of the two-point
       *  correlation function, &xi;<SUB>l</SUB>(r)
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *  @param count_dd 1 &rarr; count the number of data-data
       *  opairs; 0 &rarr; read the number of data-data pairs from
       *  file
       *  @param count_rr 1 &rarr; count the number of random-random
       *  opairs; 0 &rarr; read the number of random-random pairs from
       *  file
       *  @param count_dr 1 &rarr; count the number of data-random
       *  opairs; 0 &rarr; read the number of data-random pairs from
       *  file
       *  @param tcount 1 &rarr; activate the time counter; 0 &rarr;
       *  don't activate the time counter; 
       *  @return none
       */
      void measure (const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1);
    
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
      { ErrorMsg("Error in TwoPointCorrelation_multipoles::read of TwoPointCorrelation_multipoles.h: work in progress!"); }

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

    };
  }
}

#endif
