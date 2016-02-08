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
 *  @file Headers/Lib/TwoPointCorrelation_deprojected.h
 *
 *  @brief The class TwoPointCorrelation_deprojected
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation_deprojected, used to measure the deprojected
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTDEPROJ__
#define __TWOPOINTDEPROJ__


#include "TwoPointCorrelation_projected.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
  
    /**
     *  @class TwoPointCorrelation_deprojected TwoPointCorrelation_deprojected.h
     *  "Headers/Lib/TwoPointCorrelation_deprojected.h"
     *
     *  @brief The class TwoPointCorrelation_deprojected
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation_deprojected </EM>. It is used to measure
     *  the deprojected two-point correlation function,
     *  \f$\xi(r)=-\frac{1}{\pi}\int^{r_{\rm max}}_r dr_p'
     *  \frac{dw_p(r_p')/dr_p}{\sqrt{r_p'^2-r^2}}\f$.
     */
    class TwoPointCorrelation_deprojected : public TwoPointCorrelation_projected {
    
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation_deprojected
       */
      TwoPointCorrelation_deprojected () { m_twoPType = _1D_deprojected_; }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param rpMin minimum perpendicular separation used to count
       *  the pairs
       *  @param rpMax maximum perpendicular separation used to count
       *  the pairs
       *  @param nbins_rp number of bins in the perpendicular
       *  separation
       *  @param shift_rp shift parameter in the perpendicular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param piMin minimum parallel separation used to count
       *  the pairs
       *  @param piMax maximum parallel separation used to count
       *  the pairs
       *  @param nbins_pi number of bins in the parallel
       *  separation
       *  @param shift_pi shift parameter in the parallel
       *  separation, i.e. the radial shift is binSize*shift
       *  @return object of class TwoPointCorrelation_deprojected
       */
      TwoPointCorrelation_deprojected (catalogue::Catalogue data, catalogue::Catalogue random, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi)
	: TwoPointCorrelation_projected(data, random, _logarithmic_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi)
	{ m_twoPType = _1D_deprojected_; }
      
      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param rpMin minimum perpendicular separation used to count
       *  the pairs
       *  @param rpMax maximum perpendicular separation used to count
       *  the pairs
       *  @param binSize_rp bin size in the perpendicular separation
       *  @param shift_rp shift parameter in the perpendicular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param piMin minimum parallel separation used to count
       *  the pairs
       *  @param piMax maximum parallel separation used to count
       *  the pairs
       *  @param binSize_pi bin size in the parallel separation
       *  @param shift_pi shift parameter in the parallel
       *  separation, i.e. the radial shift is binSize*shift
       *  @return object of class TwoPointCorrelation_2D_deprojected
       */
      TwoPointCorrelation_deprojected (catalogue::Catalogue data, catalogue::Catalogue random, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi)
	: TwoPointCorrelation_projected(data, random, _logarithmic_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi)
	{ m_twoPType = _1D_deprojected_; }
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation_deprojected () = default;

      ///@}

      
      /**
       *  @name Member functions to count the number of pairs and measure the two-point correlation function
       */
      ///@{

      /**
       *  @brief measure the monopole of the two-point correlation
       *  function, &xi;(r)
       *  @param piMax_integral upper limit of the integral
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
      void measure (const double piMax_integral=50., const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1);
    
      ///@}

    
      /**
       *  @name Input/Output methods
       */
      ///@{

      /**
       *  @brief read the deprojected two-point correlation function
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read (const string dir, const string file) override;
      
      /**
       *  @brief write the deprojected two-point correlation function
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