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
 *  @file Headers/Lib/TwoPointCorrelation2D.h
 *
 *  @brief The class TwoPointCorrelation2D
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D, used to measure the 2D two-point
 *  correlation function in Cartesian and polar coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2D__
#define __TWOPOINT2D__


#include "TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
    
    /**
     *  @class TwoPointCorrelation2D TwoPointCorrelation2D.h
     *  "Headers/Lib/TwoPointCorrelation2D.h"
     *
     *  @brief The class TwoPointCorrelation2D
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation2D </EM>. It is used to measure the 2D
     *  two-point correlation function, in Cartesian and polar
     *  coordinates, \f$\xi(r_p,\pi)\f$ or \f$\xi(r,\mu)\f$, that is as
     *  a function of perpendicular, \f$r_p\f$, and parallel, \f$\pi\f$,
     *  line-of-sight separations, and as a function of absolute
     *  separation, \f$r=\sqrt{r_p^2+\pi^2}\f$, and the cosine of the
     *  angle between the separation vector and the line of sight,
     *  \f$\mu\equiv\cos\theta=s_\parallel/s\f$, respectively.
     */
    class TwoPointCorrelation2D : public TwoPointCorrelation {

    protected :
    
      /**
       *  @name Internal input/output methods
       */
      ///@{
    
      /**
       *  @brief write the number of pairs
       *  @param PP pointer to an object of class Pair
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      void write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const override;

      /**
       *  @brief read the number of pairs
       *  @param [out] PP pointer to an object of class Pair
       *  @param [in] dir vector of input directories
       *  @param [in] file input file
       *  @return none
       */
      void read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const override;

      /**
       *  @brief write the number of pairs
       *  @param PP pointer to a vector of objects of class Pair
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      void write_pairs (const vector<shared_ptr<pairs::Pair> >  PP, const string dir, const string file) const override;

      /**
       *  @brief read the number of pairs
       *  @param [out] PP pointer to a vector of objects of class Pair
       *  @param [in] dir vector of input directories
       *  @param [in] file input file
       *  @return none
       */
      void read_pairs (vector<shared_ptr<pairs::Pair> > PP, const vector<string> dir, const string file) const override;
      
      ///@}


    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation2D
       */
      TwoPointCorrelation2D () { m_dataset = Data::Create(dataType::_2D_data_); }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @return object of class TwoPointCorrelation2D
       */
      TwoPointCorrelation2D (const catalogue::Catalogue data, const catalogue::Catalogue random) 
	: TwoPointCorrelation(data, random) { m_dataset = Data::Create(dataType::_2D_data_); }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation2D () = default;

      ///@}

        /**
       *  @name Member functions to get the private/protected members
       */
      ///@{

      /**
       *  @brief get the protected member Data::m_x
       *  @return the x coordinates
       */
      vector<double> xx () const override { return m_dataset->xx(); }

      /**
       *  @brief get the protected member m_y
       *  @return the y coordinates
       */
      vector<double> yy () const override { return m_dataset->yy(); }

      /**
       *  @brief get the protected member m_fxy
       *  @return the binned correlation function 
       */
      vector<vector<double> > xi2D () const override { return m_dataset->fxy(); }

      /**
       *  @brief get the protected member m_error_fxy
       *  @return the error on the binned correlation function
       *  function
       */
      vector<vector<double> > error2D () const override { return m_dataset->error_fxy(); }
      
      ///@}

      /**
       *  @name Member functions to count measure the two-point correlation function
       */
      ///@{

      /**
       *  @brief measure the xi with Poisson error using measured pairs
       *  using the Natural Estimator
       *  
       *  @param dd pointer to an object of type Pair containing the
       *  data-data pairs
       *
       *  @param rr pointer to an object of type Pair containing the
       *  random-random pairs
       *
       *  @param nData number of objects in the data catalogue
       *
       *  @param nRandom number of objects in the random catalogue
       *
       *  @return pointer to an object of type Data
       */
      shared_ptr<Data> NaturalEstimatorTwoP (shared_ptr<pairs::Pair> dd, shared_ptr<pairs::Pair> rr, const int nData, const int nRandom) override;

      /**
       *  @brief measure the xi with Poisson error using measured pairs
       *  using the Landy-Szalay estimator
       *  
       *  @param dd pointer to an object of type Pair containing the
       *  data-data pairs
       *
       *  @param rr pointer to an object of type Pair containing the
       *  random-random pairs
       *
       *  @param dr pointer to an object of type Pair containing the
       *  data-random pairs
       *
       *  @param nData number of objects in the data catalogue
       *
       *  @param nRandom number of objects in the random catalogue
       *
       *  @return pointer to an object of type Data
       */
      shared_ptr<Data> LandySzalayEstimatorTwoP (shared_ptr<pairs::Pair> dd, shared_ptr<pairs::Pair> rr, shared_ptr<pairs::Pair> dr, int nData, int nRandom) override;

      ///@}

    };
  }
}

#endif
