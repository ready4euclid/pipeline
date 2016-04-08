/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_projected.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_projected
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_projected, used for modelling projected 2pcf
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLINGPROJ__
#define __MODELLINGPROJ__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of functions and classes used for modelling
   *  
   * The \e modelling namespace contains all the functions and classes
   * used to model any kind of measurements
   */
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_projected
     *  Modelling_TwoPointCorrelation_projected.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_projected.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_projected
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation_projected, used to model the
     *  projected of the two-point correlation function
     *
     */
    class Modelling_TwoPointCorrelation_projected : public Modelling_TwoPointCorrelation {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class ModellingTwoPointCorrelation_projected
	 */
	Modelling_TwoPointCorrelation_projected () {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation_projected () {}

	/**
	 *  @brief constructor of the Modelling_TwoPointCorrelation_projected
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @param redshift the redshift of the two-point correlation
	 *  signal
	 *
	 *  @param cosmology the fiducial cosmology
	 *
	 *  @return object of type Modelling_TwoPointCorrelation_projected
	 */
	Modelling_TwoPointCorrelation_projected(const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop, const double redshift, const Cosmology cosmology);
	
	///@}

	void fit_bias(const string LikelihoodType, const vector<double> xlimits, const double bias_value, const statistics::Prior bias_prior, const int nChains, const int chain_size, const string dir_output, const double start=0.5, const double stop=1, const int thin=1) override;

	void fit_bias_cosmology(const string LikelihoodType, const vector<double> xlimits, const double bias_value, const statistics::Prior bias_prior, const vector<CosmoPar> CosmoPars, const vector<statistics::Prior> prior_CosmoPars, const int nChains, const int chain_size, const string dir_output, const double start=0.5, const double stop=1, const int thin=1) override;

    };
  }
}

#endif
