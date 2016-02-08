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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation.h
 *
 *  @brief The class Modelling_TwoPointCorrelation
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation, used for modelling two-point
 *  correlation functions of any type
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLINGTWOPT__
#define __MODELLINGTWOPT__


#include "Likelihood.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation
     *  Modelling_TwoPointCorrelation.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation, used for modelling two-point
     *  correlation functions of any type
     *
     */
    class Modelling_TwoPointCorrelation : public Modelling {
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Modelling_TwoPointCorrelation
       */
      Modelling_TwoPointCorrelation () {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation () {}

      ///@}
      
    };
  }
}

#endif
