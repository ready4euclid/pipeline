/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Lib/Prior.h
 *
 *  @brief The class Prior 
 *
 *  This file defines the interface of the class Prior
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PRIOR__
#define __PRIOR__

#include "Chain.h"


// ============================================================================================


namespace cosmobl {

  typedef function<double(double, shared_ptr<void>, vector<double>)> prior_func;

  /**
   *  @class Prior Prior.h "Headers/Lib/Prior.h"
   *
   *  @brief The class Prior
   *
   *  This class is used to define the prior
   */
  class Prior {
      
  protected:
      
    /// the prior function
    prior_func m_func;
    
    /// parameters of the prior func
    vector<double> m_prior_func_pars;
    
    /// the prior lower limit
    double m_xmin;
    
    /// the prior upper limit
    double m_xmax;

    
  public:
    
    /**
     *  @name Constructors/destructors
     */
    ///@{

    /**
     *  @brief default constructor
     *
     *  @return object of class Prior
     */
    Prior (); 

    /**
     *  @brief constructor
     * 
     *  @param Limits: limits of the Prior
     *
     *  @return object of class Prior
     */
    Prior (vector<double>);

    /**
     *  @brief constructor
     * 
     *  @param mean the gaussian prior mean
     *  @param sigma gaussian prior standard deviation
     *
     *  @return object of class Prior
     */
    Prior (double , double );

    /**
     *  @brief constructor
     * 
     *  @param Limits limits of the Prior
     *  @param mean the gaussian prior mean
     *  @param sigma gaussian prior standard deviation
     *
     *  @return object of class Prior
     */
    Prior (vector<double>, double , double );

    /**
     *  @brief constructor
     * 
     *  @param Limits limits of the Prior
     *  @param func prior function
     *  @param prior_func_params prior function parameters
     *
     *  @return object of class Prior
     */
    Prior (vector<double> Limits, prior_func func, vector<double> prior_func_params);

    /**
     *  @brief default destructor
     *
     *  @return none
     */
    ~Prior () {} 

    ///@}
    
    /**
     * @brief evaluate prior 
     *
     * @param xx the value for prior calculation
     *
     * @return the prior value
     */
    double operator() (double xx) {
      shared_ptr<void> pp = NULL;
      if (xx<m_xmin || xx>m_xmax) return 1.e30;
      else return m_func(xx, pp, m_prior_func_pars);
    }
    
    /**
     * @brief set the prior limits 
     *
     * @param Limits the prior limits
     *
     * @return none
     */
    void set_limits (vector<double>);
    
    /**
     * @brief set parameters for gaussian prior 
     *
     * @param mean the gaussian prior mean
     * @param sigma the gaussian prior standard deviation
     *
     * @return none
     */
    void set_gaussian_parameters (double mean, double sigma);
    
    /**
     * @brief set prior function 
     *
     * @param func the prior function
     * @param func_parameters the prior function parameters
     *
     * @return none
     */       
    void set_func_parameters (prior_func func, vector<double> func_parameters);
    
    /**
     * @brief return the private member m_xmin
     *
     * @return the prior lower limit
     */
    double xmin() {return m_xmin;}
    
    /**
     * @brief return the private member m_xmax
     *
     * @return the prior upper limit
     */
    double xmax() {return m_xmax;}
    
    /**
     * @brief check if a value is included in the prior limits
     *
     * @return 0 &rarr; not included in prior range; 1 &rarr; included in prior range
     */
    bool isIncluded (double);
  };
}

#endif
