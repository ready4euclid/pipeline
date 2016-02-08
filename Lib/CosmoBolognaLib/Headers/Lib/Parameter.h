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
 *  @file Headers/Lib/Parameter.h
 *
 *  @brief The class Parameter
 *
 *  This file defines the interface of the class Parameter
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PARAM__
#define __PARAM__

#include "Prior.h"

namespace cosmobl {

  /**
   * @class Parameter Parameter.h "Headers/Lib/Parameter.h"
   *
   *  @brief The class Parameter
   *
   *  This class is used to define the parameters of models
   */
  class Parameter {

  protected:

    /// value: parameter value
    double m_value;
    /// best_value: parameter best value
    double m_best_value;
    /// proposed_value: parameter proposed value
    double m_proposed_value;

    /// name: parameter name
    string m_name;
    /// freeze: 0 &rarr; parameter is not fixed; 1 &rarr; parameter is fixed
    bool m_freeze;

    /// prioir: parameter prior
    shared_ptr<Prior> m_prior;

    /// nchains: number of chains 
    int m_nchains;
    /// chain_size: chain size
    int m_chain_size;
    /// chains: vectors containing the parameter chains
    vector< shared_ptr<Chain> > m_chains;
    
    /// mean: parameter mean value 
    double m_mean;
    /// standard deviation: parameter standard deviation 
    double m_std;
    /// median: parameter median value
    double m_median;

    /// var: binned parameter range
    vector<double> m_var;
    /// dist: parameter distribution
    vector<double> m_dist;

  public:
    /**
     *  @name Constructors/destructors
     */
    ///@{

    /**
     *  @brief default constructor
     *
     *  @return object of class Parameter
     */
    Parameter() {}

    /**
     *  @brief constructor
     *
     *  @param value: parameter value
     *  @param freeze: 0 &rarr; parameter is not fixed; 1 &rarr; parameter is fixed
     *  @param name: parameter name
     *
     *  @return object of class Parameter
     */
    Parameter(double, bool freeze=0, string name="parameter");

    /**
     *  @brief default constructor
     *
     *  @param value: parameter value
     *  @param limits: limits for the prior
     *  @param name: parameter name
     *
     *  @return object of class Parameter
     */
    Parameter(double , vector<double> , string name="parameter");

    /**
     *  @brief default constructor
     * 
     *  @param value: parameter value
     *  @param prior: prior for the parameter 
     *  @param name: parameter name
     *
     *  @return object of class Parameter
     */
    Parameter(double ,shared_ptr<Prior> , string name="parameter");

    /**
     *  @brief default destructor
     *
     *  @return none
     */
    ~Parameter() {}

    ///@}
    
    /**
     * @brief value of the prior at the parameter value
     *
     * @return prior value
     */
    double PriorProbability() {return m_prior->operator()(m_value);}

     /**
     * @brief value of the prior at a proposed value
     *
     * @param value: proporsed parameter value
     *
     * @return prior value
     */ 
    double PriorProbability(double value ) {return m_prior->operator()(value);}

    /**
     * @brief set the private attribute m_value
     *
     * @param value: parameter value
     *
     * @return none
     */
    void set_value (double value) {m_value = (m_freeze) ? m_value : value;}

    /**
     * @brief set the private attribute m_best_value
     *
     * @param value: parameter best value
     *
     * @return none
     */

    void set_best_value (double value) {m_best_value = value;}

    /**
     * @brief set the private attribute m_proposed_value
     *
     * @param value: parameter proposed value
     *
     * @return none
     */
    void set_proposed_value (double value) {m_proposed_value = value;}

    /**
     * @brief set the private attribute m_name
     *
     * @param name: parameter name
     *
     * @return none
     */
    void set_name (string name) {m_name = name;}

    /**
     * @brief lock the parameter
     *
     * @param freeze: 0 &rarr; parameter is not locked; 1 &rarr; parameter is locked
     * 
     * @return none
     */
    void set_freeze (bool freeze) {m_freeze = freeze;}

    /**
     * @brief set rectangular prior
     * 
     * @param limits: prior limits
     * 
     * @return none
     */
    void set_prior (vector<double> limits) {m_prior = make_shared<Prior>(limits);} 

    /**
     * @brief set gaussian prior
     * 
     * @param mean: gaussian prior mean value
     * @param sigma: gaussian prior standard deviation value
     * 
     * @return none
     */
    void set_prior (double mean, double sigma) {m_prior = make_shared<Prior>(mean,sigma);}

    /**
     * @brief set user defined prior
     * 
     * @param prior: user defined prior
     * 
     * @return none
     */
    void set_prior (shared_ptr<Prior> prior) {m_prior = prior;}

    /**
     * @brief set user defined prior
     * 
     * @return shared pointer to an object of class Prior
     */
    shared_ptr<Prior> prior() {return m_prior;}

    /**
     * @brief function that return the prior interval lenght
     * 
     * @return the prior interval lenght
     */
    double interval_size() {return (m_prior->xmax()-m_prior->xmin());}

    /**
     * @brief get the private member m_value
     *
     * @return value of the parameter
     */
    double value() {return m_value;}

    /**
     * @brief get the private member m_proposed_value
     *
     * @return proposed value of the parameter
     */  
    double proposed_value() {return m_proposed_value;}
    
    /**
     * @brief get the private member m_best_value
     *
     * @return best value of the parameter
     */  
    double best_value() {return m_best_value;}

    /**
     * @brief get the private member m_name
     *
     * @return return the parameter name
     */
    string name() { return m_name;}

    /**
     * @brief get the private member m_freeze
     *
     * @return return 0 &rarr; parameter is not locked; 1 &rarr; parameter is locked
     */
    bool freeze() { return m_freeze;}

    /**
     * @brief evaluate priors ratio between proposed and current parameter value
     *
     * @return prior ratio
     */
    double eval_proposed();

    /**
     * @brief evaluate prior ratio between proposed and current parameter value
     * 
     * @param proposed_value: proposed value of the parameter
     *
     * @return prior ratio
     */ 
    double eval_proposed(double );
   
    /**
     * @brief assign to the current parameter value the proposed
     *
     * @return none
     */
    void confirm_proposed_value();

    /**
     * @brief set the chains
     *
     * @param nchains: the number of chains
     * @param chain_size: the chains size
     *
     * @return none
     */
    void set_chains(int , int);

    /**
     * @brief get the i-th chain
     *
     * @param i: the i-th chain
     *
     * @return the i-th chain 
     */
    shared_ptr<Chain> chain(int i) { return m_chains[i];}

    /**
     * @brief get all the chains
     *
     * @return a vector containing all the chains
     */
    vector<shared_ptr<Chain> > chains() { return m_chains;}

    /**
     * @brief merge the chains
     *
     * @param max: maximum step
     * @param min: minimum step
     * @param thin: pick values every thin times
     *
     * @return merged chain
     */
    shared_ptr<Chain> merge_chains(int max = -1, int min=-1, int thin = 1);

    /**
     * @brief get the private member m_mean
     *
     * @return the mean value of the parameter
     */
    double mean () {return m_median;};
    
    /**
     * @brief get the private member m_median
     *
     * @return the median value of the parameter
     */   
    double median() {return m_median;};
    
    /**
     * @brief get the private member m_std
     *
     * @return the standard deviation value of the parameter
     */
    double std() {return m_std;}

    /**
     * @brief generate a random value for the parameter in the prior
     * range
     *
     * @param random: random number between 0 and 1
     * 
     * @return the parameter random value
     *
     */
    double random_value(double ); 
  };

}

#endif
