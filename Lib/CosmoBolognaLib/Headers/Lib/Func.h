/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/Lib/Func.h
 *
 *  @brief Useful generic functions
 *
 *  This file contains the prototypes of a large set of useful
 *  functions of wide used
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __FUNC__
#define __FUNC__

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
#include <limits>
#include <algorithm>
#include <memory>
#include <numeric>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <map>
#include <omp.h>

#ifdef LINUX
#include "sys/types.h"
#include "sys/sysinfo.h"
#endif

/// @cond GSLinc
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
/// @endcond

/// @cond NRinc
#include <nr3.h>
#include <ran.h>
#include <interp_1d.h>
#include <interp_linear.h>
#include <interp_2d.h>
#include <quadrature.h>
#include <romberg.h>
#include <roots.h>
#include <sort.h>
#include <gamma.h>
#include <deviates.h>
#include <mins.h>
#include <mins_ndim.h>
#include <fitmed.h>
#include <moment.h>
#include <gaussj.h>
#include <ludcmp.h>
#include <dfridr.h>
#include <fitlin.h>
/// @endcond


#include "Constants.h"

using namespace std;


// ============================================================================================


/**
 *  @example vectors.cpp 
 *
 *  This example shows how to print vectors and matrices, and how to
 *  remove elements
 */
/**
 *  @example distances.cpp  
 *
 *  This example shows how to convert redshifts into comoving
 *  distances using the class cosmobl::Cosmology
 */
/**
 *  @example fsigma8.cpp
 *
 *  This example shows how to estimate f*sigma8(z=1) using the class
 *  cosmobl::Cosmology
 */
/**
 *  @example chi2.cpp
 *
 *  This example shows how to compute the chi2 using the class
 *  cosmobl::Chi2
 */
/**
 *  @example likelihood.cpp
 *
 *  This example shows how to compute the likelihood using the class
 *  cosmobl::Likelihood
 */
/**
 * @example 2pt_monopole.cpp 
 *
 * This example shows how to measure the monopole of the two-point
 * correlation function, estimating the errors analytically, using the
 * class cosmobl::TwoPointCorrelation
 */
/**
 * @example 2pt_2D.cpp 
 *
 * This example shows how to measure the 2D two-point correlation
 * function, estimating the errors analytically, using the class
 * cosmobl::TwoPointCorrelation
 */
/**
 * @example 2pt_jackknife.cpp 
 *
 * This example shows how to measure the two-point correlation
 * function, estimating the errors with the jackknife method, using the
 * class cosmobl::TwoPointCorrelation
 */
/**
 * @example 2pt_projected_jackknife.cpp 
 *
 * This example shows how to measure the projected two-point
 * correlation function, estimating the errors with the jackknife
 * method, using the class cosmobl::TwoPointCorrelation
 */
/**
 * @example 3pt.cpp 
 *
 * This example shows how to measure the three-point correlation
 * function using the class cosmobl::ThreePointCorrelation
 */
/**
 *  @example distances.py 
 *
 *  This example shows how to convert redshifts into comoving
 *  distances using the class cosmobl::Cosmology. The CosmoBolognaLib
 *  are included as Python modules
 */
/**
 *  @example 2pt_model.py
 *
 *  This example shows how to compute power spectrum and two-point
 *  correlation function models using the class
 *  cosmobl::Cosmology. The CosmoBolognaLib are included as Python
 *  modules
 */
/**
 *  @example 2pt_monopole.py
 *
 *  This example shows how to compute the monopole of the two-point
 *  correlation function, estimating the errors analytically, using
 *  the class cosmobl::TwoPointCorrelation. The CosmoBolognaLib are
 *  included as Python modules
 */

/**
 *  @brief The namespace of the CosmoBolognaLib 
 *  
 * The \e cosmobl namespace contains all the main functions and
 * classes of the CosmoBolognaLib
 */
namespace cosmobl {

  /**
   * @enum Dim
   * @brief the dimension, used e.g. for pair and triplet vectors
   */
  enum Dim {
    
    /// 1D, used e.g. for 1D pairs, in angular or comoving separations
    _1D_,
    
    /// 2D pair, used e.g. for 2D pairs, in Cartesian or polar coordinates
    _2D_
    
  };
  

  /**
   * @enum binType
   * @brief the binning type
   */
  enum binType { 

    /// linear binning
    _linear_,
      
    /// logarithmic binning
    _logarithmic_
      
  };
    

  /**
   *  @name Functions of generic use  
   */
  ///@{

  /**
   *  @brief warning message
   *  @param msg string containing the warning message
   *  @return none
   */
  inline void WarningMsg (const string msg) { cerr << par::col_blue << msg << cosmobl::par::col_default << endl; }

  /**
   *  @brief error message
   *  @param msg string containing the Error message
   *  @return none
   *  @warning this function terminates the programs
   */
  inline void ErrorMsg (const string msg) { cerr << par::col_red << msg << cosmobl::par::col_default << endl; exit(1); }
  
  /**
   *  @brief produce a beep using the software totem
   *  @return none
   */
  inline void Beep ()
  {
    string beep = "totem "+par::DirCosmo+"Func/beep.mp3";
    if (system (beep.c_str())) {}; 
  }

  /**
   *  @brief terminate the process producing a beep with the software
   *  totem 
   *  @return none
   */
  inline void Exit ()
  {
    Beep();
    exit(1);
  }

  /**
   *  @brief check if the value of a [double] variable has already
   *  been set
   *
   *  @param var a double variable
   *  @return if var<-9.e29 &rArr; 0; else &rArr; 1
   */
  inline bool isSet (const double var) 
  {
    return (var>par::defaultDouble) ? 1 : 0;
  }
  
  /**
   *  @brief check if the values of a [double] vector have already
   *  been set
   *
   *  @param vect a vactor of double values
   *  @return if vect[i]<-9.e29 &forall; i &rArr; 0; else &rArr; 1
   */
  inline bool isSet (const vector<double> vect) 
  {
    bool is = 1;
    size_t ind = 0;
    while (is && ind<vect.size())
      if (vect[ind++]<-9.e29) is = 0;
    return is;
  }

  /**
   *  @brief convert a number to a string
   *  @param val number of any type
   *  @param fact output format 
   *  @return a string containing T
   */
  template <typename T> string conv (const T val, const char *fact) {
    char VAL[20]; sprintf (VAL, fact, val); 
    return string(VAL);
  }
  
  /**
   *  @brief the nearest integer
   *  @param val a number
   *  @return the integer value nearest to val
   */
  template <typename T> 
    int nint (const T val) 
    {
      return (val<0) ? val-0.5 : val+0.5;
    }

  /**
   *  @brief common logarithm (i.e. logarithm to base 10)
   *  @param val a number
   *  @return if val>0 &rarr; log10(val); else &rarr;
   *  par::defaultDouble
   */
  template <typename T> 
    T Log (const T val) 
    {
      return (val>0) ? log10(val) : par::defaultDouble;
    }

  /**
   *  @brief endian conversion of a short variable
   *  @param s a short variable
   *  @return the converted variable s
   */
  short ShortSwap (const short);

  /**
   *  @brief endian conversion of an integer variable
   *  @param i an integer variable
   *  @return the converted variable i
   */
  int IntSwap (const int);

  /**
   *  @brief endian conversion of a long integer variable
   *  @param i a long integer variable
   *  @return the converted variable i
   */
  long LongSwap (const long);

  /**
   *  @brief endian conversion of a float variable
   *  @param f a flot variable
   *  @return the converted variable f
   */
  float FloatSwap (const float);

  /**
   *  @brief endian conversion of a double variable
   *  @param d a double variable
   *  @return the converted variable d
   */
  double DoubleSwap (const double);
  
  /**
   *  @brief 1D interpolation
   *
   *  @param [in] _xx the point where the input function will be
   *  interpolated or extrapolated
   *
   *  @param [in] xx vector containing the binned values of x
   *
   *  @param [in] yy vector containing the binned values of the
   *  function, y(x), to be interpolated or extrapolated
   *
   *  @param [in] type the method used to interpolate or extrapolate:
   *  "Linear" &rarr; linear interpolation; "Poly" &rarr; polynomial
   *  interpolation; "Spline" &rarr; cubic spline interpolation; "Rat"
   *  &rarr; diagonal rational function interpolation; "BaryRat"
   *  &rarr; barycentric rational interpolation
   *
   *  @param [in] nPt number of points used in the polynomial,
   *  diagonal rational and barycentric rational interpolations
   *
   *  @param [out] err the estimated error in the interpolation or
   *  extrapolation
   *
   *  @return the interpolated value of the input function
   *
   *  @warning if _xx is outside the range of the input vector xx, the
   *  returned value is the extrapolation
   */
  double interpolated (const double, const vector<double>, const vector<double>, const string, const int, double &);
  
   /**
   *  @brief 1D interpolation
   *
   *  @param [in] _xx the point where the input function will be
   *  interpolated
   *
   *  @param [in] xx vector containing the binned values of x
   *
   *  @param [in] yy vector containing the binned values of the
   *  function, y(x), to be interpolated or extrapolated
   *
   *  @param [in] type the method used to interpolate or extrapolate:
   *  "Linear" &rarr; linear interpolation; "Poly" &rarr; polynomial
   *  interpolation; "Spline" &rarr; cubic spline interpolation; "Rat"
   *  &rarr; diagonal rational function interpolation; "BaryRat"
   *  &rarr; barycentric rational interpolation
   *
   *  @param [in] nPt number of points used in the polynomial,
   *  diagonal rational and barycentric rational interpolations
   *
   *  @return the interpolated or extrapolated value of the
   *  input function
   *
   *  @warning if _xx is outside the range of the input vector xx, the
   *  returned value is the extrapolation
   */
  double interpolated (const double, const vector<double>, const vector<double>, const string, const int);
  
  /**
   *  @brief 2D interpolation
   *
   *  @param [in] _x1 the point in the first dimension where the input
   *  function will be interpolated
   *
   *  @param [in] _x2 the point in the second dimension where the
   *  input function will be interpolated
   *
   *  @param [in] x1 vector containing the binned values of x in the
   *  first dimension
   *
   *  @param [in] x2 vector containing the binned values of x in the
   *  second dimension
   *
   *  @param [in] yy vector containing the binned values of the
   *  function, y(x), to be interpolated or extrapolated
   *
   *  @param [in] type the method used to interpolate or extrapolate:
   *  "Linear" &rarr; linear interpolation; "Poly" &rarr; polynomial
   *  interpolation; "Spline" &rarr; cubic spline interpolation; "Rat"
   *  &rarr; diagonal rational function interpolation; "BaryRat"
   *  &rarr; barycentric rational interpolation
   *
   *  @param [in] nPt number of points used in the polynomial,
   *  diagonal rational and barycentric rational interpolations
   *
   *  @return the interpolated or extrapolated value of the
   *  input function
   *
   *  @warning if _x1 and/or _x2 are outside the range of the input
   *  vectors x1 and/or x2, the returned value is the extrapolatation
   *
   */
  double interpolated_2D (const double, const double, const vector<double>, const vector<double>, const vector<vector<double> >, const string, const int);
  
  /**
   *  @brief check if a file can be opened
   *  @param file the file name
   *  @param isInput 1 &rarr; the file is open for input operations; 0
   *  &rarr; the file is open for output operations
   *  @return none
   */
  void checkIO (const string, const bool);

  /**
   *  @brief convert a map from a gnuplot file to a SM file 
   *  @param file_gnu input gnuplot file
   *  @param file_sm output SM file  
   *  @return none
   */
  void convert_map_gnuplot_sm (const string, const string);

  /**
   *  @brief set evironment variables
   *  @param Var vector containing the evironment variables to be set
   *  @return none
   */
  void set_EnvVar (const vector<string>);

  /**
   *  @brief check if an environment variable exists
   *  @param Var the evironment variable to be checked
   *  @return none
   */
  void check_EnvVar (const string); 

  /**
   *  @brief get the memory used by current process in kB
   *
   *  @warning this function works only on Linux systems
   *
   *  @param type 1 &rarr; Physical Memory (RAM); 2 &rarr; Virtual
   *  Memory
   *
   *  @return the Physical (RAM) or Virtual Memory used by current
   *  process in kB
   */
  int used_memory (const int);

  /**
   *  @brief check if the memory used by current process is larger
   *  than a given fraction of the available memory
   *
   *  @warning this function works only on Linux systems
   *
   *  @param frac the fraction of the available memory that is allowed
   *
   *  @param func a string that should contain the name of the
   *  function from which check_memory is called; it is used when
   *  printing the error message
   *
   *  @param exit 0 &rarr; warning message; 1 &rarr; error message;
   *  (and exit)
   *
   *  @param type 1 &rarr; Physical Memory (RAM); 2 &rarr; Virtual
   *  Memory
   *
   *  @return 0 &rarr; memory problems; 1 &rarr; no memory problems 
   */
  int check_memory (const double, const bool exit=1, const string func="", const int type=1);
  

  // ============================================================================================

  
  /**
   *  @brief extract a set of random numbers with a given probability
   *  distribution
   *  @param [in] nRan the number of random numbers
   *  @param [in] idum the random seed
   *  @param [in] xx vector containing the set of probability
   *  distribution data, x
   *  @param [in] fx vector containing the set of probability
   *  distribution data, f(x)
   *  @param [out] nRandom vector containing the random numbers
   *  @param [in] n_min the minimum number of the output set
   *  @param [in] n_max the maximum number of the output set
   *  @return none
   */
  void random_numbers (const int, const int, const vector<double>, const vector<double>, vector<double> &, const double n_min=par::defaultDouble, const double n_max=1.e30);

  
  /* ======== Alfonso Veropalumbo ======== */

  /**
   *  @brief given a number x, return the closest of two values a, b
   *  @param x the starting value
   *  @param a the first number to test
   *  @param b the second number to test
   *  @return a if x is closer to a, b if x is closer to b
   */
  template <typename T>
    T closest(T x, T a, T b)
    { 
      if (a > b) ErrorMsg("Error in closest, a>b");
      else if (a==b) return a;
      else return (fabs(x-a) < fabs(x-b)) ? a : b;
      return 1;
    }

  /**
   *  @brief given a number x, return the index of the closest element to x in vv
   *  @param x the value
   *  @param vv the vector
   *  @return the index of the closest element to x in vv
   */
  template <typename T>
    T index_closest (T x, vector<T> vv)
    { 
      if (vv.size()==0) ErrorMsg("Error in index_closest, vv is an empty vector");
      vector<double>::iterator low, up;
      low = lower_bound(vv.begin(), vv.end(), x);
      up = upper_bound(vv.begin(), vv.end(), x);
      int index = (closest(x, *low, *up)==*low) ? low-vv.begin() : up-vv.begin();
      return index;
    }

  /**
   *  @brief substitute ~ with the full path
   *  @param path the relative path
   *  @param isDir 1--> directory path, 0->otherwise
   *  @return string containing the full path
   */
  string fullpath (const string path, const bool isDir=1);

  /**
   *  @brief filter W(r/r<SUB>c</SUB>), used e.g. for filtering the
   *  correlation function
   *  @param r the scale at which calculate the filter value
   *  @param rc the characteristic filter scale
   *  @return the filter
   */
  double Filter (const double r, const double rc);

  ///@}


  /* ======== Cosimo Fedeli ======== */

  /// @cond glob
  void gauleg (const double, const double, double *, double *, const int);
  /// @endcond

  ///@}


  // ============================================================================================
  

  /**
   *  @name Functions to manipulate vectors and matrices
   */
  ///@{

  /**
   *  @brief print the elements of a vector on the screen
   *  @param vect a vector
   *  @return none
   */
  template <typename T> 
    void print (const vector<T> vect) 
    {
      for (auto &&i : vect) cout << i << endl;
    }

  /**
   *  @brief print the elements of a two vectors on the screen
   *  @param vect1 a vector
   *  @param vect2 a vector
   *  @return none
   */
  template <typename T> 
    void print (const vector<T> vect1, const vector<T> vect2) 
    {
      if (vect1.size()!=vect2.size()) ErrorMsg("Error in print of Func.h!");
      for (unsigned int i=0; i<vect1.size(); i++) cout <<vect1[i]<<"   "<<vect2[i]<<endl;
    }

  /**
   *  @brief print the elements of a matrix on the screen
   *  @param mat a matrix (i.e. a vector of vectors)
   *  @return none
   */
  template <typename T> 
    void print (const vector<vector<T> > mat) 
    {
      for (unsigned int i=0; i<mat.size(); i++) {
	for (unsigned int j=0; j<mat[i].size(); j++) 
	  cout <<mat[i][j]<<"   ";
	cout <<endl;
      }
    }

  /**
   *  @brief minimum element of a vector
   *  @param vect a vector
   *  @return the minimum element of the vector vect
   */
  template <typename T> 
    T Min (const vector<T> vect) 
    {
      if (vect.size()==0) ErrorMsg("Error in function Min of Func.h: vect.size=0!");
      return *min_element(vect.begin(), vect.end());
    }

  /**
   *  @brief maximum element of a vector
   *  @param vect a vector
   *  @return the maximum element of the vector vect
   */
  template <typename T> 
    T Max (const vector<T> vect) 
    {
      if (vect.size()==0) ErrorMsg("Error in function Max of Func.h: vect.size=0!");
      return *max_element(vect.begin(), vect.end());
    }

  /**
   *  @brief get the unique elements of a vector
   *  @param [in] vect_in a vector
   *  @return vector containing the unique elements of the input
   *  vector vect_in
   */
  template <typename T> 
    vector<T> different_elements (const vector<T> vect_in) 
    { 
      sort(vect_in.begin(), vect_in.end());
      typename vector<T>::iterator it = unique (vect_in.begin(), vect_in.end()); 
      vect_in.resize(it-vect_in.begin());    
      return vect_in;
    }

  /**
   *  @brief get the number of unique elements of a vector
   *  @param vect a vector
   *  @return the number of unique elements of the input vector
   *  vect_in
   */
  template <typename T> 
    int N_different_elements (const vector<T> vect) 
    { 
      sort (vect.begin(),vect.end());
      typename vector<T>::iterator it = unique (vect.begin(), vect.end()); 
      vect.resize(it-vect.begin());    
      return vect.size();
    }

  /**
   *  @brief erase all the equal elements of the input vector
   *  @param [in,out] vv a vector of integer values
   *  @return none
   */
  void unique_unsorted (vector<int> &);

  /**
   *  @brief erase all the equal elements of the input vector
   *  @param [in,out] vv a vector of double values
   *  @return none
   */
  void unique_unsorted (vector<double> &);

  /**
   *  @brief erase all the equal elements of the input vector
   *  @param [in,out] vv a vector of integer values
   *  @return none
   */
  void unique_unsorted (vector<string> &);

  /**
   *  @brief erase some elements of a vector
   *  @param [in,out] vv a vector
   *  @param [in] ind a vector containing the elements of the input
   *  vector vv to be erased
   *  @return none
   */
  template <typename T> 
    void Erase (vector<T> &vv, vector<int> ind) 
    {
      for (auto &&i : ind) 
	if (i>=int(vv.size())) ErrorMsg("Error in Erase of Func.h!");

      unique_unsorted(ind);
      int tt = 0;
      for (auto &&i : ind) 
	vv.erase(vv.begin()+i-(tt++));
    }

  /**
   *  @brief erase some lines of a matrix
   *  @param [in,out] Mat a matrix (i.e. a vector of vectors)
   *  @param [in] ll a vector containing the lines of the input matrix
   *  Mat to be erased
   *  @return none
   */
  template <typename T> 
    void Erase_lines (vector<vector<T> > &Mat, vector<int> ll) 
    {
      for (auto &&i : ll)
	if (i>=int(Mat.size())) ErrorMsg("Error in Erase_lines of Func.h!");

      unique_unsorted(ll);
      int tt = 0;
      for (auto &&i : ll)
	Mat.erase(Mat.begin()+i-(tt++));
    }

  /**
   *  @brief erase some columns of a matrix
   *  @param [in,out] Mat a matrix (i.e. a vector of vectors)
   *  @param [in] col a vector containing the columns of the input
   *  matrix Mat to be erased
   *  @return none
   */
  template <typename T> 
    void Erase_columns (vector<vector<T> > &Mat, vector<int> col) 
    {
      for (auto &&i : col)
	for (auto &&j : Mat)
	  if (i>=int(j.size())) ErrorMsg("Error in Erase_columns of Func.h!");

      unique_unsorted(col);
      int tt = 0;
      for (auto &&i : col) {
	for (auto &&j : Mat)
	  j.erase(j.begin()+i-tt);
	tt ++;
      }
    }
 
  /**
   *  @brief select a submatrix containing lines and columns with all
   *  elements major than \e val
   *
   *  @warning this function works only with rectangular matrices and
   *  only for a small set of specific cases: use it with care and
   *  check carefully the results!
   *
   *  @param [in,out] xx the vector x 
   *  @param [in,out] yy the vector y
   *  @param [in,out] Mat the matrix Mat(x,y) (i.e. a vector of vectors)
   *  @param [in] val a number 
   *  @return none
   */
  template <typename T> 
    void SubMatrix (vector<T> &xx, vector<T> &yy, vector<vector<T> > &Mat, T val) 
    { 
      vector<int> line, column;

      for (unsigned int i=0; i<xx.size(); i++) {
	if (i>=Mat.size()) ErrorMsg("Error in SubMatrix of Func.h!");
	bool ll = 0;

	for (unsigned int j=0; j<yy.size(); j++) {
	  if (j>=Mat[i].size()) ErrorMsg("Error in SubMatrix of Func.h!");
	  if (Mat[i][j]<val) {
	    if (j<int(yy.size()*0.5)) {line.push_back(i); ll = 1;}
	    else if (ll==0) {column.push_back(j);}
	  }
	}
      }
      
      Erase(xx, line);
      Erase_lines(Mat, line);
      
      Erase(yy, column);
      Erase_columns(Mat, column);
    }

  /**
   *  @brief check if the dimensions of two vectors are equal
   *  @param vect1 a vector
   *  @param vect2 a vector
   *  @return 0 &rarr; the dimensions are different; 1 &rarr; the
   *  dimensions are equal
   */
  template <typename T> 
    bool isDimEqual (const vector<T> vect1, const vector<T> vect2) 
    {
      return (vect1.size()==vect2.size()) ? 1 : 0;
    }

  /**
   *  @brief check if the dimensions of two matrices are equal
   *  @param mat1 a matrix
   *  @param mat2 a matrix
   *  @return 0 &rarr; the dimensions are different; 1 &rarr; the
   *  dimensions are equal
   */
  template <typename T> 
    bool isDimEqual (const vector<vector<T> > mat1, const vector<vector<T> > mat2) 
    {
      bool is = (mat1.size()==mat2.size()) ? 1 : 0;
      if (is) 
	for (unsigned int i=0; i<mat1.size(); i++) {
	  if (mat1[i].size()!=mat2[i].size()) is = 0;
	}
      return is;
    }

  /**
   *  @brief check if the dimension of a vector is equal/lower than an
   *  input value
   *  @param vect a vector
   *  @param val the input value
   *  @param vector the name of the vector (used only to write the
   *  error message)
   *  @param equal 1 &rarr; check if the dimension is equal than val;
   *  0 &rarr; check if the dimension is lower than val
   *  @return none
   */
  template <typename T> 
    void checkDim (const vector<T> vect, const int val, const string vector, bool equal=true) 
    {
      if (equal) { 
	if ((int)vect.size()!=val) {
	  string Err = "Error in checkDim of Func.h! The dimension of " + vector + " is:" + conv(vect.size(),par::fINT) + " ( != " + conv(val,par::fINT) + " )";
	  ErrorMsg(Err);
	}
      }
      else {
	if ((int)vect.size()<val) {
	  string Err = "Error in checkDim of Func.h! The dimension of " + vector + " is:" + conv(vect.size(),par::fINT) + " ( >= " + conv(val,par::fINT) + " )";
	  ErrorMsg(Err);
	}
      }
    }
  
  /**
   *  @brief check if the dimensions of a matrix are higher than two
   *  input values
   *  @param mat a matrix
   *  @param i an input value
   *  @param j an input value
   *  @param matrix the name of the matrix (using only to write the
   *  error message)
   *  @return none
   */
  template <typename T> 
    void checkDim (const vector<T> mat, const int i, const int j, const string matrix) 
    {
      if (int(mat.size())<=i) {
	string Err = "Error in checkDim of Func.h! The dimension of: " + matrix + " is:" + conv(mat.size(),par::fINT) + " <= " + conv(i,par::fINT) + "!";
	ErrorMsg(Err);
      } else {
	for (unsigned int k=0; k<mat.size(); k++)
	  if (int(mat[k].size())<=j) {
	    string Err = "Errorin checkDim of Func.h! The dimension of: " + matrix + " is:" + conv(mat[k].size(),par::fINT) + " <= " + conv(j,par::fINT) + "!";
	    ErrorMsg(Err);
	  }
      }
    }

  /**
   *  @brief fill a vector with linearly spaced values
   *  @param [in] nn the number of steps, i.e. the final dimension of
   *  vv
   *  @param [in] min the minimum value of the range of values
   *  @param [in] max the maximum value of the range of values
   *  @return none
   */
  template <typename T> 
    vector<T> linear_bin_vector (const size_t nn, const T min, const T max)
    {
      vector<T> vv(nn);
      for (size_t i = 0; i<nn; i++)
	vv[i] = min+(max-min)*T(i)/T(nn-1);
      return vv;
    }

  /**
   *  @brief fill a vector with logarithmically spaced values
   *  @param [in] nn the number of steps, i.e. the final dimension of
   *  vv
   *  @param [in] min the minimum value of the range of values
   *  @param [in] max the maximum value of the range of values
   *  @return none
   */
  template <typename T> 
    vector<T> logarithmic_bin_vector (const size_t nn, const T min, const T max)
    {
      vector<T> vv(nn);
      for (size_t i=0; i<nn; i++)
	vv[i] = exp(log(min)+(log(max)-log(min))*T(i)/T(nn-1));
      return vv;
    }

  /**
   *  @brief locate a value in a given vector
   *  @author Carlo Giocoli
   *  @author cgiocoli@gmail.com
   *  @param vv a vector of generic values
   *  @param xx a generic number
   *  @return the vector index i such that vv[i]~xx
   */
  template <typename T> 
    int locate (const vector<T> &vv, const T xx) 
    {
      size_t nn = vv.size ();
      int jl = -1;
      int ju = nn;
      bool as = (vv[nn-1] >= vv[0]);
      while (ju-jl > 1)
	{
	  int jm = (ju+jl)*0.5;
	  if ((xx >= vv[jm]) == as)
	    jl = jm;
	  else
	    ju = jm;
	}
      if (xx == vv[0])
	return 0;
      else if (xx == vv[nn-1])
	return nn-2;
      else
	return jl;
    }
  

  // sort two or more vectors at the same time
  namespace glob {
    class CL {
    public:
      vector<double> VV;
      CL(vector<double> vv) {VV = vv;};
    };
    /// @cond glob
    bool operator<(const CL &, const CL &);
    /// @endcond
  }
  
  /**
   *  @brief sort the elements of a vectors, and the elements of a
   *  second vector according to the first sorting
   *
   *  @param p1 iterator to the first vector
   *  @param p2 iterator to the second vector
   *  @param dim dimension of the two vectors 
   *  @return none
   */
  void sort_2vectors (vector<double>::iterator, vector<double>::iterator, const int);

  /**
   *  @brief sort the elements of a vectors, and the elements of two
   *  other vectors according to the first sorting
   *
   *  @param p1 iterator to the first vector
   *  @param p2 iterator to the second vector
   *  @param p3 iterator to the third vector
   *  @param dim dimension of the three vectors 
   *  @return none
   */
  void sort_3vectors (vector<double>::iterator, vector<double>::iterator, vector<double>::iterator, const int);

  /**
   *  @brief sort the elements of a vectors, and the elements of three
   *  other vectors according to the first sorting
   *
   *  @param p1 iterator to the first vector
   *  @param p2 iterator to the second vector
   *  @param p3 iterator to the third vector
   *  @param p4 iterator to the four vector
   *  @param dim dimension of the four vectors 
   *  @return none
   */
  void sort_4vectors (vector<double>::iterator, vector<double>::iterator, vector<double>::iterator, vector<double>::iterator, const int);

  /**
   *  @brief matrix multiplication
   *
   *  overloading of the * operator used to multiplicate two matrices
   *
   *  @param Mat1 STL vector of vectors, i.e. a matrix
   *  @param Mat2 STL vector of vectors, i.e. a matrix
   *  @return Mat1*Mat2
   */
  inline vector<vector<double> > operator * (const vector<vector<double> > &Mat1, const vector<vector<double> > &Mat2)
  {   
    vector<vector<double> > MatP(Mat1.size(), vector<double>(Mat2[0].size(),0.));
  
    for (unsigned int i=0; i<Mat1.size(); i++) 
      for (unsigned int j=0; j<Mat2[0].size(); j++) {
	double temp = 0.;
	for (unsigned int k=0; k<Mat1[0].size(); k++) 
	  temp += Mat1[i][k]*Mat2[k][j];
	MatP[i][j] = temp;
      }

    return MatP;
  }
  
  /**
   *  @brief matrix multiplication
   *
   *  overloading of the * operator used to multiplicate two matrices
   *
   *  @param Mat1 MatDoub matrix of double values of the Numerical libraries
   *  @param Mat2 MatDoub matrix of double values of the Numerical libraries
   *  @return Mat1*Mat2
   */
  inline MatDoub operator * (const MatDoub &Mat1, const MatDoub &Mat2)
    {   
      MatDoub MatP(Mat1.nrows(), Mat2.ncols());
  
      for (int i=0; i<Mat1.nrows(); i++) 
	for (int j=0; j<Mat2.ncols(); j++) {
	  double temp = 0.;
	  for (int k=0; k<Mat1.ncols(); k++) 
	    temp += Mat1[i][k]*Mat2[k][j];
	  MatP[i][j] = temp;
	}
  
      return MatP;
    }

  /**
   *  @brief method to invert a matrix using the GSL
   *  @param [in] mat the matrix to be inverted
   *  @param [out] mat_inv the inverted matrix
   *  @param [in] prec the precision required 
   *  @return none
   */
  void invert_matrix (const vector<vector<double> >, vector<vector<double> > &, const double prec=1.e-10); 

  /**
   *  @brief method to invert a matrix using tge GSL
   *  @param [in] mat the matrix to be inverted
   *  @param [out] mat_inv the inverted matrix
   *  @param [in] i1
   *  @param [in] i2
   *  @param [in] prec the precision required 
   *  @return none
   */
  void invert_matrix (const vector<vector<double> >, vector<vector<double> > &, const int, const int, const double prec=1.e-10); 

  /**
   *  @brief method to invert a 'small' matrix 
   *  @param [in] mat the matrix to be inverted
   *  @param [out] mat_inv the inverted matrix
   *  @param [in] prec the precision required 
   *  @return none
   */
  void invert_small_matrix (const vector<vector<double> >, vector<vector<double> > &, const double prec=1.e-10);

  /**
   *  @brief compute the covariance matrix
   *  @param [in] mat the input matrix
   *  @param [out] cov the output covariance matrix
   *  @param [in] JK 0 &rarr; normalize to 1/(n-1); 1 &rarr; normalize
   *  to n-1/n (for Jackknife)
   *  @return none
   */
  void covariance_matrix (const vector<vector<double> >, vector<vector<double> > &, const bool JK = 0);

  /**
   *  @brief compute the covariance matrix
   *  @param [in] file the vector containing the input files
   *  @param [out] rad the vector containing the binned radii
   *  @param [out] mean the vector containing the mean values
   *  @param [out] cov the output covariance matrix
   *  @return none
   */
  void covariance_matrix (const vector<string>, vector<double> &, vector<double> &, vector<vector<double> > &);


  /* ======== Alfonso Veropalumbo ======== */

  // read and invert the covariance matrix
  void read_cov (const string, vector<vector<double> > &, vector<vector<double> > &, const int, const int);

  // fill a vector from a given distribution 
  void fill_distr (const vector<double>, const double, const vector<double>, vector<double> &, const double, vector<double> &);
  void fill_distr (const vector<double>, const vector<double>, const int, vector<double> &);
  void fill_distr (const int, const vector<double>, const vector<double>, vector<double> &, const double, const double, const int);

  // find the vector index
  void find_index (const vector<double>, const double, const double, int &, int &);

  ///@}


  // ============================================================================================


  /**
   *  @name Functions for statistical analyses
   */
  ///@{

  
  /**
   *  @brief the average of a vector
   *  @param vect the input vector
   *  @return the average of vect
   */
  template <typename T> 
    T Average (const vector<T> vect) 
    {
      T aver = 0;
      if (vect.size()>0) 
	aver = accumulate(vect.begin(),vect.end(),0.)/double(vect.size());
      return aver;
    }

  /**
   *  @brief the weighted average of a vector
   *  @param vect the input vector
   *  @param weight the weight
   *  @return the weighted average of vect
   */
  template <typename T> 
    T Average (const vector<T> vect, const vector<T> weight) 
    {
      if (vect.size()!=weight.size()) ErrorMsg("Error in Average of Func.h");
      T aver = 0;

      vector<T> vw; 
      for (unsigned int i=0; i<vect.size(); i++) vw.push_back(vect[i]*weight[i]);

      if (vect.size()>0) 
	aver = accumulate(vw.begin(),vw.end(),0.)/accumulate(weight.begin(),weight.end(),0.);
      return aver;
    }

  /**
   *  @brief the standard deviation of a vector
   *  @param vect the input vector
   *  @return &sigma;
   */
  template <typename T> 
    T Sigma (const vector<T> vect) 
    {
      T aver = Average(vect);
      T sigma = 0.;
      if (vect.size()>1) {
	for (unsigned int i=0; i<vect.size(); i++)
	  sigma += pow(vect[i]-aver,2);
	sigma = sqrt(sigma/(vect.size()-1));
      }
      return sigma;
    }

  /**
   *  @brief the first, second and third quartiles of a vector
   *  @param [in] vect the input vector
   *  @param [out] first quartile
   *  @param [out] second quartile
   *  @param [out] third quartile
   *  @return none
   */
  template <typename T> 
    void Quartile (vector<T> vect, T & first, T & second, T & third) 
    {
      sort (vect.begin(),vect.end()); 
      vector<T> vect1, vect2;
      int start;
      int n = vect.size();
      if (n>0) {
	if (n==1) {
	  first = -1e10;
	  second = vect[0];
	  third = 1e10;
	}
	if (n>1) {
	  if (n % 2 == 0) { // the number of elemens is even
	    start = int(vect.size()*0.5);
	  } else {
	    start = int(vect.size()*0.5)+1;
	  }
	  for (unsigned int i=0; i<vect.size()*0.5; i++)
	    vect1.push_back(vect[i]);
	  for (unsigned int i=start; i<vect.size(); i++)
	    vect2.push_back(vect[i]);

	  // first quartile
	  n = vect1.size();
	  if (n % 2 == 0) { 
	    first = (vect1[n*0.5-1]+vect1[(n*0.5)])*0.5;
	  } else {
	    first = vect1[(n+1)*0.5-1];
	  }

	  // second quartile = median
	  n = vect.size();
	  if (n % 2 == 0) { 
	    second = (vect[n*0.5-1]+vect[(n*0.5)])*0.5;
	  } else {
	    second = vect[(n+1)*0.5-1];
	  }

	  // third quartile
	  n = vect2.size();
	  if (n % 2 == 0) { 
	    third = (vect2[n*0.5-1]+vect2[(n*0.5)])*0.5;
	  } else {
	    third = vect2[(n+1)*0.5-1];
	  }
	}
      }
    }

  /**
   *  @brief compute the moments of a set of data
   *  @param [in] data the vector containing the input data
   *  @param [out] ave the mean
   *  @param [out] adev the average deviation
   *  @param [out] sdev the standard deviation
   *  @param [out] var the variance
   *  @param [out] skew the skewness
   *  @param [out] curt the kurtosis
   *  @return none
   */
  void Moment (const vector<double>, double &, double &, double &, double &, double &, double &);

  ///@}


  // ============================================================================================
  
  
  /**
   *  @name Special functions
   */
  ///@{

  /**
   *  @brief the quadratic function 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a vector containing the coefficients
   *  @return the quadratic function: par[0]*x<SUP>2</SUP>+par[1]*x+par[2]
   *  @warning pp is not used, it is necessary only for GSL operations
   */
  template <typename T> 
    T Pol2 (T xx, shared_ptr<void> pp, vector<double> par)
    {
      return par[0]*pow(xx,2)+par[1]*xx+par[2];
    }

  /**
   *  @brief the cubic function 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a vector containing the coefficients
   *  @return the cubic function:
   *  par[0]*x<SUP>3</SUP>+par[1]*x<SUP>2</SUP>+par[2]*x+par[3]
   *  @warning pp is not used, it is necessary only for GSL operations
   */
  template <typename T> 
    T Pol3 (T xx, void *pp, vector<double> par) 
    {
      return par[0]*pow(xx,3)+par[1]*pow(xx,2)+par[2]*xx+par[3];
    }

  /**
   *  @brief linear function
   *  @param xx the coordinate x
   *  @return a vector of the Numerical libraries containing
   *  [1,xx]
   */
  inline VecDoub linearfit (const Doub xx) 
  {
    VecDoub vect(2);
    vect[0] = 1.;
    for (Int i=1; i<2; i++) vect[i] = xx*vect[i-1];
    return vect;
  }

  /**
   *  @brief quadratic function
   *  @param xx the coordinate x
   *  @return a vector of the Numerical libraries containing
   *  [1,xx,xx<SUP>2</SUP>]
   */
  inline VecDoub quadratic (const Doub xx) 
  {
    VecDoub vect(3);
    vect[0] = 1.;
    for (Int i=1; i<3; i++) vect[i] = xx*vect[i-1];
    return vect;
  }

  /**
   *  @brief cubic function
   *  @param xx the coordinate x
   *  @return a vector of the Numerical libraries containing
   *  [1,xx,xx<SUP>2</SUP>,xx<SUP>3</SUP>]
   */
  inline VecDoub cubicfit (const Doub xx) 
  {
    VecDoub vect(4);
    vect[0] = 1.;
    for (Int i=1; i<4; i++) vect[i] = xx*vect[i-1];
    return vect;
  }

  /**
   *  @brief the Identity function 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a vector 
   *  @return 1
   *
   *  @warning pp and par are not used, but they are necessary in the
   *  function template
   */
  template <typename T> 
    T identity (T xx, shared_ptr<void> pp, vector<double> par)
    {
      return 1;
    }

  /**
   *  @brief the Gaussian function 
   *  @param xx the variable x
   *  @param pp a void pointer 
   *  @param par a vector containing the coefficients: par[0]=mean,
   *  par[1]=&sigma;
   *  @return the Gaussian function
   *  @warning pp is not used, it is necessary only for GSL operations
   */
  template <typename T> 
    T gaussian (T xx, shared_ptr<void> pp, vector<double> par)
    {
      T gauss = 1./(par[1]*sqrt(2.*par::pi))*exp(-pow(xx-par[0],2)/(2.*par[1]*par[1]));
      return (par.size()==2) ? gauss : gauss*par[2];
    }

  /**
   *  @brief the Maxwellian distribution
   *  @param vel velocity 
   *  @param sigma &sigma;
   *  @return the Maxwellian distribution, P(vel)
   */
  template <typename T> 
    T maxwellian_distr (const T vel, const T sigma) 
    {
      return sqrt(54./par::pi)*pow(vel/sigma,2)*exp(-1.5*pow(vel/sigma,2))/sigma;
    }

  /**
   *  @brief the power-law function
   *  @param xx the variable x
   *  @param x0 the normalization, x<SUB>0</SUB>
   *  @param gamma the slope, &gamma;
   *  @return the power-law function: (x/x<SUB>0</SUB>)<SUP>&gamma;</SUP>
   */
  template <typename T> 
    T powerlaw (const T xx, const T x0, const T gamma)
    {
      return pow(xx/x0,-gamma);
    }

  /**
   *  @brief the double power-law function
   *  @param xx the variable x
   *  @param x0 the normalization, x<SUB>0</SUB>
   *  @param alpha the slope, &alpha;
   *  @param beta the slope, &beta;
   *  @return the double power-law function
   */
  template <typename T> 
    T double_powerlaw (const T xx, const T x0, const T alpha, const T beta) 
    {
      return pow(2.,beta-alpha)/(pow(xx/x0,alpha)*pow(1.+xx/x0,beta-alpha));
    }

  /**
   *  @brief the top-hat window function
   *  @param kR the variable k*R
   *  @return the top-hat window function
   */
  template <typename T> 
    T TopHat_WF (const T kR) 
    {
      return 3.*(sin(kR)-kR*cos(kR))/pow(kR,3);
    }
  
  /**
   *  @brief the radius of a sphere of a given mass and density
   *  @param Mass the mass
   *  @param Rho the density
   *  @return the radius
   */
  template <typename T> 
    T Radius (const T Mass, const T Rho) 
    {
      return pow(3.*Mass/(4.*par::pi*Rho),1./3.); 
    }

  /**
   *  @brief the mass of a sphere of a given radius and density
   *  @param RR the radius
   *  @param Rho the density
   *  @return the mass
   */
  template <typename T> 
    T Mass (const T RR, const T Rho) 
    {
      return 4./3.*par::pi*Rho*pow(RR,3);
    }
  
  /**
   *  @brief the radial velocity
   *  @param vx the velocity along the x direction
   *  @param vy the velocity along the y direction
   *  @param vz the velocity along the z direction
   *  @param ra the Right Ascension
   *  @param dec the Declination
   *  @return the radial velocity
   */
  template <typename T> 
    T radial_velocity (const T vx, const T vy, const T vz, const T ra, const T dec)
    {
      return vx*cos(dec)*sin(ra)+vy*cos(dec)*cos(ra)+vz*sin(dec);
    }

  /**
   *  @brief the Legendre polynomial P<SUB>2</SUB>
   *  @param x the variable x
   *  @return P<SUB>2</SUB>
   */
  template <typename T> 
    T P_2 (const T x) 
    {
      return (3.*x*x-1.)*0.5;
    }
      
  /**
   *  @brief the Legendre polynomial P<SUB>4</SUB>
   *  @param x the variable x
   *  @return P<SUB>4</SUB>
   */
  template <typename T> 
    T P_4 (const T x) 
    {
      return (35.*x*x*x*x-30.*x*x+3.)*0.125;
    }
  
  /**
   *  @brief the Legendre polynomial P<SUB>6</SUB>
   *  @param x the variable x
   *  @return P<SUB>6</SUB>
   */
  template <typename T> 
    T P_6 (const T x) 
    {
      return (231.*x*x*x*x*x*x-315.*x*x*x*x+105.*x*x-5.)*0.0625;
    }

  ///@}
  

  // ============================================================================================


  /**
   *  @name Generic operations on functions 
   */
  ///@{

  /**
   *  @brief measure the var function (where "var" could be the mass,
   *  the luminosity, the radius, etc.)
   *  @param [in] var vector containing the set of data
   *  @param [in] bin the number of bin used
   *  @param [in] V_min the minimum value of the range
   *  @param [in] V_max the maximum value of the range
   *  @param [in] Volume the volume
   *  @param [out] Var vector containing the binned values of "var" 
   *  @param [out] Phi vector containing the binned values of the var
   *  function
   *  @param [out] err vector containing the Poisson errors
   *  @return none
   */
  void measure_var_function (const vector<double>, const int, const double, const double, const double, vector<double> &, vector<double> &, vector<double> &);

  /**
   *  @brief fit a given set of data with a quadratic function 
   *  @param [in] xx vector containing the set of data, x
   *  @param [in] fx vector containing the set of data, f(x)
   *  @param [in] err vector containing the errors, error[f(x)]
   *  @param [out] AA the best-fit parameter A
   *  @param [out] BB the best-fit parameter B
   *  @param [out] CC the best-fit parameter C
   *  @return none
   */ 
  void quad_fit (const vector<double>, const vector<double>, const vector<double>, double &, double &, double &);

  /**
   *  @brief fit a given set of data with a Gaussian function 
   *  @param [in] xx vector containing the set of data, x
   *  @param [in] fx vector containing the set of data, f(x)
   *  @param [out] mean the best-fit value of the mean
   *  @param [out] sigma the best-fit value of &sigma;
   *  @return none
   */
  void gaussian_fit (const vector<double>, const vector<double>, double &, double &);

  /**
   *  @brief convolve a given set of data with a Gaussian function
   *  @param xx vector containing the set of data, x
   *  @param fx vector containing the set of data, f(x)
   *  @param mean the mean 
   *  @param sigma &sigma;
   *  @return the Gaussian convolution
   */
  double gaussian_convolution (const vector<double>, const vector<double>, const double, const double);

  /**
   *  @brief derive and store the number distribution of a given
   *  vector 
   *  @param [out] xx vector containing the binned values of the
   *  variable 
   *  @param [out] fx vector containing the binned values of the
   *  distribution
   *  @param [in] FF vector containing the given set of data
   *  @param [in] WW vector containing the weights
   *  @param [in] nbin the number of bins
   *  @param [in] linear 1 &rarr; linear binning; 0 &rarr; logarithmic binning 
   *  @param [in] file_out the output file where the distribution is
   *  stored
   *  @param [in] fact factor used to normalized the distribution
   *  @param [in] V1 the minimum limit of the distribution
   *  @param [in] V2 the maximum limit of the distribution
   *  @param [in] bin_type 1 &rarr; dn/dvar; 0 &rarr; dn/dlogvar
   *  @param [in] conv 1 &rarr; compute the Gaussian convolvolution of the distribution; 0 &rarr; do not convolve
   *  @param [in] sigma &sigma; of the Gaussian kernel
   *  @return none
   */
  void distribution (vector<double> &, vector<double> &, const vector<double>, const vector<double>, const int, const bool linear=1, const string file_out=par::defaultString, const double fact=1., const double V1=par::defaultDouble, const double V2=par::defaultDouble, const bool bin_type=1, const bool conv=0, const double sigma=0);

  /**
   *  @brief simple Monte Carlo integration of f(x)
   *  @param func the function f(x)
   *  @param x1 minimum limit of the integral
   *  @param x2 maximum limit of the integral
   *  @return \f$\int_{x1}^{x2} f(x)dx\f$
   */
  double MC_Int (double func(const double), const double x1, const double x2); 

  /**
   *  @brief simple Monte Carlo integration of f(x,A)
   *  @param func the function f(x,A)
   *  @param AA parameter of the function, A
   *  @param x1 minimum limit of the integral
   *  @param x2 maximum limit of the integral
   *  @return \f$\int_{x1}^{x2} f(x)dx\f$
   */
  double MC_Int (double func(const double, const double AA), const double AA, const double x1, double x2); 
  
  /**
   *  @brief simple Monte Carlo integration of f(x,A,B,C,D,E)
   *  @param func the function f(x,A,B,C,D,E)
   *  @param AA parameter of the function, A
   *  @param BB parameter of the function, B
   *  @param CC parameter of the function, C
   *  @param DD parameter of the function, D
   *  @param EE parameter of the function, E
   *  @param x1 minimum limit of the integral
   *  @param x2 maximum limit of the integral
   *  @return \f$\int_{x1}^{x2} f(x,A,B,C,D,E)dx\f$
   */
  double MC_Int (double func(const double, const double AA, const double BB, const double CC, const double DD, const double EE), const double AA, const double BB, const double CC, const double DD, const double EE, const double x1, const double x2); 

  // ============================================================================================
 
  /**
   *  @brief compute the first derivative of a given function
   *  @param XX the point where to compute the derivative
   *  @param xx vector containing the set of data, x
   *  @param yy vector containing the function, y(x)
   *
   *  @param interpType type the method used to interpolate or
   *  extrapolate: "Linear" &rarr; linear interpolation; "Poly" &rarr;
   *  polynomial interpolation; "Spline" &rarr; cubic spline
   *  interpolation; "Rat" &rarr; diagonal rational function
   *  interpolation; "BaryRat" &rarr; barycentric rational
   *  interpolation
   *
   *  @param Num number of points used in the polynomial, diagonal
   *  rational and barycentric rational interpolations
   *
   *  @param stepsize the binning size
   *  @return the first derivative of the function y(x)
   */
  double D1 (const double, const vector<double>, const vector<double>, const string, const int, const double);

  /**
   *  @brief compute the second derivative of a given function
   *  @param XX the point where to compute the derivative
   *  @param xx vector containing the set of data, x
   *  @param yy vector containing the function, y(x)
   *
   *  @param interpType type the method used to interpolate or
   *  extrapolate: "Linear" &rarr; linear interpolation; "Poly" &rarr;
   *  polynomial interpolation; "Spline" &rarr; cubic spline
   *  interpolation; "Rat" &rarr; diagonal rational function
   *  interpolation; "BaryRat" &rarr; barycentric rational
   *  interpolation
   *
   *  @param Num number of points used in the polynomial, diagonal
   *  rational and barycentric rational interpolations
   *
   *  @param stepsize the binning size
   *  @return the second derivative of the function y(x)
   */
  double D2 (const double, const vector<double>, const vector<double>, const string, const int, const double);

  /**
   *  @brief compute the n-th order derivative of a given function
   *  @param nd the order of the derivative
   *  @param XX the point where to compute the derivative
   *  @param xx vector containing the set of data, x
   *  @param yy vector containing the function, y(x)
   *
   *  @param interpType type the method used to interpolate or
   *  extrapolate: "Linear" &rarr; linear interpolation; "Poly" &rarr;
   *  polynomial interpolation; "Spline" &rarr; cubic spline
   *  interpolation; "Rat" &rarr; diagonal rational function
   *  interpolation; "BaryRat" &rarr; barycentric rational
   *  interpolation
   *
   *  @param Num number of points used in the polynomial, diagonal
   *  rational and barycentric rational interpolations
   *
   *  @param stepsize the binning size
   *  @return the nd derivative of the function y(x)
   */
  double Deriv (const int, const double, const vector<double>, const vector<double>, const string, const int Num=-1, const double stepsize=1.);

  /**
   *  @brief compute the n-th order derivative of a given function
   *  @param nd the order of the derivative
   *  @param XX the point where to compute the derivative
   *  @param Func the function to be derived
   *
   *  @param interpType type the method used to interpolate or
   *  extrapolate: "Linear" &rarr; linear interpolation; "Poly" &rarr;
   *  polynomial interpolation; "Spline" &rarr; cubic spline
   *  interpolation; "Rat" &rarr; diagonal rational function
   *  interpolation; "BaryRat" &rarr; barycentric rational
   *  interpolation
   *
   *  @param Num number of points used in the polynomial, diagonal
   *  rational and barycentric rational interpolations
   *
   *  @param stepsize the binning size
   *  @return the nd derivative of the function y(x)
   */
  template<typename T> double Deriv(const int nd, const double XX, const T &Func, const string interpType, const int Num=-1, const double stepsize=1.) {
    double err = -1.;
    double DD = dfridr(Func, XX, stepsize, err);
      
    if (nd==1) DD = dfridr(Func, XX, stepsize, err);

    else {
    
      double xmin = XX*0.1; // check!!!
      double xmax = XX*10.; // check!!!
      int step = 100; // check!!!
      vector<double> xg, yg;
    
      for (int i=0; i<step; i++)
	xg.push_back(xmin+(xmax-xmin)*double(i)/double(step-1));

      for (unsigned int i=0; i<xg.size(); i++) 
	yg.push_back(Func(xg[i]));

      DD = Deriv(nd,XX,xg,yg,interpType,Num,stepsize);
    
    } 
  
    if (err>fabs(DD)*1.e-1) { 
      double errR = err/fabs(DD)*100.;
      string Warn = "Attention: the error in the derivative is = " + conv(errR,par::fDP3) + " % !"; 
      WarningMsg(Warn);
    }

    return DD;
  }

  /**
   *  @brief create a 1D grid given an input function
   *  @param file_grid the file with the input function
   *  @param [in] func the input function
   *  @param [in] par the function parameters
   *  @param [in] bin the number of bins in the grid
   *  @param [in] x_min the minimum limit of the grid
   *  @param [in] x_max the maximum limit of the grid
   *  @param [in] binning the binning type: it can be "lin", "loglin" or
   *  "log"
   *  @param [in,out] xx vector containing the grid points
   *  @param [in,out] yy vector containing the values of the function at the
   *  grid points
   *  @return none
   */
  void bin_function (const string, double func(double, void*), void *, const int, const double, const double, const string, vector<double> &, vector<double> &);

  /**
   *  @brief create a 2D grid given an input function
   *  @param [in] file_grid the file with the input function
   *  @param [in] func the input function
   *  @param [in] par the function parameters
   *  @param [in] bin the number of bins in the grid
   *  @param [in] x1_min the minimum limit of the grid in one direction
   *  @param [in] x1_max the maximum limit of the grid in one direction
   *  @param [in] x2_min the minimum limit of the grid in one direction
   *  @param [in] x2_max the maximum limit of the grid in one direction
   *  @param [in] binning the binning type: it can be "lin", "loglin" or
   *  "log"
   *  @param [in,out] xx1 vector containing the grid points in one direction
   *  @param [in,out] xx2 vector containing the grid points in one direction
   *  @param [in,out] yy vector containing the values of the function at the
   *  grid points
   *  @return none
   */
  void bin_function_2D (const string, double func(double *, size_t, void *), void *, const int, const double, const double, const double, const double, const string, vector<double> &, vector<double> &, vector<vector<double> > &);

  /// @cond glob
  double func_grid_lin (double, void *);
  double func_grid_loglin (double, void *);
  double func_grid_log (double, void *);
  double func_grid_lin_2D (double *, size_t, void *);
  double func_grid_loglin_2D (double *, size_t, void *);
  double func_grid_log_2D (double *, size_t, void *);
  /// @endcond

  
  /* ======== Alfonso Veropalumbo ======== */

  // convolution of two function
  void convolution (vector<double>, vector<double>, vector<double> &, double);

  ///@}
  

  // ============================================================================================


  /**
   *  @name Functions to calculate distances (angles must be in radians) 
   */
  ///@{

  /**
   *  @brief conversion from Cartesian coordinates to polar
   *  coordinates
   *
   *  @param [in] XX the Cartesian coordinate x
   *  @param [in] YY the Cartesian coordinate y
   *  @param [in] ZZ the Cartesian coordinate z
   *  @param [out] ra the Right Ascension
   *  @param [out] dec the Declination
   *  @param [out] dd the comoving distance
   *  @return none
   */
  void polar_coord (const double, const double, const double, double &, double &, double &); 

  /**
   *  @brief conversion from polar coordinates to Cartesian
   *  coordinates
   *
   *  @param [in] ra the Right Ascension
   *  @param [in] dec the Declination
   *  @param [in] dd the comoving distance
   *  @param [out] XX the Cartesian coordinate x
   *  @param [out] YY the Cartesian coordinate y
   *  @param [out] ZZ the Cartesian coordinate z
   *  @return none
   */
  void cartesian_coord (const double, const double, const double, double &, double &, double &);

  /**
   *  @brief conversion from Cartesian coordinates to polar
   *  coordinates used for a set of objects
   *
   *  @param [in] XX vector containing the Cartesian coordinates x
   *  @param [in] YY vector containing the Cartesian coordinates y
   *  @param [in] ZZ vector containing the Cartesian coordinates z
   *  @param [out] ra vector containing the Right Ascension values
   *  @param [out] dec vector containing the Declination values
   *  @param [out] dd vector containing the comoving distances
   *  @return none
   */
  void polar_coord (const vector<double>, const vector<double>, const vector<double>, vector<double> &, vector<double> &, vector<double> &); 

  /**
   *  @brief conversion from polar coordinates to Cartesian
   *  coordinates used for a set of objects
   *
   *  @param [in] ra vector containing the Right Ascension values
   *  @param [in] dec vector containing the Declination values
   *  @param [in] dd vector containing the comoving distances
   *  @param [out] XX vector containing the Cartesian coordinates x
   *  @param [out] YY vector containing the Cartesian coordinates y
   *  @param [out] ZZ vector containing the Cartesian coordinates z
   *  @return none
   */
  void cartesian_coord (const vector<double>, const vector<double>, const vector<double>, vector<double> &, vector<double> &, vector<double> &);

  /**
   *  @brief the Euclidean distance in 3D relative to the origin
   *  (0,0,0), i.e. the Euclidean norm
   *  @param x1 x coordinate of the first object
   *  @param x2 x coordinate of the second object
   *  @param y1 y coordinate of the first object
   *  @param y2 y coordinate of the second object
   *  @param z1 z coordinate of the first object
   *  @param z2 z coordinate of the second object
   *  @return the Euclidean distance
   */
  double Euclidean_distance (const double x1, const double x2, const double y1, const double y2, const double z1, const double z2);

  /**
   *  @brief the perpendicular separation, r<SUB>p</SUB>
   *  @param ra1 the Right Ascension of the first object
   *  @param ra2 the Right Ascension of the second object
   *  @param dec1 the Declination of the first object
   *  @param dec2 the Declination of the second object
   *  @param d1 the comoving distance of the first object
   *  @param d2 the comoving distance of the second object
   *  @return the perpendicular separation, r<SUB>p</SUB>
   */
  double perpendicular_distance (const double ra1, const double ra2, const double dec1, const double dec2, const double d1, const double d2);
  
  /**
   *  @brief the angular separation in 3D 
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   *  @param x1 x coordinate of the first object
   *  @param x2 x coordinate of the second object
   *  @param y1 y coordinate of the first object
   *  @param y2 y coordinate of the second object
   *  @param z1 z coordinate of the first object
   *  @param z2 z coordinate of the second object
   *  @return the angular separation
   */
  double angular_distance (const double x1, const double x2, const double y1, const double y2, const double z1, const double z2);

  /**
   *  @brief the haversine angular separation in 3D
   *  @author Alfonso Veropalumbo
   *  @author alfonso.veropalumbo@unibo.it
   *  @param ra1 the Right Ascension of the first object
   *  @param ra2 the Right Ascension of the second object
   *  @param dec1 the Declination of the first object
   *  @param dec2 the Declination of the second object
   *  @return the haversine angular separation
   */
  double haversine_distance (const double ra1, const double ra2, const double dec1, const double dec2);
  

  /* ======== Alfonso Veropalumbo ======== */

  // sdss<->eq transformation
  void eq2sdss (const vector<double>, const vector<double>, vector<double> &, vector<double> &); 
  void sdss2eq (const vector<double>, const vector<double>, vector<double> &, vector<double> &);
  void sdss_stripe (const vector<double>, const vector<double>, vector<int> &, vector<int> &);


  ///@}


  // ============================================================================================
  
  
  /**
   *  @name Functions to model the correlation function
   */
  ///@{
  
 
  /**
   *  @brief the two-point correlation function computed from the
   *  Fourier transform of the power spectrum
   * 
   *  @param rr the comoving separation, r
   *
   *  @param lgkk vector containing the logarithm of the wave vectors,
   *  log<SUB>10</SUB>k
   *
   *  @param lgPk vector containing the logarithm of the power
   *  spectrum, log<SUB>10</SUB>P(k)
   *
   *  @param k_min the minimum value of the wave vector used in the
   *  integral of the Fourier transform
   *
   *  @param k_max the maximum value of the wave vector used in the
   *  integral of the Fourier transform
   *
   *  @param aa parameter used to smooth the integrand, given by the
   *  eq. 24 of Anderson et al. 2012
   *
   *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
   *  the GSL libraries are used
   *
   *  @param prec accuracy of the GSL integration 
   *
   *  @return the two-point correlation function, &xi;(r)
   */
  double xi_from_Pk (const double, const vector<double>, const vector<double>, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2);

  /**
   *  @brief the two-point correlation function computed from the
   *  Fourier transform of the power spectrum read from a file
   * 
   *  @param rr the comoving separation, r
   *
   *  @param file name of the file where the power spectrum is stored
   *
   *  @param c1 the column of the file corresponding to the wave
   *  vector, k
   *
   *  @param c2 the column of the file corresponding to the power
   *  spectrum, P(k)
   *
   *  @param k_min the minimum value of the wave vector used in the
   *  integral of the Fourier transform
   *
   *  @param k_max the maximum value of the wave vector used in the
   *  integral of the Fourier transform
   *
   *  @param aa parameter used to smooth the integrand, given by the
   *  eq. 24 of Anderson et al. 2012
   *
   *  @param GSL 0 &rarr; the Numerical libraries are used; 1 &rarr;
   *  the GSL libraries are used
   *
   *  @param prec accuracy of the GSL integration 
   *
   *  @return the two-point correlation function, &xi;(r)
   */
  double xi_from_Pk (const double, const string, const int c1=1, const int c2=2, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2);

  /**
   *  @brief the power spectrum computed from the Fourier transform of
   *  the two-point correlation function
   * 
   *  @param kk the wave vector, k
   *
   *  @param lgrr vector containing the logarithm of the comoving
   *  separations, log<SUB>10</SUB>r
   *
   *  @param lgxi vector containing the logarithm of the two-point
   *  correlation function, log<SUB>10</SUB>&xi;(r)
   *
   *  @param r_min the minimum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @return the power spectrum, P(k)
   */
  double Pk_from_xi (const double, const vector<double>, const vector<double>, const double r_min=0.03, const double r_max=100.); 

  /**
   *  @brief the power spectrum computed from the Fourier transform of
   *  the two-point correlation function read from a file
   * 
   *  @param kk the wave vector, k
   *
   *  @param file name of the file where the two-point correlation
   *  function is stored
   *
   *  @param c1 the column of the file corresponding to the comoving
   *  separation, r
   *
   *  @param c2 the column of the file corresponding to the two-point
   *  correlation function, &xi;(r)
   *
   *  @param r_min the minimum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral of the Fourier transform
   *
   *  @return the power spectrum, P(k)
   */
  double Pk_from_xi (const double, const string, const int c1=1, const int c2=2, const double r_min=0.03, const double r_max=100.); 

  namespace glob {
    /// @cond glob
    double func_xi_GSL (double, void *);
    double func_SSM_GSL (double, void *);
    /// @endcond
  }

  /**
   *  @brief the projected two-point correlation function
   * 
   *  this function estimates the projected correlation function by
   *  integrating a given two-point correlation function as follows:
   *  \f[
   *  w_p(r_p)=\int_{r_p}^{r_{max}}\frac{\xi(r)}{\sqrt{r^2-r_p^2}}rdr
   *  \f]
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param rr vector containing the central values of the binned
   *  comoving separations, r
   *
   *  @param xi vector containing the central values of the binned
   *  two-point correlation function, &xi;(r)
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral 
   *
   *  @return the projected correlation function, w(r<SUB>p</SUB>)
   */
  double wp (const double, const vector<double>, const vector<double>, const double r_max=100.); 

  /**
   *  @brief the projected two-point correlation function
   * 
   *  this function estimates the projected correlation function by
   *  integrating a given two-point correlation function as follows:
   *  \f[
   *  w_p(r_p)=\int_{r_p}^{r_{max}}\frac{\xi(r)}{\sqrt{r^2-r_p^2}}rdr
   *  \f]
   *
   *  the two-point correlation function is read from a file
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param file name of the file where the two-point correlation
   *  function is stored
   *
   *  @param r_max the maximum value of the comoving separation used
   *  in the integral 
   *
   *  @return the projected correlation function,
   *  w<SUB>p</SUB>(r<SUB>p</SUB>)
   */
  double wp (const double, const string, const double r_max=100.); 

  /**
   *  @brief the rms mass fluctuation within a given radius
   * 
   *  computed with the spherically averaged correlation function as
   *  follows:
   *
   *  \f[
   *  \sigma_R^2=\int_0^2dy\,y^2\xi(yR)\left(3-\frac{9y}{4}+\frac{3y^3}{16}\right)
   *  \f]
   *
   *  or with the projected correlation function as follows:
   *
   *  \f[
   *  \sigma_R^2=\frac{1}{R^3}\int_0^\infty dr_p\,r_pw_p(r_p)g(r_p/R)
   *  \f]
   *
   *  where g(x) is:
   *
   *  \f[
   *  \left\{
   *  \begin{array}{ll}
   *  \frac{1}{2\pi}\left(3\pi-9x+x^3\right) & \mbox{for}\,x\le2 \\
   *  \frac{1}{2\pi}\left(\frac{-x^4+11x^2-28}{\sqrt{x^2-4}}+x^3-9x+6\sin^{-1}(2/x)\right) & \mbox{for}\,x>2 \\
   *  \end{array}
   *  \right.
   *  \f] 
   *  see e.g. Zehavi et al. 2005, ApJ, 621, 22
   *
   *  @param RR the radius R [Mpc/h]
   *
   *  @param corrType 1 &rarr; the spherically averaged correlation
   *  function is used; 2 &rarr; the projected correlation function is
   *  used
   *
   *  @param rr vector containing comoving separations
   *
   *  @param corr vector containing the two-point correlation function
   *
   *  @return &sigma;<SUB>R</SUB>: the rms mass fluctuation within a radius R [Mpc/h] 
   */
  double sigmaR (const double, const int, const vector<double>, const vector<double>);

  /**
   *  @brief the projected correlation function,
   *  w<SUB>p</SUB>(r<SUB>p</SUB>), computed by modelling the
   *  two-point correlation function, &xi;(r), as a power-law
   * 
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param r0 r<SUB>0</SUB>: the clustering normalization
   *
   *  @param gamma &gamma;: the clustering slope
   *
   *  @return the projected correlation function
   *  w<SUB>p</SUB>(r<SUB>p</SUB>)
   */
  double xi_projected_powerlaw (const double, const double, const double); 

  /**
   *  @brief the ratio between the redshift-space and real-space
   *  correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:
   * 
   *  \f[ \frac{\xi(s)}{\xi(r)} = 1 + \frac{2\beta}{3} +
   *  \frac{\beta^2}{5} \f]
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @return &xi;(s)/&xi;(r)
   */
  double xi_ratio (const double);                         

  /**
   *  @brief the ratio between the redshift-space and real-space
   *  correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:
   * 
   *  \f[ \frac{\xi(s)}{\xi(r)} = 1 + \frac{2\beta}{3} +
   *  \frac{\beta^2}{5} \f]
   *
   *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
   *  
   *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
   *
   *  @return &xi;(s)/&xi;(r)
   */
  double xi_ratio (const double, const double);                

  /// @cond glob
  /**
   *  @brief the ratio between the redshift-space and real-space
   *  correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:

   *  \f[ \frac{\xi(s)}{\xi(r)} = 1 + \frac{2\beta}{3} +
   *  \frac{\beta^2}{5} \f]
   *
   *  @param xx coordinate x
   *  @param pp void pointer
   *  @param par vector containing one or two parameters
   *  @return &xi;(s)/&xi;(r)
   */
  double xi_ratio (double, shared_ptr<void> , vector<double>); // for the chi2 function
  /// @endcond

  /**
   *  @brief error on the ratio between the redshift-space and
   *  real-space correlation functions
   *
   *  as predicted by the large-scale limit of the Kaiser/Hamilton
   *  model:
   * 
   *  \f[ \delta\left[\frac{\xi(s)}{\xi(r)}\right] = \left(\frac{2}{3}
   *  + \frac{2\beta}{5}\right)\cdot\delta(\beta) \f]
   *
   *  @param beta &beta;=f/b
   *  
   *  @param error_beta error on &beta;
   *
   *  @return &delta;[&xi;(s)/&xi;(r)]
   */
  double error_xi_ratio (const double, const double); 

  /**
   *  @brief the barred correlation function
   *
   *  computed with a simple rectangular integration 
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\xi}(r) = \frac{3}{r^3}\int^r_0dr'\xi(r')r'{^2}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr vector containing the input comoving separations
   *
   *  @param xi vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\xi}(r) \f$
   */
  double barred_xi_direct (const double, const vector<double>, const vector<double>, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief the double barred correlation function
   *
   *  computed with a simple rectangular integration 
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\overline{\xi}}(r) = \frac{5}{r^5}\int^r_0dr'\xi(r')r'{^4}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr vector containing the input comoving separations
   *
   *  @param xi vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\overline{\xi}}(r) \f$
   */
  double barred_xi__direct (const double, const vector<double>, const vector<double>, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief the barred correlation function
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\xi}(r) = \frac{3}{r^3}\int^r_0dr'\xi(r')r'{^2}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr vector containing the input comoving separations
   *
   *  @param xi vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\xi}(r) \f$
   */
  double barred_xi_ (const double, const vector<double>, const vector<double>, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief the double barred correlation function
   *
   *  (see e.g. Hamilton 1992)
   *
   *  \f[ \overline{\overline{\xi}}(r) = \frac{5}{r^5}\int^r_0dr'\xi(r')r'{^4}
   *  \f]
   *
   *  @param RR the comoving separation, r, where the barred
   *  correlation function is computed
   *
   *  @param rr vector containing the input comoving separations
   *
   *  @param xi vector containing the input two-point correlation
   *  function
   *
   *  @param rAPP comoving scale below which a power-law model for the
   *  two-point correlation function is assumed in the integral
   *
   *  @param r0 the power-law normalization, r<SUB>0</SUB>  
   *
   *  @param gamma the power-law slope, &gamma;
   *
   *  @return \f$ \overline{\overline{\xi}}(r) \f$
   */
  double barred_xi__ (const double, const vector<double>, const vector<double>, const double rAPP=0., const double r0=-1., const double gamma=1.); 

  /**
   *  @brief xi<SUB>0</SUB>(s) from &xi;(r,&mu;)
   *
   *  \f[ \xi_0(s) = \frac{1}{2}\int_{-1}^1\xi(s,\mu)d\mu \f]
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu vector containing the angle between the separation
   *  vector and the line of sight
   *
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *
   *  @return xi<SUB>0</SUB>(s)
   */
  double multipole_xi0 (const int, const vector<double>, const vector<vector<double> >);
  
  /**
   *  @brief xi<SUB>2</SUB>(s) from &xi;(r,&mu;)
   *
   *  \f[ \xi_2(s) = \frac{5}{2}\int_{-1}^1\xi(s,\mu)L_2(\mu)d\mu \f]
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu vector containing the angle between the separation
   *  vector and the line of sight
   *
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *
   *  @return xi<SUB>2</SUB>(s)
   */
  double multipole_xi2 (const int, const vector<double>, const vector<vector<double> >);
 
  /**
   *  @brief xi<SUB>4</SUB>(s) from &xi;(r,&mu;)
   *
   *  \f[ \xi_4(s) = \frac{9}{2}\int_{-1}^1\xi(s,\mu)L_4(\mu)d\mu \f]
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu vector containing the angle between the separation
   *  vector and the line of sight
   *
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *
   *  @return xi<SUB>4</SUB>(s)
   */
  double multipole_xi4 (const int, const vector<double>, const vector<vector<double> >);
  
  /**
   *  @brief error on xi<SUB>0</SUB>(s) from &xi;(r,&mu;)
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu vector containing the angle between the separation
   *  vector and the line of sight
   *
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *
   *  @return error on xi<SUB>0</SUB>(s)
   */
  double error_multipole_xi0 (const int, const vector<double>, const vector<vector<double> >);

  /**
   *  @brief error on xi<SUB>2</SUB>(s) from &xi;(r,&mu;)
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu vector containing the angle between the separation
   *  vector and the line of sight
   *
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *
   *  @return error on xi<SUB>2</SUB>(s)
   */
  double error_multipole_xi2 (const int, const vector<double>, const vector<vector<double> >);

  /**
   *  @brief error on xi<SUB>4</SUB>(s) from &xi;(r,&mu;)
   *
   *  @param indexR index correspondent to the comoving separation
   *  where the multipole is computed
   *
   *  @param mu vector containing the angle between the separation
   *  vector and the line of sight
   *
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *
   *  @return error on xi<SUB>4</SUB>(s)
   */
  double error_multipole_xi4 (const int, const vector<double>, const vector<vector<double> >);

  /**
   *  @brief xi<SUB>0</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  \f[ \xi_0(s) = \frac{1}{2}\int_{-1}^1\xi(s,\mu)d\mu \f]
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp vector containing the values of r<SUB>p</SUB>
   *  @param pi vector containing the values of &pi;
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return xi<SUB>0</SUB>(s)
   */
  double multipole_xi0 (const double, const vector<double>, const vector<double>, const vector<vector<double> >, const double);

  /**
   *  @brief xi<SUB>2</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  \f[ \xi_2(s) = \frac{5}{2}\int_{-1}^1\xi(s,\mu)L_2(\mu)d\mu \f]
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp vector containing the values of r<SUB>p</SUB>
   *  @param pi vector containing the values of &pi;
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return xi<SUB>2</SUB>(s)
   */
  double multipole_xi2 (const double, const vector<double>, const vector<double>, const vector<vector<double> >, const double);

  /**
   *  @brief xi<SUB>4</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  \f[ \xi_4(s) = \frac{9}{2}\int_{-1}^1\xi(s,\mu)L_4(\mu)d\mu \f]
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp vector containing the values of r<SUB>p</SUB>
   *  @param pi vector containing the values of &pi;
   *  @param xi matrix containing the values of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return xi<SUB>4</SUB>(s)
   */
  double multipole_xi4 (const double, const vector<double>, const vector<double>, const vector<vector<double> >, const double);

  /**
   *  @brief error on xi<SUB>0</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp vector containing the values of r<SUB>p</SUB>
   *  @param pi vector containing the values of &pi;
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return error on xi<SUB>0</SUB>(s)
   */
  double error_multipole_xi0 (const double, const vector<double>, const vector<double>, const vector<vector<double> >, const double);

  /**
   *  @brief error on xi<SUB>2</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp vector containing the values of r<SUB>p</SUB>
   *  @param pi vector containing the values of &pi;
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return error on xi<SUB>2</SUB>(s)
   */
  double error_multipole_xi2 (const double, const vector<double>, const vector<double>, const vector<vector<double> >, const double);

  /**
   *  @brief error on xi<SUB>4</SUB>(s) from &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param ss comoving scale where the multipole is computed
   *  @param rp vector containing the values of r<SUB>p</SUB>
   *  @param pi vector containing the values of &pi;
   *  @param error matrix containing the errors of &xi;(r,&mu;)
   *  @param delta_s bin size 
   *  @return error on xi<SUB>4</SUB>(s)
   */
  double error_multipole_xi4 (const double, const vector<double>, const vector<double>, const vector<vector<double> >, const double);

  /// @cond glob
  /**
   *  @brief multipoles (&xi;<SUB>0</SUB> + &xi;<SUB>2</SUB>) of the
   *  two-point correlation function used in the &chi;<SUP>2</SUP>
   *
   *  @param rr the comoving separation
   *  @param pp a void pointer 
   *  @param par a vector containing the coefficients
   *
   *  @return multipoles (&xi;<SUB>0</SUB> + &xi;<SUB>2</SUB>) of the
   *  two-point correlation
   */
  double multipoles (double, shared_ptr<void> , vector<double>);
  /// @endcond

  /**
   *  @brief the model multipole &xi;<SUB>0</SUB> of the two-point
   *  correlation function
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias 
   *
   *  @param xi_real the real-space two-point correlation function
   *
   *  @return the multipole &xi;<SUB>0</SUB>
   */
  double multipole_xi0_model (const double, const double);

  /**
   *  @brief the model multipole &xi;<SUB>0</SUB> of the two-point
   *  correlation function
   *
   *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
   *  
   *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
   *
   *  @param sigma8z &sigma;<SUB>8</SUB>
   *  
   *  @param xi_DM the real-space two-point correlation function of
   *  the dark matter
   *
   *  @return the multipole &xi;<SUB>0</SUB>
   */
  double multipole_xi0_model (const double, const double, const double, const double);

  /// @cond glob
  /**
   *  @brief the model multipole &xi;<SUB>0</SUB> of the two-point
   *  correlation function
   *
   *  @param xx coordinate x
   *  @param pp void pointer
   *  @param par vector containing one or two parameters
   *
   *  @return the multipole &xi;<SUB>0</SUB>
   */
  double multipole_xi0_model (double, shared_ptr<void> , vector<double>);
  /// @endcond

  /**
   *  @brief the model multipole &xi;<SUB>2</SUB> of the two-point
   *  correlation function
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias 
   *
   *  @param xi_real the real-space two-point correlation function
   *
   *  @param xi_ \f$ \overline{\xi}(r) \f$
   *
   *  @return the multipole &xi;<SUB>2</SUB>
   */
  double multipole_xi2_model (const double, const double, const double); 

  /**
   *  @brief the model multipole &xi;<SUB>4</SUB> of the two-point
   *  correlation function
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias 
   *
   *  @param xi_real the real-space two-point correlation function
   *
   *  @param xi_ \f$ \overline{\xi}(r) \f$
   *
   *  @param xi__ \f$ \overline{\overline{\xi}}(r) \f$
   *
   *  @return the multipole &xi;<SUB>4</SUB>
   */
  double multipole_xi4_model (const double, const double, const double, const double);

  /// @cond glob
  // theoretical model for the linear xi(rp,pi)
  double xi2D_lin_model (double, double, shared_ptr<void> , vector<double>);
  /// @endcond

  /**
   *  @brief the linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias
   *
   *  @param xi_real the real-space correlation function
   *
   *  @param xi_ \f$ \overline{\xi}(r) \f$
   *
   *  @param xi__ \f$ \overline{\overline{\xi}}(r) \f$
   *
   *  @param P_2 the Legendre polynomial P<SUB>2</SUB>
   *
   *  @param P_4 the Legendre polynomial P<SUB>4</SUB>
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_lin_model (const double, const double, const double, const double, const double, const double, const double);

  /**
   *  @brief the linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias  
   *
   *  @param rad_real vector containing the binnend values of the
   *  comoving separations
   *
   *  @param xi_real vector containing the binnend values of the
   *  real-space correlation function 
   *
   *  @param xi_ vector containing the binnend values of \f$
   *  \overline{\xi}(r) \f$
   *
   *  @param xi__ vector containing the binnend values of \f$
   *  \overline{\overline{\xi}}(r) \f$
   *
   *  @param index index for internal use
   *
   *  @param bias_nl 0 &rArr; linear bias; &rArr; 1 non-linear bias 
   *
   *  @param bA the parameter b<SUB>A</SUB> used to model the bias
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_lin_model (const double, const double, const double, const double, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const int index=-1, const bool bias_nl=0, const double bA=0.);

  /// @cond glob
  // dispersion model for xi(rp,pi)
  double xi2D_model (double, double, shared_ptr<void> , vector<double>);
  /// @endcond

  /**
   *  @brief the non-linear dispersion model for
   *  &xi;(r<SUB>p</SUB>,&pi;)
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param beta &beta;=f/b, where f is the linear growth rate and b
   *  is the bias
   *
   *  @param bias the bias
   *
   *  @param sigma12 &sigma;<SUB>12</SUB>
   *
   *  @param rad_real vector containing the binnend values of the
   *  comoving separations
   *
   *  @param xi_real vector containing the binnend values of the
   *  real-space correlation function 
   *
   *  @param xi_ vector containing the binnend values of \f$
   *  \overline{\xi}(r) \f$
   *
   *  @param xi__ vector containing the binnend values of \f$
   *  \overline{\overline{\xi}}(r) \f$
   *
   *  @param var 1/[H(z)a(z)]
   *
   *  @param FV 0 &rArr; exponential; &rArr; 1 gaussian 
   *
   *  @param index index for internal use
   *
   *  @param bias_nl 0 &rArr; linear bias; &rArr; 1 non-linear bias 
   *
   *  @param bA the parameter b<SUB>A</SUB> used to model the bias
   *
   *  @param v_min the minimum value of the velocity used in the
   *  convolution
   *
   *  @param v_max the maximum value of the velocity used in the
   *  convolution
   *
   *  @param step_v the step of the convolution integral
   *
   *  @return &xi;(r<SUB>p</SUB>,&pi;)
   */
  double xi2D_model (const double, const double, const double, const double, const double, const vector<double>, const vector<double>, const vector<double>, const vector<double>, const double, const int, int index=-1, const bool bias_nl=0, const double bA=0., const double v_min=-3000., const double v_max=3000., const int step_v=500);


  /**
   *  @brief pairwise velocity distribution
   *  @param vel comoving velocity
   *  @param sigma12 &sigma;<SUB>12</SUB>
   *  @param FV 0 &rArr; exponential; &rArr; 1 gaussian 
   *  @return f(v)
   */
  double f_v (const double, const double, const int);

  /**
   *  @brief pairwise velocity distribution
   *
   *  (see Chuang&Wang 2012 (1209.0210), eq. 25)
   *
   *  @param vel comoving velocity
   *
   *  @param rp r<SUB>p</SUB>: comoving separation perpendicular to
   *  the line-of-sight
   *
   *  @param pi &pi;: comoving separation parallel to the
   *  line-of-sight
   *
   *  @param var 1/[H(z)a(z)]
   *  @param sigmav0 &sigma;<SUB>v,0</SUB>
   *  @param cmu C<SUB>&mu;</SUB>
   *  @param cs1 C<SUB>&sigma;1</SUB>
   *  @param cs2 C<SUB>&sigma;2</SUB>
   *
   *  @return f(v)
   */
  double f_v (const double, const double, const double, const double, const double, const double, const double, const double);

  /**
   *  @brief velocity distribution used to model BAO
   *  
   *  (see Chuang&Wang 2012, 1209.0210)
   *
   *  @param xx s
   *  @param f_g f<SUB>g</SUB>
   *  @param k_star k<SUB>*</SUB>
   *
   *  @return f<SUB>*</SUB>
   */
  double f_star (const double, const double, const double);

  /**
   *  @brief a possible parameterization of the non-linear bias
   *  
   *  the non-linear bias (see Chuang&Wang 2012 (1209.0210),
   *  eqs. 20-21): \f$ b(r) = r^{\frac{b_A}{1+(\frac{r}{b_B})^{b_C}}}
   *  \f$
   *
   *  @param rr the comoving scale
   *  @param bA b<SUB>A</SUB>
   *  @param bB b<SUB>B</SUB>
   *  @param bC b<SUB>C</SUB>
   *
   *  @return b(r)
   */
  double b_nl (const double, const double, const double bB=10., const double bC=4.);

  /**
   *  @brief estimated relative error on \f$\beta=f/b\f$
   *
   *  \f$
   *  \delta\beta/\beta\sim Cb^{0.7}V^{-0.5}
   *  \exp\left[n_0/(b^2n)\right]
   *  \f$
   *
   *  where \f$n_0=1.7\times10^{-4} h^3 Mpc^{-3}\f$ and
   *  \f$C=4.9\times10^2 h^{-1.5} Mpc^{1.5}\f$ 
   *
   *  see Eq. 20 of Bianchi et al. 2012, http://arxiv.org/abs/1203.1545
   *
   *  @param bias the linear galaxy bias, b
   *  @param Volume the survey volume, V
   *  @param density the galaxy density, n
   *  @return \f$\delta\beta/\beta\f$
   */
  double relative_error_beta (const double, const double, const double); 

  ///@}


  // ============================================================================================


  /**
   *  @brief The namespace of the global classes, structures and
   *  functions of the CosmoBolognaLib of internal auxiliary use
   *  
   *  The \e glob namespace contains all the global classes,
   *  structures and functions structures of the CosmoBolognaLib of
   *  internal auxiliary use
   */
  namespace glob {

    struct STR_grid
    {
      vector<double> _xx, _yy;
    };

    struct STR_grid_2D
    {
      vector<double> _xx1, _xx2;
      vector<vector<double> > _yy;
    };

    struct STR_xi0_model
    {
      double bias_sigma8;
      double sigma8z;
      vector<double> xi_DM;
    };

    struct STR_xi2D_model
    {
      vector<double> rp, pi, xi_real, xi_, xi__, P2, P4, vel;
      vector<int> lim_index_fit, type;
      int step_v, FV, dim;
      double delta_v;
      bool bias_nl;
    };

    struct STR_xi 
    {
      double rr, aa;
      vector<double> lgkk, lgPk;
    };

    struct STR_SSM 
    {
      int unit;
      double hh, mass, rho, n_spec;
      vector<double> lgkk, lgPk;
    };
  }

}

#include "FuncClassFunc.h"
#include "Chi2ClassFunc.h"

#endif
