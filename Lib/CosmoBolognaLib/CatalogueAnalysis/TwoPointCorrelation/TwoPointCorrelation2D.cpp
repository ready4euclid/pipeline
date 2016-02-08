/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation2D used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation2D used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "TwoPointCorrelation2D.h"
#include "Data2D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ===================================================================================================


void cosmobl::twopt::TwoPointCorrelation2D::write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  
  fout.setf(ios::fixed);

  for (int i=0; i<PP->nbins_D1(); i++)
    for (int j=0; j<PP->nbins_D2(); j++) 
      fout << setprecision(0) << PP->PP2D(i, j) << endl;
  
  fout.clear(); fout.close(); cout << "I wrote the file " << file_out << endl;
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const
{
  if (dir.size()==0)
    ErrorMsg ("Error in cosmobl::twopt::TwoPointCorrelation2D::read_pairs of TwoPointCorrelation2D.cpp! dir.size()=0!");
      
  for (size_t dd=0; dd<dir.size(); dd++) {
        
    string file_in = dir[dd]+file; 
    cout << "I'm reading the pair file: " << file_in << endl;
    
    ifstream fin(file_in.c_str()); checkIO(file_in, 1);
   
    double pp;
    for (int i=0; i<PP->nbins_D1(); i++) 
      for (int j=0; j<PP->nbins_D2(); j++) {
	fin >>pp;
	PP->add_PP2D(i, j, pp);
      }
    
    fin.clear(); fin.close(); cout << "I read the file " << file_in << endl;
  }
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::write_pairs (const vector<shared_ptr<pairs::Pair> > PP, const string dir, const string file) const 
{  

  int nRegions = 0.5*(-1+sqrt(1+PP.size()*8));

  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;  
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  fout.setf(ios::fixed);

  for (int i=0; i<nRegions; i++) {
   
    for (int j=i; j<nRegions; j++) {

      int index = i*nRegions-(i-1)*i/2+j-i;
      
      for (int r1=0; r1<PP[index]->nbins_D1(); r1++)
	for (int r2=0; r2<PP[index]->nbins_D2(); r2++)
	  if(PP[index]->PP2D(r1,r2)>0)
	    fout << i << " " << j << " " << r1 << " " << r2 << " " << setprecision(0) << PP[index]->PP2D(r1,r2) << endl;

    }
  }

  fout.clear(); fout.close();
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation2D::read_pairs (vector<shared_ptr<pairs::Pair> > PP, const vector<string> dir, const string file) const
{

  int nRegions = 0.5*(-1+sqrt(1+PP.size()*8));
  for (unsigned int dd=0; dd<dir.size(); dd++) {

    string ff = dir[dd]+file; 
    cout << "I'm reading the pair file: " << ff << endl;
    ifstream fin(ff.c_str()); checkIO(ff, 1);

    int I, J, bin1, bin2, index;
    double pairs;

    while (fin >> I >> J >> bin1 >> bin2 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
      PP[index]->add_PP2D(bin1,bin2,pairs);
    }
    fin.clear(); fin.close();

    fin.clear(); fin.close(); cout << "I read the file " << ff << endl;

  }
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation2D::NaturalEstimatorTwoP (shared_ptr<pairs::Pair> dd, shared_ptr<pairs::Pair> rr, const int nData, const int nRandom) {

  vector<double> scale_D1,scale_D2;
  vector<vector<double> > xi,error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));
  error.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));

  double norm = double(nRandom)*double(nRandom-1)/(nData*double(nData-1));

  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);

      xi[i][j] = -1.;
      error[i][j] = 1.e30;

      if (dd->PP2D(i,j)>0 && rr->PP2D(i,j)>0) {

	xi[i][j] = max(-1., norm*dd->PP2D(i,j)/rr->PP2D(i,j)-1.);

	error[i][j]= PoissonError(dd->PP2D(i,j), rr->PP2D(i,j), 0); 
      }
    }
  }

  return move(unique_ptr<Data2D>(new Data2D(scale_D1,scale_D2,xi,error)));
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation2D::LandySzalayEstimatorTwoP (shared_ptr<pairs::Pair> dd, shared_ptr<pairs::Pair> rr, shared_ptr<pairs::Pair> dr, int nData, int nRandom) {

  vector<double> scale_D1,scale_D2;
  vector<vector<double> > xi,error;
  scale_D1.resize(m_dd->nbins_D1()); 
  scale_D2.resize(m_dd->nbins_D2()); 

  xi.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));
  error.resize(m_dd->nbins_D1(),vector<double>(m_dd->nbins_D2(),0));
  
  double norm = double(nRandom)*double(nRandom-1)/(nData*double(nData-1));
  double norm1 = double(nRandom-1)/nData;

  for (int i=0; i<dd->nbins_D1(); i++) {
    scale_D1[i] = dd->scale_D1(i);
    for (int j=0; j<dd->nbins_D2(); j++) {
      scale_D2[j] = dd->scale_D2(j);

      xi[i][j] = -1.;
      error[i][j] = 1.e30;

      if (dd->PP2D(i,j)>0 && rr->PP2D(i,j)>0) {

	xi[i][j] = max(-1., norm*dd->PP2D(i,j)/rr->PP2D(i,j)-norm1*dr->PP2D(i,j)/rr->PP2D(i,j)+1.);

	error[i][j]= PoissonError(dd->PP2D(i,j), rr->PP2D(i,j), dr->PP2D(i,j)); 
      }
    }
  }

  return move(unique_ptr<Data2D>(new Data2D(scale_D1,scale_D2,xi,error)));
}
