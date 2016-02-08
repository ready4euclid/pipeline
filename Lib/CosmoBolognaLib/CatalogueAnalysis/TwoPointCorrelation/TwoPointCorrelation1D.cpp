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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D.h"
#include "Data1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_pairs (const shared_ptr<pairs::Pair> PP, const string dir, const string file) const 
{  
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  
  fout.setf(ios::fixed);

  for (int i=0; i<PP->nbins(); i++) 
    fout << setprecision(0) << PP->PP1D(i) << endl;
  
  fout.clear(); fout.close(); cout << "I wrote the file " << file_out << endl;
  
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_pairs (shared_ptr<pairs::Pair> PP, const vector<string> dir, const string file) const
{
  if (dir.size()==0)
    ErrorMsg ("Error in cosmobl::twopt::TwoPointCorrelation1D::read_pairs of TwoPointCorrelation1D.cpp! dir.size()=0!");
      
  for (size_t dd=0; dd<dir.size(); dd++) {
        
    string file_in = dir[dd]+file; 
    cout << "I'm reading the pair file: " << file_in << endl;
    
    ifstream fin(file_in.c_str()); checkIO(file_in, 1);
   
    double pp;
    for (int i=0; i<PP->nbins(); i++) {
      fin >>pp;
      PP->add_PP1D(i, pp);
    }
    
    fin.clear(); fin.close(); cout << "I read the file " << file_in << endl;
  }
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::write_pairs (const vector<shared_ptr<pairs::Pair> > PP, const string dir, const string file) const 
{  

  int nRegions = m_data->get_region_list().size();
  string MK = "mkdir -p "+dir; if (system (MK.c_str())) {};
  
  string file_out = dir+file;  
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);
  fout.setf(ios::fixed);

  for (int i=0; i<nRegions; i++) {
   
    for (int j=i; j<nRegions; j++) {

      int index = i*nRegions-(i-1)*i/2+j-i;

      for (int r1=0; r1<PP[index]->nbins(); r1++)
	if(PP[index]->PP1D(r1)>0)
	  fout << i << " " << j << " " << r1 << " " << setprecision(0) << PP[index]->PP1D(r1) << endl;

    }
  }

  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation1D::read_pairs (vector<shared_ptr<pairs::Pair> > PP, const vector<string> dir, const string file) const
{

  //int nRegions = 0.5*(-1+sqrt(1+PP.size()*8));
  int nRegions = m_data->get_region_list().size();
  for (unsigned int dd=0; dd<dir.size(); dd++) {

    string ff = dir[dd]+file; 
    cout << "I'm reading the pair file: " << ff << endl;
    ifstream fin(ff.c_str()); checkIO(ff, 1);

    int I, J, bin1, index;
    double pairs;

    while (fin >> I >> J >> bin1 >> pairs) {
      index = I*nRegions-(I-1)*I/2+J-I;
    //  cout << I << " " << J << " " << index << " " << bin1 << " " << pairs << endl;
      PP[index]->add_PP1D(bin1,pairs);
    }
    fin.clear(); fin.close(); cout << "I read the file " << ff << endl;
  }

}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation1D::NaturalEstimatorTwoP (shared_ptr<pairs::Pair> dd, shared_ptr<pairs::Pair> rr, int nData, int nRandom){

  vector<double> rad,xi,error;
  rad.resize(m_dd->nbins()); xi.resize(m_dd->nbins()); error.resize(m_dd->nbins());
  
  double norm = double(nRandom)*double(nRandom-1)/(nData*double(nData-1));

  for (int i=0; i<dd->nbins(); i++) {
  
    rad[i] = dd->scale(i);
    xi[i] = -1.;
    error[i] = 1.e30;
    
    if (dd->PP1D(i)>0 && rr->PP1D(i)>0) {
      
      xi[i] = max(-1., norm*dd->PP1D(i)/rr->PP1D(i)-1.);

      error[i] = PoissonError(dd->PP1D(i), rr->PP1D(i), 0); 
    }
  }

  return move(unique_ptr<Data1D>(new Data1D(rad,xi,error)));
}


// ============================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation1D::LandySzalayEstimatorTwoP (shared_ptr<pairs::Pair> dd, shared_ptr<pairs::Pair> rr, shared_ptr<pairs::Pair> dr, int nData, int nRandom){

  vector<double> rad,xi,error;
  rad.resize(m_dd->nbins()); xi.resize(m_dd->nbins()); error.resize(m_dd->nbins());
  
  double norm = double(nRandom)*double(nRandom-1)/(nData*double(nData-1));
  double norm1 = double(nRandom-1)/nData;

  for (int i=0; i<dd->nbins(); i++) {
  
    rad[i] = dd->scale(i);
    xi[i] = -1.;
    error[i] = 1.e30;
    
    if (dd->PP1D(i)>0 && rr->PP1D(i)>0) {
      
      xi[i] = max(-1., norm*dd->PP1D(i)/rr->PP1D(i)-norm1*dr->PP1D(i)/rr->PP1D(i)+1.);

      error[i] = PoissonError(dd->PP1D(i), rr->PP1D(i), dr->PP1D(i)); 
    }
  }

  return move(unique_ptr<Data1D>(new Data1D(rad,xi,error)));

}
