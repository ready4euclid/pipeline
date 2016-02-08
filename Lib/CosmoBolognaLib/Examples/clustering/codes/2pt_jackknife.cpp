// ===================================================================================
// How to measure the two-point correlation function and estimate the jackknife errors
// ===================================================================================

#include "RandomCatalogue.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "GlobalFunc.h"

using namespace cosmobl;
using namespace catalogue;
using namespace twopt;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;

int main () {
  
  // -----------------------------------------------------------------------------------------------
  // ---------------- set the cosmological parameters and create the object 'cosmology' ------------
  // -----------------------------------------------------------------------------------------------


  Cosmology cosmology;

  
  // --------------------------------------------------------------------
  // ---------------- Input/Output files and directories ----------------
  // --------------------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";

  string dir_output = par::DirLoc+"../output3/";
  string dir_pairs = dir_output+"pairs/";
  string dir_random_cat = dir_output;
  string dir_covariance = dir_output+"covariance/";
  
  string MK = "mkdir -p "+dir_output+" "+dir_pairs+" "+dir_covariance; if (system(MK.c_str())) {};

  
  // -------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue and create the object 'catalogue'----------------
  // -------------------------------------------------------------------------------------------

  cout << "I'm reading the input catalogue..." << endl;

  Catalogue catalogue {file_catalogue, cosmology, _Galaxy_};

  
  // ----------------------------------------------------------------
  // ---------------- construct the random catalogue ----------------
  // ----------------------------------------------------------------

  double N_R = 1.; // random/object ratio
  Catalogue random_catalogue {catalogue, N_R};

  int nx=5,ny=5,nz=5;
  set_ObjectRegion_SubBoxes(catalogue,random_catalogue,nx,ny,nz);

  // --------------------------------------------------------------------------------------------
  // ---------------- measure the monopole of the two-point correlation function ----------------
  // --------------------------------------------------------------------------------------------

  // binning parameters

  double rMin = 1.;   // minimum separation 
  double rMax = 50.;  // maximum separation 
  int nbins = 20;     // number of bins
  double shift = 0.5; // spatial shift used to set the bin centre 
      
  // measure the monopole

  TwoPointCorrelation1D_monopole TwoP {catalogue, random_catalogue, _logarithmic_, rMin, rMax, nbins, shift};

  TwoP.measure(dir_output);

  TwoP.write(dir_output,"xi_Poisson.dat");

  TwoP.measure(dir_output,{},ErrorType::_Jackknife_);

  TwoP.write(dir_output,"xi_Jackknife.dat");

  TwoP.measure(dir_output,{},ErrorType::_Bootstrap_,"NULL",100);

  TwoP.write(dir_output,"xi_Bootstrap.dat");

  return 0;
}

