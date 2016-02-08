// =================================================================================================
// Example code: how to measure the angle-averaged two-point correlation function, i.e. the monopole 
// =================================================================================================

#include "RandomCatalogue.h"
#include "TwoPointCorrelation1D_monopole.h"

using namespace cosmobl;
using namespace catalogue;
using namespace twopt;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // -----------------------------------------------------------------
  // ---------------- use default cosmological parameters ------------
  // -----------------------------------------------------------------

  Cosmology cosmology;

  
  // ----------------------------------------------------------
  // ---------------- read the input catalogue ----------------
  // ----------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";

  Catalogue catalogue {file_catalogue, cosmology, _Galaxy_};

  
  // ----------------------------------------------------------------
  // ---------------- construct the random catalogue ----------------
  // ----------------------------------------------------------------

  double N_R = 1.; // random/object ratio
  Catalogue random_catalogue {catalogue, N_R};

  
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
  TwoP.measure(par::DirLoc+"../output/");

  
  // store the output data
  
  string dir_output = par::DirLoc+"../output/";
  string file_xi = "xi.dat";
  TwoP.write(dir_output, file_xi);
  
  
  return 0;
}

