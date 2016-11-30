// ===================================================================================================================
// Example code: how to model the Cartesian 2D two-point correlation function to constrain the linear bias and sigma12
// ===================================================================================================================

#include "Modelling_TwoPointCorrelation_cartesian.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // --------------------------------------------------------------------------------
    // ---------------- use default cosmological parameters and set sigma8 ------------
    // --------------------------------------------------------------------------------
  
    cosmobl::cosmology::Cosmology cosmology;
    cosmology.set_sigma8(0.8);

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------

    const string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";
  
    const cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

  
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------
    // ---------------- measure the projected two-point correlation function ----------------
    // --------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rMin = 0.1;  // minimum separation 
    const double rMax = 50.;  // maximum separation 
    const int nbins = 25;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre 

    const string dir = cosmobl::par::DirLoc+"../output/";
    const string file = "xi2D.dat";

  
    // measure the 2D Cartesian two-point correlation function and estimate Poissonian errors

    const auto TwoP = cosmobl::twopt::TwoPointCorrelation::Create(cosmobl::twopt::TwoPType::_2D_Cartesian_, catalogue, random_catalogue, cosmobl::_linear_, rMin, rMax, nbins, shift, cosmobl::_linear_, rMin, rMax, nbins, shift);

    TwoP->measure(cosmobl::twopt::_Poisson_, dir, {dir});
    TwoP->write(dir, file);
  

    // ----------------------------------------------------------------------------------------------------------------------------------
    // ----------------- model the Cartesian 2D two-point correlation function and estimate the linear bias and sigma12 ----------------- 
    // ----------------------------------------------------------------------------------------------------------------------------------

    cosmobl::modelling::Modelling_TwoPointCorrelation_cartesian model_twop(TwoP); // object used for modelling

  
    // flat prior for f*sigma8
    const vector<double> fsigma8_limits = {0., 1.}; 
    const cosmobl::statistics::Prior fsigma8_prior {cosmobl::statistics::PriorType::_UniformPrior_, fsigma8_limits[0], fsigma8_limits[1]}; 
  
    // flat prior for b*sigma8
    const vector<double> bsigma8_limits = {0.8*cosmology.sigma8(), 3.*cosmology.sigma8()}; 
    const cosmobl::statistics::Prior bsigma8_prior {cosmobl::statistics::PriorType::_UniformPrior_, bsigma8_limits[0], bsigma8_limits[1]}; 
  
    // flat prior for sigma12
    const vector<double> sigma12_limits = {1., 1000.}; 
    const cosmobl::statistics::Prior sigma12_prior {cosmobl::statistics::PriorType::_UniformPrior_, sigma12_limits[0], sigma12_limits[1]}; 

    const vector<double> model_scales = cosmobl::linear_bin_vector(100, 0.1, 50.); // scales used in the fit
    const double redshift = 1.; // mean redshift of the sample
  
    model_twop.set_parameters_xi2D_DM(model_scales, cosmology, redshift); // set the model parameters

    model_twop.set_model_dispersionModel(fsigma8_prior, bsigma8_prior, sigma12_prior); // set the model for the redshift-space 2D two-point correlation 

  
    // sample the likelihood and store the chains
  
    const int nChains = 10; // number of chains
    const int chain_size = 100; // size of the chains
    const string chain_file = "chain_cartesian_bias_sigma12.dat";

    model_twop.sample_likelihood(1., 40., 1., 40., cosmobl::statistics::_GaussianLikelihood_Error_, nChains, chain_size, 32113, dir, chain_file); 

  }
  
  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

