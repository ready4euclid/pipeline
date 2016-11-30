// =================================================================================
// Example code: how to use the NormalRandomNumbers class to generate random numbers
// =================================================================================

#include "Func.h" 

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {
  
    // choose a seed
    int seed = 1231;

    // choose a range 
    double xmin = 0.;
    double xmax = 10.;

    // choose the parameters for the Normal distribution
    double mean = 5.;
    double sigma = 2.;

    // construct the object used to extract random numbers from a Normal distribution
    cosmobl::random::NormalRandomNumbers ran(mean, sigma, seed, xmin, xmax);

    // extract the random numbers
    vector<double> numbers(1000000);
    for (auto && num : numbers) num = ran();

    // compute the distribution of the extracted random numbers
    vector<double> xx, fx, error;
    cosmobl::distribution(xx, fx, error, numbers, {}, 10);
  
    // show the results
    for (size_t i=0; i<xx.size(); ++i)
      cout << setprecision(1) << fixed << xx[i] << "  " << string(fx[i]*100./numbers.size(), '*') << endl;

  }
  
  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}
