# pipeline for cosmological analyses of galaxy cluster catalogs

The implemented C++/Python software (will) perform the following analyses :
- construction of the cluster catalog(s) 
- calibration of scaling relations
- measurement of cluster number counts and associated errors
- measurement of 2pt statistics and associated errors
- full likelihood analysis and extraction of cosmological parameters

Github folders: 
Lib/  :  the internal libraries
Codes/  :  the pipeline, that is a single python code: pipeline.py
Input/  :  the input data and cluster catalog(s)
Output/  :  all the output files
 
At present, the only internal libraries are the CosmoBolognaLib, and a slightly modified version of the Numerical recipes in C++.
The CosmoBolognaLib are implemented in C++, and then converted to python using SWIG.


- the required external libraries are the following:
GSL (GNU Scientific Library)
FFTW (Fastest Fourier Transform in the West)
OpenMP (Open Multiprocessing)

