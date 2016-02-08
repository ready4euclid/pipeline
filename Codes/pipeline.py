#==========================================#
# The pipeline of the ready4euclid project #
#==========================================#


###############################################################################

#------------------------------------#
# import the required Python modules #
#------------------------------------#

import os
import numpy as np
import CosmoBolognaLib as cbl 


###############################################################################

#-------------------#
# input/output data #
#-------------------#

HOME = os.getenv("HOME")
file_data = HOME+"/Input/data.dat"
file_catalogue = HOME+"/Input/clusterCatalog.dat"
dir_output = HOME+"/Output/"


###############################################################################

#-------------------------------#
# define the cosmological model #
#-------------------------------#

cosmology = cbl.Cosmology()


###############################################################################

#-----------------#
# detect clusters #
#-----------------#

#cbl.detectClusters(file_data, file_catalogue) [to be implemented/included]


###############################################################################

#----------------------------#
# read the cluster catalogue # 
#----------------------------#

catalogue = cbl.Catalogue(file_catalogue, cosmology, cbl._Galaxy_)


###############################################################################

#-----------------------------------#
# measure the cluster mass function #
#-----------------------------------#

#MF  = cbl.measure_MF(catalogue) [to be included]


###############################################################################

#--------------------------------------------#
# measure the two-point correlation function #
#--------------------------------------------#

# binnig parameters #
rMin = 1.   # minimum separation 
rMax = 50.  # maximum separation 
nbins = 20  # number of bins
shift = 0.5 # spatial shift used to set the bin centre 

# construct the random catalogue #
random_catalogue = cbl.Catalogue(catalogue, 100)

# create the object used to measure the two-point correlation function #
TwoP = cbl.TwoPointCorrelation1D_monopole(catalogue, random_catalogue, cbl._logarithmic_, rMin, rMax, nbins, shift)

# measure the two-point correlation function #
TwoP.measure(dir_pairs)
  
# store the output data #
TwoP.write(dir_output, "xi.dat")


###############################################################################

#---------------------------------#
# extract cosmological parameters #
#---------------------------------#

#[to be implemented/included]
