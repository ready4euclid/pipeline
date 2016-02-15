#!/usr/bin/env python

from setuptools import setup, Extension
import os

from distutils.sysconfig import get_config_vars

(opt,) = get_config_vars('OPT')
os.environ['OPT'] = " ".join(flag for flag in opt.split() if flag != '-Wstrict-prototypes')

def readme():
    with open('README.rst') as f:
        return f.read()

HOME = os.getenv("HOME")

dirLib       = "CosmoBolognaLib/" 
dirNumerical = HOME+"/Numerical/"
dirH         = dirLib+"Headers/Lib/"
dirO         = dirLib+"Headers/Objects/"
dirEH        = dirLib+"Cosmology/EH/"

FLAGS = ["-std=c++11", "-fopenmp", "-w"]
FLAGSL = ["-Wl,-rpath,CosmoBolognaLib/", "-LCosmoBolognaLib/"]

include_dirs = [dirLib, dirH, dirO, dirEH, dirNumerical]

libraries = ["gomp", "gsl", "gslcblas", "m", "fftw3", "fftw3_omp"]

sources = ["CBL_wrap.cxx",
           dirLib+"Func/Func.cpp",
           dirLib+"Func/FuncXi.cpp",
           dirLib+"Func/FuncMultipoles.cpp",
           dirLib+"Func/GSLfunction.cpp",
           dirLib+"Statistics/Chi2.cpp",
           dirLib+"Statistics/Chain.cpp",
           dirLib+"Statistics/Data1D_collection.cpp",
           dirLib+"Statistics/Data1D.cpp",
           dirLib+"Statistics/Data2D.cpp",
           dirLib+"Statistics/Data.cpp",
           dirLib+"Statistics/Likelihood.cpp",
           dirLib+"Statistics/Model.cpp",
           dirLib+"Statistics/Parameter.cpp",
           dirLib+"Statistics/Prior.cpp",
           dirLib+"Cosmology/EH/power_whu.cpp",
           dirLib+"Cosmology/Lib/Cosmology.cpp",
           dirLib+"Cosmology/Lib/Sigma.cpp",
           dirLib+"Cosmology/Lib/PkXi.cpp",
           dirLib+"Cosmology/Lib/PkXizSpace.cpp",
           dirLib+"Cosmology/Lib/Bias.cpp",
           dirLib+"Cosmology/Lib/RSD.cpp",
           dirLib+"Cosmology/Lib/Velocities.cpp",
           dirLib+"Cosmology/Lib/MassGrowth.cpp",
           dirLib+"Cosmology/Lib/NG.cpp",
           dirLib+"Cosmology/Lib/BAO.cpp",
           dirLib+"Cosmology/Lib/MassFunction.cpp",
           dirLib+"Cosmology/Lib/SizeFunction.cpp",
           dirLib+"ChainMesh/ChainMesh.cpp",
           dirLib+"Catalogue/Object.cpp",
           dirLib+"Catalogue/Catalogue.cpp",
           dirLib+"Catalogue/ChainMesh_Catalogue.cpp",
           dirLib+"LogNormal/LogNormal.cpp",
           dirLib+"CatalogueAnalysis/RandomCatalogue/RandomCatalogue.cpp",     
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/Pair.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_angular.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_monopole.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_filtered.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D_cartesian.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D_polar.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_projected.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_deprojected.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_multipoles.cpp",
           dirLib+"CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_wedges.cpp",
           dirLib+"CatalogueAnalysis/ThreePointCorrelation/Triplet.cpp",
           dirLib+"CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation.cpp",
           dirLib+"CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_angular_connected.cpp",
           dirLib+"CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_angular_reduced.cpp",
           dirLib+"CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_comoving_connected.cpp",
           dirLib+"CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_comoving_reduced.cpp",
           dirLib+"Modelling/Modelling.cpp",
           dirLib+"GlobalFunc/FuncCosmology.cpp",
           dirLib+"GlobalFunc/Func.cpp",
           dirLib+"GlobalFunc/SubSample.cpp"
       ]



CosmoBolognaLib = Extension(  "_CosmoBolognaLib",
                              language             = "c++",
                              sources              = sources,
                              include_dirs         = include_dirs,
                              libraries            = libraries,
                              extra_compile_args   = FLAGS,
                              extra_link_args      = FLAGSL )

setup(  name             = "CosmoBolognaLib",
        version          = "1.14",
        description      = "C++ libraries for cosmological calculations",
        long_description = readme(),
        author           = "Federico Marulli",
        author_email     = "federico.marulli3@unibo.it",
        url              = "http://github.com/federicomarulli/CosmoBolognaLib",
        license          = "GNU General Public License",
        ext_modules      = [CosmoBolognaLib],
        packages         = ["CosmoBolognaLib"],
        zip_safe         = False )

