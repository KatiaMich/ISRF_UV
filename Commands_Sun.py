import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import FluxCalculator as fc #python file FluxCalculator.py should be in the same directory

##############these are the only lines you have to change##########################
wavelength = 1565 #Angstrom 
Catalog_cd = 'all_Catalogs/' #path to the Catalog directory
cwd = '/home/katia/Computation/'  #current directory, where also the input directory is located, and where you should have created an empty directory for the results: 'result'
###################################################################################


x0, y0, z0 = 0,0,0 
co = np.array([x0,y0,z0]) #galactic heliocentric coordinates of the Sun, where we will compute the ISRF

def DefineMatrix_dustMap():
    X = np.arange(-3000,3001,5) #dimension of the dust map, one value every 5 pc
    Y = np.arange(-3000,3001,5)
    Z = np.arange(-400,401,5)
    matrixDust = np.load(cwd+'inputs/matrixDust.npy')
    interpolated_dustMap = RegularGridInterpolator((X,Y,Z),matrixDust,bounds_error=False,fill_value=0) 
    return interpolated_dustMap    
interpolated_dustMap = DefineMatrix_dustMap()


###############Stars in Hipparcos, not in StarHorse 
HIPvsSH = 'HIP'
filename = cwd+Catalog_cd+'HIPL_2023_cleanedBinaries'
pdtable = pd.read_pickle(filename)
fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
pdFluxes = pd.DataFrame({'flux':fluxes})
df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
df.to_pickle('result/HIPL_2023_Sun')


###############Stars is StarHorse, first the ones with a known spectral type from an external catalog
HIPvsSH = 'SH'
filename = cwd+Catalog_cd+'SH_2023_external_SpType'
pdtable = pd.read_pickle(filename)
fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
pdFluxes = pd.DataFrame({'flux':fluxes})
df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
df.to_pickle('result/SH_extSpType_2023_Sun')

###############Stars is StarHorse, the ones without a spectral type, where we use StarHorse temperature.
HIPvsSH = 'SH'
filename = cwd+Catalog_cd+'SH_2023_noSpType_above20000'
pdtable = pd.read_pickle(filename)
fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
pdFluxes = pd.DataFrame({'flux':fluxes})
df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
df.to_pickle('result/SH_noSpType_2023_Sun')
