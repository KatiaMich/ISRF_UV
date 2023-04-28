import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import FluxCalculator as fc

x0, y0, z0 = 0,0,0 
co = np.array([x0,y0,z0])

Catalog_cd = 'all_Catalogs/'
cwd = '/home/katia/Computation/'
wavelength = 1565

def DefineMatrix():
    X = np.arange(-3000,3001,5) #dimension of the dust map, one value every 5 pc
    Y = np.arange(-3000,3001,5)
    Z = np.arange(-400,401,5)
    matrixDust = np.load(cwd+'inputs/matrixDust.npy')
    interpolated_dustMap = RegularGridInterpolator((X,Y,Z),matrixDust,bounds_error=False,fill_value=0) 
    return interpolated_dustMap

    
interpolated_dustMap = DefineMatrix()


#Stars in Hipparcos, not in StarHorse 
HIPvsSH = 'HIP'

filename = cwd+Catalog_cd+'HIPL_2023_cleanedBinaries'
pdtable = pd.read_pickle(filename)
fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
pdFluxes = pd.DataFrame({'flux':fluxes})
df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
df.to_pickle('result_halfext/HIPL_2023_Sun')


#Stars is StarHorse, first the ones also in Hipparcos, where we use the Hipparcos temperature, then the ones with StarHorse temperature.
HIPvsSH = 'SH'

filename = cwd+Catalog_cd+'SH_2023_external_SpType'
pdtable = pd.read_pickle(filename)
fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
pdFluxes = pd.DataFrame({'flux':fluxes})
df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
df.to_pickle('result_halfext/SH_extSpType_2023_Sun')

filename = cwd+Catalog_cd+'SH_2023_noSpType_above20000'
pdtable = pd.read_pickle(filename)
fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
pdFluxes = pd.DataFrame({'flux':fluxes})
df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
df.to_pickle('result_halfext/SH_noSpType_2023_Sun')
