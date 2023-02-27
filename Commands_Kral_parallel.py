import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from multiprocessing import Pool

import FluxCalculator as fc

Catalog_cd = 'all_Catalogs/'
cwd = '/home/katia/Computation/'
wavelength = 1565
KralCoord = np.load(cwd+'inputs/KralCoord.npy')

def DefineMatrix():
    X = np.arange(-3000,3001,5) #dimension of the dust map, one value every 5 pc
    Y = np.arange(-3000,3001,5)
    Z = np.arange(-400,401,5)
    matrixDust = np.load(cwd+'inputs/matrixDust.npy')
    interpolated_dustMap = RegularGridInterpolator((X,Y,Z),matrixDust,bounds_error=False,fill_value=0) 
    return interpolated_dustMap
   
interpolated_dustMap = DefineMatrix()


def Kral_Comput(k):
    co = np.array(KralCoord[k])
    
    #Stars in Hipparcos, not in StarHorse 
    HIPvsSH = 'HIP'

    filename = cwd+Catalog_cd+'HIPL_2023_OBA'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    df.to_pickle('result/HIPL_2023_Kral'+str(k))


    #Stars is StarHorse, first the ones also in Hipparcos, where we use the Hipparcos temperature, then the ones with StarHorse temperature.
    HIPvsSH = 'SH'

    filename = cwd+Catalog_cd+'SH_HIPtemp_OBA'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    df.to_pickle('result/SH_HIPtemp_2023_Kral'+str(k))

    filename = cwd+Catalog_cd+'SH_2023_noSpType'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    df.to_pickle('result/SH_noSpType_2023_Kral'+str(k))
    return None           
    
if __name__ == "__main__":
    Pool().map(Kral_Comput, range(0,30))
    