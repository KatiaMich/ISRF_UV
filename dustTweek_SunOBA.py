import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import FluxCalculator_OBA_dustTweek as fc

x0, y0, z0 = 0,0,0 
co = np.array([x0,y0,z0])

Catalog_cd = 'all_Catalogs/'
cwd = '/home/katia/Computation/'
wavelength = 1565
dustTweeks = np.linspace(0,2,11)

flux_lambda = []

def DefineMatrix():
    X = np.arange(-3000,3001,5) #dimension of the dust map, one value every 5 pc
    Y = np.arange(-3000,3001,5)
    Z = np.arange(-400,401,5)
    matrixDust = np.load(cwd+'inputs/matrixDust.npy')
    interpolated_dustMap = RegularGridInterpolator((X,Y,Z),matrixDust,bounds_error=False,fill_value=0) 
    return interpolated_dustMap

interpolated_dustMap = DefineMatrix()

def Compute_summed_flux(dustTweek):
    #Stars in Hipparcos, not in StarHorse 
    HIPvsSH = 'HIP'

    filename = cwd+Catalog_cd+'HIPL_2023_cleanedBinaries'
    pdtable = pd.read_pickle(filename)
    fluxesO, fluxesB, fluxesA = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd,dustTweek) #this is a long array
    flux0O = np.sum(fluxesO)
    flux0B = np.sum(fluxesB)
    flux0A = np.sum(fluxesA)
    
    #Stars is StarHorse, first the ones also in Hipparcos, where we use the Hipparcos temperature, then the ones with StarHorse temperature.
    HIPvsSH = 'SH'

    filename = cwd+Catalog_cd+'SH_2023_external_SpType'
    pdtable = pd.read_pickle(filename)
    fluxesO,fluxesB, fluxesA = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd,dustTweek) #this is a long array
    flux1O = np.sum(fluxesO)
    flux1B = np.sum(fluxesB)
    flux1A = np.sum(fluxesA)
    

    filename = cwd+Catalog_cd+'SH_2023_noSpType_above20000'
    pdtable = pd.read_pickle(filename)
    fluxesO,fluxesB, fluxesA = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd,dustTweek) #this is a long array
    flux2O = np.sum(fluxesO)
    flux2B = np.sum(fluxesB)
    flux2A = np.sum(fluxesA)
    
    fluxTotO = flux0O+flux1O+flux2O
    fluxTotB = flux0B+flux1B+flux2B
    fluxTotA = flux0A+flux1A+flux2A
    
    return np.array([fluxTotO, fluxTotB, fluxTotA])

for k in range(len(dustTweeks)):
    dustTweek = dustTweeks[k]
    flux_lambda.append(np.array(Compute_summed_flux(dustTweek)))
flux_lambda = np.transpose(flux_lambda)

df = pd.DataFrame({'dustTweek':dustTweeks,'summed_fluxO':flux_lambda[0],'summed_fluxB':flux_lambda[1],'summed_fluxA':flux_lambda[2], 'summed_flux':np.sum(flux_lambda,axis=0)})
df.to_pickle('resultdustTweek/dustTweekSun2')