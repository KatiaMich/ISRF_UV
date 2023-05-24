import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

import FluxCalculator as fc

x0, y0, z0 = 0,0,0 
co = np.array([x0,y0,z0])

Catalog_cd = 'all_Catalogs/'
cwd = '/home/katia/Computation/'
wavelengthsTD1 = [912,1565, 1965, 2365, 2740]
wavelength00 = np.arange(920, 1500, 10)
wavelengths1  = np.arange(1500,2000,10)
#wavelengths2  = np.arange(2000,2900,100)
wavelengths0 = np.concatenate((wavelengthsTD1,wavelength00, wavelengths1))#,wavelengths2))
wavelengths = np.sort(wavelengths0)



flux_lambda = np.zeros(len(wavelengths))

def DefineMatrix():
    X = np.arange(-3000,3001,5) #dimension of the dust map, one value every 5 pc
    Y = np.arange(-3000,3001,5)
    Z = np.arange(-400,401,5)
    matrixDust = np.load(cwd+'inputs/matrixDust.npy')
    interpolated_dustMap = RegularGridInterpolator((X,Y,Z),matrixDust,bounds_error=False,fill_value=0) 
    return interpolated_dustMap

interpolated_dustMap = DefineMatrix()

def Compute_summed_flux(wavelength):
    #Stars in Hipparcos, not in StarHorse 
    HIPvsSH = 'HIP'

    filename = cwd+Catalog_cd+'HIPL_2023_OBA'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    flux0 = np.sum(fluxes)

    #Stars is StarHorse, first the ones also in Hipparcos, where we use the Hipparcos temperature, then the ones with StarHorse temperature.
    HIPvsSH = 'SH'

    filename = cwd+Catalog_cd+'SH_HIPtemp_OBA'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    flux1 = np.sum(fluxes)

    filename = cwd+Catalog_cd+'SH_2023_noSpType'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    flux2 = np.sum(fluxes)
    fluxTot = flux0+flux1+flux2
    return fluxTot

for k in range(len(wavelengths)):
    wavelength = wavelengths[k]
    flux_lambda[k] = Compute_summed_flux(wavelength)
    
df = pd.DataFrame({'wavelength':wavelengths,'summed_flux':flux_lambda})
df.to_pickle('spectrum/spectrum_Sun')