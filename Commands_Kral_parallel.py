import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from multiprocessing import Pool

import FluxCalculator as fc

############################################### these are the only lines you should change #########################################
Catalog_cd = 'all_Catalogs/'
cwd = '/home/katia/Computation/' #you should create a 'resultKral' directory to store the results.
wavelength = 1000
####################################################################################################################################

KralCatalog = pd.read_pickle(cwd+'inputs/KralCatalog2023')
x0,y0,z0 = KralCatalog['x(pc)'],KralCatalog['y(pc)'],KralCatalog['z(pc)']
x0,y0,z0 = np.array(x0),np.array(y0),np.array(z0)
KralCoord = np.array([x0,y0,z0]).transpose()
KralHIP = KralCatalog['HIP'] #all O,B,A stars in the Kral Catalog are from Hipparcos. If KralHIP is not None, then we should remove the contribution from the host star to its ISRF.

def DefineMatrix_dustMap():
    X = np.arange(-3000,3001,5) #dimension of the dust map, one value every 5 pc
    Y = np.arange(-3000,3001,5)
    Z = np.arange(-400,401,5)
    matrixDust = np.load(cwd+'inputs/matrixDust.npy')
    interpolated_dustMap = RegularGridInterpolator((X,Y,Z),matrixDust,bounds_error=False,fill_value=0) 
    return interpolated_dustMap
interpolated_dustMap = DefineMatrix_dustMap()

def Kral_Comput(k):
    co = np.array(KralCoord[k])
    
    ###############Stars in Hipparcos, not in StarHorse 
    HIPvsSH = 'HIP'
    filename = cwd+Catalog_cd+'HIPL_2023_cleanedBinaries'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    if KralHIP[k]!='None':
        df = df.loc[df['HIP']!=int(KralHIP[k])]
        print('Kral remove' + str(k)+'\n')
    flux_HIP = np.sum(df['flux'])
    #df.to_pickle('resultKral/HIPL_Kral'+str(k)) #uncomment if you want to save the details of all contributors (500 MB for each Kral coordinate)

    ###############Stars is StarHorse, first the ones with spectral type from external surveys
    HIPvsSH = 'SH'
    filename = cwd+Catalog_cd+'SH_2023_external_SpType'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    flux_SH_extSp = np.sum(df['flux'])
    #df.to_pickle('resultKral/SH_extSpType_Kral'+str(k)) #uncomment if you want to save the details of all contributors (500 MB for each Kral coordinate)
    
    ###############Stars is StarHorse, the ones without spectral type where we use StarHorse temperatures.
    HIPvsSH = 'SH'
    filename = cwd+Catalog_cd+'SH_2023_noSpType_above20000'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    flux_SH_noSp = np.sum(df['flux'])
    #df.to_pickle('resultKral/SH_noSpType20000_Kral'+str(k)) #uncomment if you want to save the details of all contributors (500 MB for each Kral coordinate)
    
    ##############The computation is finished, let's save the summed flux (ISRF) from all stars at the coordinates of the k's Kral star.
    np.save('resultKral/summedISRF_Kral'+str(k),flux_HIP+flux_SH_extSp+flux_SH_noSp)
    return None           
    
if __name__ == "__main__":
    Pool().map(Kral_Comput, range(190))
