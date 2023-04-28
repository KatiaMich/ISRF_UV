import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from multiprocessing import Pool
import time
print(time.time())
import FluxCalculator as fc

Catalog_cd = 'all_Catalogs/'
cwd = '/home/katia/Computation/'
wavelength = 1000
#KralCoord = np.load(cwd+'inputs/KralCoord.npy')
KralCatalog = pd.read_pickle(cwd+'inputs/KralCatalog2023')
x0,y0,z0 = KralCatalog['x(pc)'],KralCatalog['y(pc)'],KralCatalog['z(pc)']
x0,y0,z0 = np.array(x0),np.array(y0),np.array(z0)
KralCoord = np.array([x0,y0,z0]).transpose()
KralHIP = KralCatalog['HIP']

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

    filename = cwd+Catalog_cd+'HIPL_2023_cleanedBinaries'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    if KralHIP[k]!='None':
        df = df.loc[df['HIP']!=int(KralHIP[k])]
        print('Kral remove' + str(k)+'\n')
    flux0 = np.sum(df['flux'])
    #df.to_pickle('resultK_Vienna/HIPL_Kral'+str(k))


    #Stars is StarHorse, first the ones also in Hipparcos or Tycho, where we use the sptypes, then the ones with StarHorse temperature.
    HIPvsSH = 'SH'

    filename = cwd+Catalog_cd+'SH_2023_external_SpType'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    #df.to_pickle('resultK/SH_extSpType_Kral'+str(k))
    flux1 = np.sum(df['flux'])

    filename = cwd+Catalog_cd+'SH_2023_noSpType_above20000'
    pdtable = pd.read_pickle(filename)
    fluxes = fc.getFlux(HIPvsSH, filename ,wavelength,co,interpolated_dustMap,cwd) #this is a long array
    pdFluxes = pd.DataFrame({'flux':fluxes})
    df = pd.concat((pdtable.reset_index(drop=True),pdFluxes.reset_index(drop=True)),axis=1)
    #df.to_pickle('resultK/SH_noSpType20000_Kral'+str(k))
    flux2 = np.sum(df['flux'])
    np.save('Kral1000A/1000A_Kral'+str(k),flux0+flux1+flux2)
    return None           
    
if __name__ == "__main__":
    Pool().map(Kral_Comput, range(190))

print(time.time())    