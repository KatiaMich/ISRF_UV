import pandas as pd
import numpy as np
import scipy as sc
from scipy import interpolate

import ReadDustMap_pd as dust

def colorFunction(g,Teff,wavelength): 
    X = np.array(g['teff50'])
    Y = np.array(g['ratio_'+str(wavelength)+'A'])
    function = sc.interpolate.interp1d(X,Y, fill_value = 'extrapolate')
    return function(Teff)

def extinctionRatio(wavelength,extinction_Function):
    ext = extinction_Function.loc[extinction_Function['wavelength']==wavelength]['Alambda_over_A_V']
    return float(ext)
    
def Vegaflux(HIPvsSH):
    if HIPvsSH == 'SH':
        #flux in the gaia band pass, svo
        return 2.50386e-9
    if HIPvsSH == 'HIP': 
        #flux in the  hipparcos band pass, svo
        return 3.76938e-9 
    
def exponent(HIPvsSH,pdtable,co,interpolated_dustMap,Alambda_overAv):
    Av , norms = dust.getAv(pdtable,co,interpolated_dustMap) #dust Extinction
    Av = np.zeros(len(Av))    
    if HIPvsSH == 'SH':
        absMag = pdtable['mg0'] #in the g band
        expo = np.array( (-0.4)*(absMag - 5 + 5*np.log10(norms) + Alambda_overAv*Av  ) ) 
        
    if HIPvsSH == 'HIP':
        absMag = pdtable['absMag_HIP'] #in the hipparcos band
        expo = np.array( (-0.4)*(absMag - 5 + 5*np.log10(norms) + Alambda_overAv*Av  ) ) 
    
    return expo #this is a numpy array

def importFiles(HIPvsSH,pdtable,cwd):
    extinction_Function = pd.read_pickle(cwd+'/inputs/extinctionFunction')
    Teff = pdtable['teff50']
    
    if HIPvsSH == 'SH':
        g = pd.read_pickle(cwd+'/inputs/ColourFunctionStarHorse.pkl')
    
    if HIPvsSH == 'HIP':
        g = pd.read_pickle(cwd+'/inputs/ColourFunctionHipparcos.pkl')
    return extinction_Function, g, Teff      

def compute(HIPvsSH,pdtable, wavelength,co,interpolated_dustMap,cwd):
    extinction_Function, g, Teff = importFiles(HIPvsSH,pdtable,cwd)   

    Alambda_overAv = extinctionRatio(wavelength,extinction_Function)
    expon = exponent(HIPvsSH, pdtable,co,interpolated_dustMap,Alambda_overAv)

    g_lambda = colorFunction(g,Teff,wavelength)
    F_Vega_detect = Vegaflux(HIPvsSH) 
    F_lambda = F_Vega_detect * g_lambda * np.exp( np.log(10) * expon )
    return F_lambda #this is a numpy array

def getFlux(HIPvsSH,filename, wavelength,co,interpolated_dustMap,cwd):
    pdtable = pd.read_pickle(filename)
    
    pdO = pdtable.loc[pdtable['teff50']>29000]
    pdnotO = pdtable.loc[pdtable['teff50']<=29000]
    pdB = pdnotO.loc[pdnotO['teff50']>10500]
    pdA = pdnotO.loc[pdnotO['teff50']<=10500]
    F_O = compute(HIPvsSH,pdO, wavelength,co,interpolated_dustMap,cwd)
    F_B = compute(HIPvsSH,pdB, wavelength,co,interpolated_dustMap,cwd)
    F_A = compute(HIPvsSH,pdA, wavelength,co,interpolated_dustMap,cwd)
    return F_O, F_B, F_A
    
