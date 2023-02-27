import pandas as pd
import numpy as np
import scipy as sc
from scipy import interpolate

import ReadDustMap as dust

def colorFunction(g,Teff,wavelength): 
    X = np.array(g['teff50'])
    Y = np.array(g['ratio_'+str(wavelength)+'A'])
    function = sc.interpolate.interp1d(X,Y, fill_value = 'extrapolate')
    return function(Teff)

def extinctionRatio(wavelength,extinction_Function):
    return float(extinction_Function.loc[extinction_Function['wavelength']==wavelength]['Alambda_over_A_V'])
    
def Vegaflux(HIPvsSH):
    if HIPvsSH == 'SH':
        #flux in the gaia band pass, svo
        return 2.50386e-9
    if HIPvsSH == 'HIP': 
        #flux in the  hipparcos band pass, svo
        return 3.76938e-9 
    
def exponent(HIPvsSH,filename,co,interpolated_dustMap,Alambda_overAv):
    pdtable = pd.read_pickle(filename)
    Av , norms = dust.getAv(filename,co,interpolated_dustMap) #dust Extinction
    
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

def getFlux(HIPvsSH,filename, wavelength,co,interpolated_dustMap,cwd):
    pdtable = pd.read_pickle(filename)
    extinction_Function, g, Teff = importFiles(HIPvsSH,pdtable,cwd)   

    Alambda_overAv = extinctionRatio(wavelength,extinction_Function)
    expon = exponent(HIPvsSH, filename,co,interpolated_dustMap,Alambda_overAv)

    g_lambda = colorFunction(g,Teff,wavelength)
    F_Vega_detect = Vegaflux(HIPvsSH) 
    F_lambda = F_Vega_detect * g_lambda * np.exp( np.log(10) * expon )
    return F_lambda #this is a numpy array