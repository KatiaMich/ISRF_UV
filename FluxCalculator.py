import pandas as pd
import numpy as np
import scipy as sc
from scipy import interpolate

import ReadDustMap as dust

def colorFunction(g,Teff,wavelength): #the amount of flux in the UV in units of visible flux varies with temperature. The function describing this variation is called colour function.
    temperature = np.array(g['teff50'])
    flux_ratio = np.array(g['ratio_'+str(wavelength)+'A'])
    color_function = sc.interpolate.interp1d(temperature,flux_ratio, fill_value = 'extrapolate')
    return color_function(Teff)

def extinctionRatio(wavelength,extinction_Function): #the extinction depends on wavelength. Since extinction by dust is provided in the visual band, in it convenient to work with the ratio.
    ext = extinction_Function.loc[extinction_Function['wavelength']==wavelength]['Alambda_over_A_V']
    return float(ext)
    
def Vegaflux(HIPvsSH):
    if HIPvsSH == 'SH':
        #flux from Vega in the gaia band pass, from svo
        return 2.50386e-9
    if HIPvsSH == 'HIP': 
        #flux from Vega in the hipparcos band pass, from svo
        return 3.76938e-9 
    
def exponent(HIPvsSH,filename,co,interpolated_dustMap,Alambda_overAv):
    pdtable = pd.read_pickle(filename)
    Av , norms = dust.getAv(filename,co,interpolated_dustMap) #dust Extinction
    Av = 1/2*Av
    if HIPvsSH == 'SH':
        absMag = pdtable['mg0'] #in the g band
        expo = np.array( (-0.4)*(absMag - 5 + 5*np.log10(norms) + Alambda_overAv*Av  ) ) 
        
    if HIPvsSH == 'HIP':
        absMag = pdtable['absMag_HIP'] #in the hipparcos band
        expo = np.array( (-0.4)*(absMag - 5 + 5*np.log10(norms) + Alambda_overAv*Av  ) ) 
    
    return expo #this is a numpy array

def importFiles(HIPvsSH,catalogue,cwd):
    extinction_Function = pd.read_pickle(cwd+'/inputs/extinctionFunction')
    Teff = catalogue['teff50']
    
    if HIPvsSH == 'SH':
        g = pd.read_pickle(cwd+'/inputs/ColourFunctionStarHorse.pkl')
    
    if HIPvsSH == 'HIP':
        g = pd.read_pickle(cwd+'/inputs/ColourFunctionHipparcos.pkl')
    return extinction_Function, g, Teff      

def getFlux(HIPvsSH,filename, wavelength,co,interpolated_dustMap,cwd):
    catalogue = pd.read_pickle(filename)
    extinction_Function, g, Teff = importFiles(HIPvsSH,catalogue,cwd)   

    Alambda_overAv = extinctionRatio(wavelength,extinction_Function)
    expon = exponent(HIPvsSH, filename,co,interpolated_dustMap,Alambda_overAv)

    g_lambda = colorFunction(g,Teff,wavelength)
    F_Vega_detect = Vegaflux(HIPvsSH) 
    F_lambda = F_Vega_detect * g_lambda * np.exp( np.log(10) * expon )
    return F_lambda #this is a numpy array
