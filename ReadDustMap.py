import numpy as np
import pandas as pd

def averagedExtinction(ListOfVects,co,interpolated_dustMap): #gives a list of extinction between the target location and all stars in the catalogue
    numberPoints = 300 
    dl = np.linspace(0,1,numberPoints)**2 #quadratic spacing of the sample points
    VectOfSight = ListOfVects- co[None,:] #co are OriginPoints
    norms = np.linalg.norm(VectOfSight,axis=1)
    extinctions = np.zeros((numberPoints,len(VectOfSight)),dtype=float)
    differential_averagedExtinctions = np.zeros(len(VectOfSight),dtype=float)
    for k in range(numberPoints):
        point = dl[k]*VectOfSight + co[None,:] 
        differential_extinctions[k] = interpolated_dustMap(point) #extinctions[k] has as many elements as stars in the catalog
        differential_averagedExtinctions += differential_extinctions[k]*dl[k]/np.sum(dl) 
    averaged_extinctions = np.array(differential_averagedExtinctions*norms, dtype=floats)
    return differential_averagedExtinctions,norms

def ListOfVect(catalogue_filename): #gives the coordinates of the star in the frame of the dustmap, which is galactic, centered at the sun.
    catalogue = pd.read_pickle(catalogue_filename)
    x,y,z = catalogue['x(pc)'], catalogue['y(pc)'], catalogue['z(pc)']
    x,y,z = np.array(x),np.array(y),np.array(z) 
    Coord = np.array([x,y,z])
    return np.transpose(Coord)

def getAv(catalogue_filename,co,interpolated_dustMap): 
    ListOfVects = ListOfVect(catalogue_filename)
    Av, norms = averagedExtinction(ListOfVects,co,interpolated_dustMap)
    return Av , norms 
