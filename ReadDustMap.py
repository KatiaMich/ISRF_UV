import numpy as np
import pandas as pd

def averagedExtinction(ListOfVects,co,interpolated_dustMap):
    numberPoints = 300 
    dl = np.linspace(0,1,numberPoints)**2 #quadratic spacing of the sample points
    VectOfSight = ListOfVects- co[None,:] #co are OriginPoints
    norms = np.linalg.norm(VectOfSight,axis=1)
    extinctions = np.zeros((numberPoints,len(VectOfSight)),dtype=float)
    averagedExtinctions = np.zeros(len(VectOfSight),dtype=float)
    for k in range(numberPoints):
        LineOfSight = dl[k]*VectOfSight + co[None,:] 
        extinctions[k] = interpolated_dustMap(LineOfSight) #extinctions[k] has as many elements as stars in the catalog
        averagedExtinctions += extinctions[k]*dl[k]/np.sum(dl) 
    return averagedExtinctions,norms

def ListOfVect(pdtable): #gives the coordinates of the star in the frame of the dustmap, which is galactic, centered at the sun.
    x,y,z = pdtable['x(pc)'], pdtable['y(pc)'], pdtable['z(pc)']
    x,y,z = np.array(x),np.array(y),np.array(z) 
    Coord = np.array([x,y,z])
    return np.transpose(Coord)

def getAv(filename,co,interpolated_dustMap): 
    pdtable = pd.read_pickle(filename)
    ListOfVects = ListOfVect(pdtable)
    dA_dD, norms = averagedExtinction(ListOfVects,co,interpolated_dustMap)
    Av = np.array(dA_dD * norms,dtype=float)
    return Av , norms 
