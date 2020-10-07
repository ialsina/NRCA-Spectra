import os
import sys
from tqdm import tqdm
from datetime import datetime
import time
import shutil

import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as plt

import spectra_Basics as basic
from spectra_InitSettings import cf, err

def Interpolate(arrayin,value):
    """Function that locates the value between to elements in an array and returns the linear interpolation."""
    arrx = arrayin[0]
    arry = arrayin[1]
    if value in arrx: return arry[basic.GetIndex(arrx,value)]
    x0,x1 = basic.InBetween(arrx,value,False)
    i0,i1 = basic.InBetween(arrx,value,True)
    coef = (value-x0)/(x1-x0)
    return arry[i1]*coef + arry[i0]*(1-coef)


def GetWeighted(Dict,element,suf,composition):
    """Function that weighs elements from a dictionary according to its composition.
    inputs:
        - Dict: dictionary to look up ingredients (cupboard).
        - element: name of the new mix
        - suf: suffix. Should be '_n-tot' or '_n-g'
        - composition: dictionary of abundances to weight the components (cookbook)."""
    from spectra_Objects import Isotope, Element, Compound
    components = dict()

    #The new mesh (x values) goes from the maximum x value amongst the first x values
    #to the minimum x value amongst the last x values.
    firsts = [Dict[c+suf].spectrum[0,0] for c in composition]
    lasts = [Dict[c+suf].spectrum[0,-1] for c in composition]
    startmesh,endmesh = max(firsts),min(lasts)

    #We start it and compute the mathematical union of the individual meshes. But only the bits that fall within startmesh and endmesh
    mesh = np.empty((0))
    for c in composition:
            mesh = np.union1d(mesh, Dict[c+suf].spectrum[0, basic.InBetween(Dict[c+suf].spectrum[0,:],startmesh,True)[1] : basic.InBetween(Dict[c+suf].spectrum[0,:],endmesh,True)[0] ])

    #The y values are computed by interpolating each of the y values in the mesh and stored in the components dictionary.
    for component in composition:
        time.sleep(.5)
        valy = np.empty((0))
        for xx in mesh:#tqdm(mesh,desc=component,leave=False):
            valy = np.append(valy,Interpolate(Dict[component+suf].spectrum,xx))
        components[component] = valy
    
    #We join together the mesh and the sum of the components weighted by their correspoinding value, and return it
    return np.array([mesh,np.sum([components[comp]*composition[comp] for comp in composition],axis=0)])