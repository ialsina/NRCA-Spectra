import os
import sys
from tqdm import tqdm
from datetime import datetime
import time
import shutil

import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as plt

from spectra_InitSettings import cf,err

# HANDY VARIABLES
path = os.getcwd()
isd = lambda d: os.path.isdir(os.path.join(path,d))
isfx= lambda f: os.path.isfile(os.path.join(path,f))
isf = lambda f: isfx(f+'.txt')


catalog_volumes = ('isotopes','elements','compounds','samples')


#Take functions defined in a tricky place for them to be used whenever
L0, E2t, t2E, dt2dE, dE2dt = cf._funcs()


def InterpretName(name):
    """A function that takes a string with the typical format we are using and splits its components.
    The string format should be: %d-%s-%d_%s (where %d stands for integer and %s for string).
    inputs:
        - name: string to be interpreted
    outputs:
        - znumber: atomic number
        - element: chemical symbol
        - nnumber: atomic mass
        - mode: should be either 'n-tot' or 'n-g'
    Warning: some tricky inputs are considered, but not comprehensively. Avoid entering strings
    that aren't in the proper format"""
    if name.count('_') > 1 or name.count('-') > 3: return '',name,'',''
    iso_name = name.split('_')[0]
    mode = name.split('_')[1] if name.count('_') > 0 else ''
    znumber = iso_name.split('-')[0] if iso_name.count('-') > 0 else ''
    element = iso_name.split('-')[1] if iso_name.count('-') > 0 else iso_name.split('-')[0]
    nnumber = iso_name.split('-')[2] if iso_name.count('-') > 1 else ''
    return znumber, element, nnumber, mode


def IndMaxima(arr,s=1):
    """Takes an array and looks for the indices that correspond either to local maxima or to local minima
    inputs:
        - arr: numpy array to consider
        - s:    1: means local maxima
               -1: means local minima
    outputs:
        - iextr: array of local extrema (maxima or minima) indices"""
    dersgn = np.zeros((len(arr)-1))
    iextr = np.empty((0),dtype=np.int64)
    for i in range(len(arr)-1):
        dersgn[i] = np.sign(arr[i+1]-arr[i])
    for i in range(0,len(arr)-2):
        if dersgn[i] == s and dersgn[i+1] == -s:
            iextr = np.append(iextr,i+1)
        if dersgn[i] == s and dersgn[i+1] == 0:
            for j in range(i+1, len(arr)-1):
                if dersgn[j] != 0:
                    if dersgn[j] == -s:
                        iextr = np.append(iextr,np.arange(i+1,j+1))
                    break
    return iextr

# A couple little functions to ask for parameters. Just for code reusability's sake.
def AskAxis():
    return True if input('x-axis: (1 eV; [2] ToF) >').lower() in ['1','ev'] else False
def AskLim():
    return True if input('Show detection limits? ([y]/n) >').lower() not in ['n','no'] else False

def GetIndex(arr,elements,single_as_int=True):
    """Given an array and an element of it, gives back the index array where the element is found.
    Works for ordered iterables, too. Then it gives back an array of indices.
    inputs:
        - arr (np.ndarray):
            numpy array to be searched in
        - elements (int, float, list, tuple, np.ndarray):
            single element or list/tuple/array of elements to search in arr
        - single_as_int (bool):
            True: if the output is a single value, it is given as int
            False: if the output is a single value, it is given as an array of shape (1,)
    outputs:
        - outp (int/array)
            index/indices of elements in arr
    Warnings:
        - if the element occurs more than once in the array, the first indexs is returned.
        - if the element doesn't belong to the array, returns None"""
    outp = None
    if type(elements) in [int,float,np.float_]: elements = np.array([elements])
    if type(elements) in [list,tuple]: elements = np.array(elements)
    outp = np.array([np.nonzero(arr==el)[0][0] for el in elements])
    if not outp is None and np.size(outp) == 1 and single_as_int: outp = outp[0]
    return outp


def InBetween(arr,val,outp_i=False):
    """Given an ordered array and a value, returns the elements in between of which the value would stand.
    inputs:
        - arr (np.ndarray):
            numpy array to be searched in
        - val (int, float):
            single element to 'sandwich' in arr
        - outp_i (bool):
            True: the indices of the elements that 'sandwich' the input value are returned.
            False: the actual elements that 'sandwich' the input value are returned
    outputs:
        - (2,)-shaped ndarray() with either the array elements that 'sandwich' the input value or
            their indices.
    Warning:
        - if the element is below the minimum value of the array or above the maximum, returns a
            (2,)-shaped ndarray of None.
            """
    if val in arr: return np.array([val,val]) if not outp_i else GetIndex(arr,np.array([val,val]))
    for i in range(len(arr)-1):
        if val>arr[i] and val<arr[i+1] or val<arr[i] and val>arr[i+1]:
            return np.array([arr[i],arr[i+1]]) if not outp_i else np.array([i,i+1])
    return np.array([None,None])

def Closest(arr,val,outp_i=False):
    """Given an ordered array and a value, returns the closest element of the array to the value.
    inputs:
        - arr (np.ndarray):
            numpy array to be searched in
        - val (int, float):
            single element to compare in arr
        - outp_i (bool):
            True: the index of the closest element in arr is returned.
            False: the actual closest element in arr is returned.
    outputs:
        - either element in arr which is closest to val or its index.
    Warning:
        - if the element is below the minimum value of the array or above the maximum, returns None"""
    sides = InBetween(arr,val,False)
    if sides[0] is None: return None
    res = sides[np.argmin(np.abs(sides-val))]
    return res if not outp_i else GetIndex(arr,res)

def Chunk(arr,tup):
    """Given an ordered array and a tuple of values, gives the array that is best bounded by the
    values in the tuple.
    inputs:
        - arr (np.ndarray):
            numpy array to be searched in. It can have one or two axes. If it has two, the bounds are
            selected from the 0-th row.
        - tup (dim-2 tuple):
            tuple of lower and upper limits aiming to bound the desired array.
    outputs:
        - np.array with the same number of dimensions as the input one, bounded by elements in tup."""
    layers = np.shape(arr)[0] > 1
    warr = arr[0] if layers else arr[:]
    imin = Closest(warr,tup[0],1)
    imax = Closest(warr,tup[1],1)
    return np.vstack(warr[:,imin:imax+1]) if layers else warr[imin:imax+1]

def Smooth(arr,it):
    """Smooths an array. If it has two axes, (x being 0th and y being 1st row), y axes is smoothed.
    The smoothing mesh can be outlined as {1 2 1}/4
    inputs:
        - arr: (np.ndarray):
            numpy array to smooth.
        -it: (int)
            number of iterations
    outputs:
        - smoothed array"""
    if len(np.shape(arr)) > 1:
        arr0 = np.copy(arr)[1]
    else:
        arr0 = np.copy(arr)
    for j in range(it):
        temp = [0.]
        for i in range(1,np.size(arr0)-1):
            temp.append((arr0[i-1]+2*arr0[i]+arr0[i+1])/4)
        temp.append(0.)
        arr0 = np.copy(np.array(temp))
    if len(np.shape(arr)) > 1 and np.shape(arr)[0] > 1:
        arr0 = [arr[0],arr0]
        for i in range(2,np.shape(arr)[0]):
            arr0.append(arr[i])
    return np.array(arr0)

def FitBoxes(array,dboxes):
    """Forces all the y values of an array to be fit within a made-up mesh of a certain amount of points.
    This amount of points is actually given by as a density, so that num_points = (y_max - y_min)*density
    inputs:
        -array: (np.ndarray):
            numpy array to fit
        - dboxes:
            density of points/boxes
    outputs:
        -fitted data as an array of the same shape, box width"""
    try:
        b0 = np.min(array)
        b1 = np.max(array)
        nboxes = dboxes*(b1-b0)
        norm = (array-b0)/(b1-b0)
        fitn = np.around(norm*(nboxes-1))/(nboxes-1)
        fit = (b1-b0)*fitn+b0
        return fit, (b1-b0)/(2*nboxes)
    except Exception as e:
        err.add('FitBoxes',isotope,e)
        return None, None