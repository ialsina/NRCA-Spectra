#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:48:50 2020

@author: ivan
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint

import spectra_Basics as basic
import spectra_Plotters as plotter
import spectra_Finders as finder
from spectra_InitSettings import cf, peakattr, err

def definepeak(self,npeak,prange,params):
    """Given an Data instance, defines the npeak-th peak in x position.
    How it works:
        - prange defines the neighbourhood we are interested in.
        - the derivative in this nhood is considered.
        - it finds the outerslope, i.e. the slope far away from the peak.
            This is done by making a histogram with FitBoxes and looking at the highest occurrence.
            If the outerslope is higher in abs value to a certain amount (maxouterslope), it is set to 0
        - For each direction, sign -1 (+1) for left (right):
            Starting from the center, we set the derivative to 0 and we go left (right).
            When we start decreasing (increasing) in derivative, we trigger a 'lock' variable and start counting.
            At this point, we record the derivative (dermax) and check its value at each step further, as it
            approaches the outerslope value.
            Once the derivative has dropped a certain amount (slopedrop) between dermax and outerslope, we set the limit.
            At any moment, if the derivative changes sign, we set a limit as well.
            If the two limits are set because the derivative doesn't hit the end of the range, the peak is successfully bounded.
    inputs:
        - self: class
            an instance of Isotope, Element or Compound.
        - npeak: int
            ranked peak position (x).
        - prange: int
            range in number of points to the left and to the right. Total range should be 2*prange + 1 (1 for the central value, peak)
        - params: dictionary that has tro bring information in:
            -dboxes
            -maxouterslope
            -slopedrop
    outputs:
        - outp_e: tuple: (float, float)
            x boundary values
        - outp_i: tuple: (int, int)
            x boundary indices
        -         tuple: (float, float)
            outerslope and box width (uncertainty, in a way)
        - peakr: tuple: (int, int)
            for each boundary: 1: stopped because the derivative dropped enough
                               2: stopped because the derivative changed sign. This is set as an 'emergency stopping'
                                    and happens in very fast changes in slope.
                               3: user-defined (not implemented here, but nice to know)

    Warning: In case of error, when one or both of the peak boundaries aren't defined, same thing is returned,
                only filled with zeros."""
    try:
        #Center and redder (reduced derivative), i.e. smoothed derivative in a neighbourhood defined by prange
        center = self.mai[npeak]
        redder = self.sder[center-prange:center+prange+1]
        for i in range(1,10):
            fit,boxwidth = basic.FitBoxes(redder,params['dboxes']*i)
            temp1,temp2 = np.unique(fit,return_counts=True)
            if not (temp2==1).all(): break
        outerslope = temp1[np.argmax(temp2)]
        if abs(outerslope) > params['maxouterslope']: outerslope = 0
        outp_i = []
        outp_e = []
        peakr = []

        for sign in [-1,1]:
            decreasing, increasing, lock = False, False, False
            dermax = None
            for i in range(center,center+sign*(prange+1),sign):
                decreasing = self.sder[i+sign]<self.sder[i]
                increasing = self.sder[i+sign]>self.sder[i]
                if not lock:
                    if (decreasing if sign==-1 else increasing):
                        lock = True
                        dermax = self.sder[i]
                        if dermax == outerslope: raise Exception('Non-standing slope',prange)
                if lock:
                    if abs((self.der[i]-outerslope)/(dermax-outerslope)) <= params['slopedrop']:
                        peakr.append(1)
                        outp_i.append(i)
                        outp_e.append(self.spectrum[0,i])
                        break
                    if self.sder[i+sign]*self.sder[i] <= 0:
                        peakr.append(2)
                        outp_i.append(i)
                        outp_e.append(self.spectrum[0,i])
                        break

        if np.size(outp_e)!=2 or np.size(outp_i)!=2: raise Exception('Unable to Define Peak',prange)
        if outp_e[0] == outp_e[1]: raise Exception('Zero-width peak',prange)
        return tuple(outp_e), tuple(outp_i), tuple([outerslope,boxwidth]), tuple(peakr)
    except Exception as e:
        err.add(e,'definepeak',self.fullname,npeak)
        return (0.,0.), (0,0), (0.,0.), (0,0)

def Integrate(array,iedges):
    """Integrates an x-y array between two given indices with simpson, and subtracting the background
    assuming a trapezoid.
    inputs:
        - array: np.ndarray
            array to be integrated
        - iedges: tuple, (ind, ind)
            edges between which we want the integration.
    outputs:
        Integral minus background."""
    listx = array[0,iedges[0]:iedges[1]+1]
    listy = array[1,iedges[0]:iedges[1]+1]
    raw_int = spint.simps(listy,listx)
    background = 0.5*(array[0,iedges[1]]-array[0,iedges[0]])*(array[1,iedges[1]]+array[1,iedges[0]])
    return raw_int-background

def Fwhm(array,ic,coords,ilims):
    """Computes the Full Width at Half of Maximum (FWHM) of a peak.
    inputs:
        - array:
        - ic: int
            peak center (index)
        - coords: tuple (float, float)
            peak coordinates
        - ilims:
            peak limits
    outputs:
        - fwhm: float
            FWHM. Zero in case of miscalculation."""
    yhm = coords[1]/2
    i0,i1 = ilims[0],ilims[1]
    iredhm1 = basic.InBetween(array[1,i0:ic],yhm,True)[0]
    iredhm2 = basic.InBetween(array[1,ic:i1+1],yhm,True)[1]
    xhm1 = array[0,iredhm1+i0] if not iredhm1 is None else None
    xhm2 = array[0,iredhm2+ic] if not iredhm2 is None else None
    return xhm2 - xhm1 if (( not xhm1 is None ) and (not xhm2 is None)) else 0

def computepeak(self,peakpos,prange,params,setx=None):
    """Computes every peak parameter, sets up the tuple of them, and creates the
    actual peak instance. For more info on this parameters go to the peakattr documentation.
    inputs:
        - self: object instance
        - peakpos: rank in peak position
        - prange: prange
        - params: parameters
        - setx: tuple with the peak boundaries. Used for the peak editing tool.
    outputs:
        - Peak instance."""
        
    from spectra_Objects import Peak
    
    try:
        center = self.ma[0,peakpos]
        coords = tuple(self.ma[:,peakpos])
        coords_tof = tuple((basic.E2t(coords[0]), coords[1]))
        icenter= self.mai[peakpos]
        center_tof = basic.E2t(center,self.mode)

        # If setx is None, we are computing. Otherwise, we com from the editing function.
        if setx is None:
            xlims,ilims,outerslope,peakreason = definepeak(self,peakpos,prange,params)
            user_edited = False
        else:
            xlims = tuple(setx)
            ilims = tuple(basic.GetIndex(self.spectrum[0], setx))
            outerslope = (0., 0.)
            peakreason = (3, 3)
            user_edited = True
        successful = True if xlims != (0.,0.) else False
        integral = Integrate(self.spectrum,ilims) if successful else 0
        #integral_tof0 = integral*dE2dt(1,center,self.mode) if successful else 0
        integral_tof = -Integrate(self.spectrum_tof,ilims) if successful else 0
        yvals = tuple(self.spectrum[1,ilims]) if successful else 0
        width = xlims[1]-xlims[0] if successful else -1
        height = self.ma[1,peakpos]-(yvals[1]+yvals[0])/2 if successful else -1
        fwhm = Fwhm(self.spectrum,icenter,coords,ilims) if successful else 0
        #Keep calm: big tuple.
        #Rank entries are set to -1 to show that are still unknown.
        #0 values are stored for computation problems, e.g. properties of unbounded peak.
        info = tuple((
                self.fullname,            #00 #fullname
                -1,                         #01 #num
                peakpos,                    #02 #center_
                center,                     #03 #center
                icenter,                    #04 #icenter
                center_tof,                 #05 #center_tof
                integral,                   #06 #integral
                -1,                         #07 #integral_
                integral_tof,               #08 #integral_tof
                width,                      #09 #width
                -1,                         #10 #width_
                height,                     #11 #height
                -1,                         #12 #height_
                fwhm,                       #13 #fwhm
                -1,                         #14 #fwhm_
                integral/(height**2),       #15 #ahh
                -1,                         #16 #ahh_
                integral/(height*width),    #17 #ahw
                -1,                         #18 #hw_
                xlims,                      #19 #xlims
                ilims,                      #20 #ilims
                outerslope,                 #21 #outerslope
                peakreason,                 #22 #peakreason
                yvals,                      #23 #yvals
                coords,                     #24 #coords
                coords_tof,                 #25 #coords_tof
                prange,                     #26 #prange
                user_edited,                #27 #user_edited
                False,                      #28 #user defined     
                ))
        return Peak(info)
    except Exception as e:
        err.add(e,'computepeak',self.fullname,-1)
        return Peak(info)


def maxima(arr,xbounds=None,ybounds=None,smoothing=0):
    """Looks for local maxima in the array, restricted to the bounds set.
    input:
        - arr:
            Numpy array to search
        - xboudns:
            If None, there aren't. If tuple, peaks are restricted to those values.
        - ybounds:
            If None, there isn't. If float, peaks are restricted to be above this value.
        - smoothing: int
            Pre-smoothing of the y values.
    output:
        - cmaxima: array (horizontal) with the x-y values of the maxima
        - imaxima: array 1-d with the peak indices of the maxima"""

    arr = basic.Smooth(arr,smoothing)

    # First we look for maxima candidates
    imax = basic.IndMaxima(arr[1,:])
    cmaxima = arr[:,imax]

    #Array of booleans. Do the peaks in cmaxima accomplish the condition for x and y?
    condx = np.array(cmaxima[0]>=xbounds[0]) * np.array(cmaxima[0]<=xbounds[1]) if not xbounds is None else np.ones(np.shape(cmaxima)[1])
    condy = cmaxima[1]>=ybounds if not ybounds is None else np.ones(np.shape(cmaxima)[1])

    # We peak the ones that do from the maxima candidates array.
    cmaxima = cmaxima[:,np.nonzero(condx*condy)[0]]
    imaxima = basic.GetIndex(arr[0],cmaxima[0],False)
    return cmaxima,imaxima

def minima(arr,xbounds=None,ybounds=None,smoothing=0):
    """Looks for local minima in the array, restricted to the boudns set.
    Exactly the same as maxima(). Look for that documentation."""
    imin = basic.IndMaxima(arr[1,:], -1)
    cminima = arr[:,imin]
    condx = np.array(cminima[0]>=xbounds[0]) * np.array(cminima[0]<=xbounds[1]) if not xbounds is None else np.ones(np.shape(cminima)[1])
    condy = cminima[1]<=ybounds if not ybounds is None else np.ones(np.shape(cminima)[1])
    cminima = cminima[:,np.nonzero(condx*condy)[0]]
    iminima = basic.GetIndex(arr[0],cminima[0],False)
    return cminima,iminima

def propsisot(self,params,setx=dict()):
    """Computes all the peaks for a given instance of Data.
    inputs:
        - self: instance
        - params: dictionary with needed settings.
        - setx:
            dictionary of peak position ranks as keys and tuples with bounds as values.
            Used when editing peaks. For the peaks that have a value already set, it is
            passed to computepeak()
    outputs:
        - sorted dictionary of peaks."""
    peaks_pos = {}
    for i in range(np.size(self.mai)):
        try:
            nleft = self.mai[i]
            nright = np.shape(self.spectrum)[1] - self.mai[i] - 1
            if np.size(self.mai) > 1:
                if i == 0:
                    prange = min(params['prangemax'], nleft, nright)
                elif i == np.size(self.mai)-1:
                    prange = min(params['prangemax'], nleft, nright)
                else:
                    prange = min(params['prangemax'],self.mai[i+1]-self.mai[i-1],nleft,nright)
            elif np.size(self.mai) == 1:
                prange = min(self.mai[0]-1,np.shape(self.spectrum[1])-self.mai[-1]-1,params['prangemax'])
            else:
                prange = 0
            peaks_pos[i] = computepeak(self,i,prange,params,setx.get(i))
        except Exception as e:
            err.add(e,'propsisot',self.fullname,i)
            continue

    return sorting(peaks_pos)
        
def sorting(inp):
    """Function that ranks the rankable parameters, and orders the peaks so that they become label-ranked by intensity.
    input:
        -inp: dictionary of unsorted peaks.
    output:
        - dictionary of sorted peaks with rank parameters set and ready."""

    peaks = list(inp.values())

    #Build and array where each column correspond to an instance and each row to a property.
    array = np.empty((6,0))
    for el in peaks:
        array = np.append(array, np.array([el.integral, el.width, el.height, el.fwhm, el.ahh, el.ahw]).reshape(6,1),axis=1)
    
    #Index i points at the position in peaks whose rank is the #i
    integral_sorting = (-array).argsort()[0]
    #For each row (parameter): index i gives the #rank that peak position i has.
    ranks = (-array).argsort().argsort()

#   Alternative (probably faster) method (problem: within same integral val., order is backwards in peak position):
#   integral_sorting = np.flip(array.argsort(),axis=1)[0]
#   ranks = np.flip(array.argsort(),axis=1).argsort()
    
    for el in range(len(peaks)):
        peaks[el].integral_ = ranks[0,el]
        peaks[el].width_ = ranks[1,el]
        peaks[el].height_ = ranks[2,el]
        peaks[el].fwhm_ = ranks[3,el]
        peaks[el].ahh_ = ranks[4,el]
        peaks[el].ahw_ = ranks[5,el]
        peaks[el].num = peaks[el].integral_
    
    #Build and return the new dictionary of peaks.
    return {i: peaks[integral_sorting[i]] for i in range(len(peaks))}

def sampprocess(self):
    """It is a method.
    Strips peaks from the sample spectra by joining the minima a certain number of iterations.
    Stores convenient parameters as well."""
    try:
        #Omit the points with positive slope at the outermost left.
        for i in range(np.shape(self.spectrum_tof)[1]):
            if self.der[i] < 0: break
        processed = self.spectrum_tof[:,i:]
        processed = basic.Smooth(processed, 1)[:,1:]
        
        #Get minima
        self.mi_tof, self.mii_tof = minima(processed)
        self.stripped = np.copy(processed)

        #Iterate stripping over the minima
        for it in range(cf.iterspeaks):
            for i in range(np.shape(processed)[1]):
                processed = np.copy(self.stripped)
                _, tempmii = minima(processed)
                iends = basic.InBetween(tempmii, i, False)
                if iends[0] == 0: continue
                if iends[0] == iends[1]: continue
                if iends[0] is None or iends[1] is None: continue
                x0, y0 = processed[0:2,iends[0]]
                x1, y1 = processed[0:2,iends[1]]
                A = (y1-y0)/(x1-x0)
                B = (y0*x1-y1*x0)/(x1-x0)
                self.stripped[1,i] = A*processed[0,i] + B

        self.stripped = basic.Smooth(self.stripped, 1)[:,1:]

        #Coefficients of the chosen degree in the settings file. Biggest order first.
        #i.e. [A, B, C, D] means A*x**3+B*x**2+C*x+D
        #Also fitted background.
        self.coeffs = np.polyfit(self.stripped[0],self.stripped[1],cf.fitting_coeff-1)
        self.background_tof = np.vstack((self.spectrum_tof[0], np.poly1d(self.coeffs)(self.spectrum_tof[0])))
    except Exception as e:
        err.add(e,'sampprocess',self.fullname,-1)


def RankNearest(Dict,xx):
    closest_peak = np.empty((0,3))
    for isotname in Dict:
        isot = Dict[isotname]
        dist = np.abs(isot.get_from_peaks('center_tof')-xx)
        if np.size(dist) == 0: continue
        closest_dist = np.min(dist)
        nisot = np.argwhere(dist == closest_dist)[0,0]
        closest_peak = np.append(closest_peak, np.array([isot, closest_dist, nisot]).reshape(1,-1), axis=0)
    closest_isot = closest_peak[np.argsort(closest_peak[:,1])]
    return closest_isot


def MatchPeaks(self,distmax=cf.max_match,samp=None):
    try:
        if not samp: samp = self.get(finder.Select(self.get_as_dict(isotopes=False, elements=False, compounds=False), recursive=False, ask_if_one=False))

        plotter.PlotBars(samp,dict())

        while True:
            ind = AskPeak(np.size(samp.ma_tof[0]))
            if ind == -2: break
            xx = samp.ma_tof[0,ind]
            closest_isot = RankNearest(self.Datas(samp.mode), xx)
            closest_isot = closest_isot[np.where(closest_isot[:,1]<distmax)[0]]
            for i in range(len(closest_isot)):
                print('{:3d}: {:>20s} ({:>4d} - {:6.3f})'.format(i, closest_isot[i][0].fullname, closest_isot[i][2], closest_isot[i][1]))
    except Exception as e:
        err.add(e,'MatchPeaks','',-1)


def AskPeak(firstno):
    p = input('Enter peak: >')
    if p in ['q', 'quit', '']: return -2
    if not p.isnumeric():
        print('Invalid input.')
        return -1
    ind = int(p)
    if ind not in np.arange(0,firstno):
        print('Invalid input.')
        return -1
    return int(p)


def ComparePeaks(self):
    from spectra_Objects import Summer
    try:
        samp, Dict = self.smart_select()
        plotter.PlotBars(samp,Dict)
        intratios = Summer()
        while True:
            ind = AskPeak(np.size(samp.ma_tof[0]))
            if ind == -2: break
            if ind == -1: continue
            inp = input('\tEnter peak boundaries: >')
            if not inp.replace(',','').replace('.','').isnumeric():
                print('Invalid input')
                continue
            inp = inp.split(',')
            try:
                newlims_ = (float(inp[0]), float(inp[1]))
            except:
                print('Invalid input')
                continue
            xx = samp.ma_tof[0,ind]
            closest_isot = RankNearest(Dict, xx)
            newlims = (basic.Closest(samp.spectrum_tof[0], newlims_[0], True), basic.Closest(samp.spectrum_tof[0], newlims_[1], True))
            print('\tPeak center: {}. Peak newlims: {}'.format(xx, newlims))
            integral = Integrate(samp.spectrum_tof, newlims) - Integrate(samp.stripped, newlims)
            print('\tClosest component peaks:')
            for i in range(len(closest_isot)):
                print('\t{}: {} ({})'.format(i, closest_isot[i][0].fullname, closest_isot[i][1]))
            inp2 = input('\tEnter one of the above components: >')
            if inp2.isnumeric():
                inp2 = int(inp2)
                if inp2 < len(closest_isot):
                    ichoose = closest_isot[inp2,0]
                    pchoose = closest_isot[inp2,2]
                else:
                    ichoose = closest_isot[0,0]
                    pchoose = closest_isot[0,2]

            else:
                ichoose = closest_isot[0,0]
                pchoose = closest_isot[0,2]

            intratio = integral / ichoose.peaks[pchoose].integral_tof
            intratios.append(ichoose.fullname, intratio)
            print('Closest component: {}. Intensities ratio: {:3f}'.format(ichoose.fullname,intratio))

        print('Composition results for {}:'.format(samp.fullname))
        outcome = intratios.percentage(True,2)
        for el in outcome:
            print('{:>12s}: {:>5.2f}%'.format(el, outcome[el]))

        return intratios
    except Exception as e:
        err.add(e,'ComparePeaks','',-1)


def EditPeaks(self):
    try:
        dictout = {}
        listnum = []
        through = False
        self.plot()
        while True:
            ind = AskPeak(len(self.peaks))
            if ind == -2:
                break
            if ind == -1:
                continue
            if through: plt.close()
            self.plotsingle(ind)
            print('\tOld peak boundaries: {}'.format(self.peaks[ind].xlims))
            inp = input('\tEnter new peak boundaries: >').replace(' ','')
            if not inp.replace(',','').replace('.','').isnumeric():
                print('Invalid input')
                continue
            inp = inp.split(',')
            try:
                newlims_ = (float(inp[0]), float(inp[1]))
            except:
                print('Invalid input')
                continue
            newlims = (basic.Closest(self.spectrum[0], newlims_[0], False), basic.Closest(self.spectrum[0], newlims_[1], False))
            print('\tNew peak boundaries: {}'.format(newlims))
            plt.axvline(x=newlims[0], color='green')
            plt.axvline(x=newlims[1], color='green')
            plt.show()
            dictout[self.peaks[ind].center_] = newlims
            listnum.append(ind)
            through = True
        print('This action will edit the following peaks:\n',listnum)
        if not input('Continue? ([y]/n) >') in ['n','no']:
            plt.close('all')
            return dictout
        else:
            plt.close()
            return dict()

    except Exception as e:
        err.add(e,'EditPeaks',self.fullname,-1)


def DeletePeaks(self):
    try:
        listnum = []
        self.plot()
        while True:
            ind = AskPeak(len(self.peaks))
            if ind == -2: break
            if ind == -1: continue
            listnum.append(ind)
            print('To delete:', listnum)
        print('This action will delete the following peaks:\n', listnum)
        if input('Continue? (y/[n]) >') in ['y','yes']:
            plt.close()
            return listnum
        else:
            plt.close()
            return []
    except Exception as e:
        err.add(e,'DeletePeaks',self.fullname,-1)


def pick(self,attr,magn):
        if attr in self.__dict__ and attr+'_tof' in self.__dict__:
            if magn in [1,True,'ToF','tof','time of flight']:
                return getattr(self,attr+'_tof')
            else:
                return getattr(self,attr)
        else:
            if attr in self.__dict__:
                return getattr(self,attr)
            elif attr+'_tof' in self.__dict__:
                return getattr(self,attr+'_tof')
            else:
                return None