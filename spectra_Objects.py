import os
import sys
from tqdm import tqdm
from datetime import datetime
import time
import shutil

import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as plt

import spectra_Basics as basics
import spectra_Plotters as plotter
import spectra_Finders as finder
import spectra_ObjectsFunc as func
from spectra_InitSettings import cf, peakattr, err


class Catalog:

    def __init__(self,**kwargs):
        self.loadfiles()
        self.volumes = basics.catalog_volumes
        self.date_created = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    def _format(self, **kwargs):
        for volume in self.volumes:
            setattr(self, volume, kwargs.get(volume,dict()))

    def replace(self,**kwargs):
        for volume in self.volumes:
            setattr(self, volume, kwargs.get(volume,self.volume))

    def update(self,**kwargs):
        for volume in self.volumes:
            setattr(self, volume, dict( getattr(self,volume), **kwargs.get(volume,dict()) ) )

    def loadfiles(self):
        self._format()
        err.start()
        self.data_in()
        self.sample_in()
        self.mix_in()
        err.present()

    def _discriminate(self,Dict,mode=None):
        if not mode:
            return Dict
        else:
            return {el: Dict[el] for el in Dict if Dict[el].mode==mode}

    def Isotopes(self,mode=None):
        return self._discriminate(self.isotopes,mode)

    def Elements(self,mode=None):
        return self._discriminate(self.elements,mode)

    def Compounds(self,mode=None):
        return self._discriminate(self.compounds,mode)

    def Samples(self,mode=None):
        return self._discriminate(self.samples,mode)

    def Substances(self,mode=None):
        return dict(self.Datas(mode), **self.Samples(mode))

    def Datas(self,mode=None):
        return dict(self.Isotopes(mode),**self.Mixes(mode))

    def Mixes(self,mode=None):
        return dict(self.Elements(mode),**self.Compounds(mode))

    def get_as_dict(self,**kwargs):
        return {volume: getattr(self,volume) if kwargs.get(volume,True)==True else dict() for volume in self.volumes}

    def _unravelDatas(self):
        """Returns:
            - isots: list of unique isotopes
            - ielems: list of unique elements built from the list of isotopes
            - elems: list of unique elements
            - non_unique: list of ielems that contain more than one isotope
            - comps: list of compounds
        Note:
            The reason to differentiate between ielems and elems is that this enables to build
            the non_unique list, useful for mixing purposes"""
        
        spl = lambda a, n: '-'.join([a.split('_')[0].split('-')[i] for i in range(n)])
        
        isots = np.unique([spl(isot,3) for isot in self._discriminate(self.isotopes, cf.default_mode)],return_counts=False)
        ielems,nisots = np.unique([spl(isot,2) for isot in self._discriminate(self.isotopes, cf.default_mode)],return_counts=True)
        non_unique = ielems[nisots>1]
        elems = np.unique([spl(elem,2) for elem in self.elements],return_counts=False)
        comps = self.compounds.keys()
        
        return list(isots), list(ielems), list(elems), list(non_unique), list(comps)

    def get_isotopes(self):
        """Isotopes, without suffix."""
        return self._unravelDatas()[0]

    def get_elements(self):
        """Elements, without suffix"""
        return self._unravelDatas()[2]

    def get_elements_from_isotopes(self, non_unique=False):
        """Looking at the isotopes list, list of elements they account for."""
        return self._unravelDatas()[1] if not non_unique else self._unravelDatas()[3]

    def get_compounds(self):
        """Compounds"""
        return self._unravelDatas()[4]

    def get_samples(self):
        return list(self.samples.keys())

    def ready_to_mix(self):
        """Isotopes,element from the DB are always ready to be mixed."""
        outp = self.get_isotopes()
        outp.extend(self.get_elements())
        outp.extend(self.get_compounds())
        return outp

    def find(self,askmode=False,ask_if_one=False):
        return self.Substances().get(finder.Select(self.get_as_dict(),askmode=False,recursive=False,ask_if_one=False))

    def get(self,inp,otherwise=None):
        return self.Substances().get(inp, otherwise)

    def export(self):
        from spectra_FileHandlers import ExportProps, ExportProps2
        ExportProps(self.Datas())
        ExportProps2(self.Datas())

    def save(self):
        from spectra_FileHandlers import psave
        self.date_modified = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        psave(self)

    def data_in(self):
        from spectra_FileHandlers import ImportData
        new_isotopes, new_elements, new_compounds = ImportData()
        self.update(isotopes=new_isotopes, elements=new_elements, compounds=new_compounds)

    def sample_in(self):
        from spectra_FileHandlers import ImportSamp
        self.update(samples=ImportSamp())

    def mix_out(self):
        from spectra_FileHandlers import MixOut
        MixOut(self.Substances(), self._unravelDatas())

    def mix_in(self):
        from spectra_FileHandlers import MixIn
        new_elements, new_compounds = MixIn(self.Datas(), self.ready_to_mix())
        self.update(elements=new_elements, compounds=new_compounds)

    def smart_select(self,samp=None):
        """Makes the user decide one sample and then anything but samples."""
        print('Work on sample:')
        if not samp: samp = self.get(finder.Select(self.get_as_dict(isotopes=False, elements=False, compounds=False), recursive=False, ask_if_one=False))
        print('Select peaks from:')
        isotsl = finder.Select(self.get_as_dict(samples=False), recursive=True, restrict=samp.mode)
        if isotsl == []: return None
        isots = dict()
        for isot in isotsl:
            if self.get(isot).npeaks>1: isots[isot] = self.get(isot)
        return samp, isots

    def plotbars(self):
        samp, Dict = self.smart_select()
        plotter.PlotBars(samp,Dict)

    plot = plotter.Plot
    pmatch = func.MatchPeaks
    pcompare = func.ComparePeaks

class Summer:
    def __init__(self):
        self.values = dict()
        self.counts = dict()
        self.lists = dict()

    def append(self,key,amount):
        if key not in self.values:
            self.values[key] = 0
            self.counts[key] = 0
            self.lists[key] = []
        self.values[key] = (self.values[key]*self.counts[key] + amount)/(self.counts[key]+1)
        self.counts[key] += 1
        self.lists[key].append(amount)

    def get_as_dict(self):
        return self.values

    def get_as_lists(self):
        keys = [key for key in sorted(self.values)]
        return keys, [self.values[key] for key in keys]

    def get_all(self):
        return self.lists

    def get_stats(self, key):
        return np.mean(np.array(self.lists[key])), np.var(np.array(self.lists[key]))

    def sum(self):
        return np.sum([self.get_as_lists()[1]])

    def percentage(self, outof100=False, decimals=None):
        outp = {}
        for key in self.values:
            res = self.values[key]/self.sum()
            if outof100: res = res*100
            if not decimals is None: res = np.round(res,decimals)
            outp[key] = res
        return outp

class Substance:
    def __init__(self,namestr,array):
        self.fullname = namestr
        self.npeaks = np.shape(self.ma)[1]
        self.der = np.array([(array[1,i+1]-array[1,i])/(array[0,i+1]-array[0,i]) for i in range(np.shape(array)[1]-1)])
        i0 = func.GetIndex(np.int32(self.der<0),0)
        target = np.hstack((np.ones((i0)),np.zeros((np.size(self.der)-i0))))
        self.der = self.der*(np.int64(np.abs(self.der)<cf.maxleftslope)*target + (1-target))
        self.sder = func.Smooth(self.der,cf.itersmooth)
        self.date_created = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        
    def plot(self,showlim=True,showma=True,tof=None,peaklabs=True):
        if tof is None: tof = self.intof
        plt.figure()
        plotter.Plotter(self,self.fullname,showlim=showlim,showma=showma,tof=tof,peaklabs=peaklabs,axlabsin=False,vlines=None,ax=None)
        plt.show()

    def arr_t2E(self, arr):
        if np.shape(arr)[0] < 2:
            return basics.t2E(arr, self.mode)
        else:
            return np.vstack((basics.t2E(arr[0], self.mode), arr[1:]))

    def arr_E2t(self, arr):
        if np.shape(arr)[0] < 2:
            return basics.E2t(arr, self.mode)
        else:
            return np.vstack((basics.E2t(arr[0], self.mode), arr[1:]))

    def arr_dt2dE(self):
        pass

    pick = func.pick

class Data(Substance):
    def __init__(self,namestr,array,peaksdict=None):
        self.atom, self.symb, self.mass, self.mode = basics.InterpretName(namestr)
        self.intof = False
        self.xbounds = cf.xbounds()
        self.ybounds = cf.ybounds(self.symb, self.mode)
        self.spectrum = array
        self.spectrum_tof = self.arr_E2t(self.spectrum)
        self.xmagnitude = 'Energy (eV)'
        self.ymagnitude = 'Cross Section (b)'
        self.ma, self.mai = func.maxima(array, self.xbounds, self.ybounds, 0)
        self.ma_tof = self.arr_E2t(self.ma)
        super().__init__(namestr,array)
        self.peaks = peaksdict or func.propsisot(self, cf.pack())
        self._seterrors()
        
    def _seterrors(self):
        self.errors = [self.peaks[i] for i in self.peaks if self.peaks[i].xlims == (0.,0.)]

    def infopeaks(self):
        from spectra_FileHandlers import infoone
        print()
        print(self.fullname, 'PEAKS:','='*56)
        print('{:>4s}  {:>10s}        {:<10s} {:>10s} {:>17s} {:>17s}\n'.\
                    format('Rk.','Energy (eV)', 'TOF (us)', 'Integral','Peak width','Peak height'))
        for line in infoone(self): print(line)
        print('='*78)
        print()
    
    def plotsingle(self,num):
        plt.figure()
        plotter.plotone(self,num,title='{} #{}'.format(self.fullname, num))
        plt.show()

    def get_from_peaks(self,attr):
        if peakattr.has(attr): return np.array([getattr(self.peaks[i],attr) for i in sorted(self.peaks)])

    def getclosest(self,inp,in_tof=True):
        malist = self.get_from_peaks('ma_tof') if not in_tof else self.get_from_peaks('ma')

    def edit(self):
        editing = func.EditPeaks(self)
        if editing != dict():
            self.peaks = func.propsisot(self, cf.pack(), setx=editing)
            self._seterrors()
            self.date_edited = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    def delete(self):
        deleting = func.DeletePeaks(self)
        if deleting != []:
            self.peaks = func.sorting({i: self.peaks[i] for i in self.peaks if not i in deleting})
            self._seterrors()
            self.date_edited = datetime.now().strftime("%d/%m/%Y %H:%M:%S")


    plotpeaks = plotter.plotpeaks



class Mix(Data):
    def __init__(self,namestr,array,abund,peaksdict=None):
        super().__init__(namestr,array,peaksdict)
        self.abundances = abund
        self.components = list(abund.keys())

    def recompute(self,**kwargs):
        self.__init__(self.fullname, self.spectrum, self.abundances, None)
        #self.peaks = propsisot(self, cf.pack(**kwargs))
        #self.err = [self.peaks[i] for i in self.peaks if self.peaks[i].xlims == (0.,0.)]

class Element(Mix):
    kind = 'element'
    def __init__(self,namestr,array,abund,peaksdict=None):
        super().__init__(namestr,array,abund,peaksdict)

class Compound(Mix):
    kind = 'compound'
    def __init__(self,namestr,array,abund,peaksdict=None):
        super().__init__(namestr,array,abund,peaksdict)
    
class Isotope(Data):
    kind = 'isotope'
    def __init__(self,namestr,array,peaksdict=None):
        super().__init__(namestr,array,peaksdict)

    def recompute(self,**kwargs):
        self.__init__(self.fullname, self.spectrum, None)
        #self.peaks = propsisot(self, cf.pack(**kwargs))
        #self.err = [self.peaks[i] for i in self.peaks if self.peaks[i].xlims == (0.,0.)]

class Sample(Substance):
    kind = 'sample'
    def __init__(self,namestr,arrayin,mode=None):
        self.intof = True
        self.xbounds = None
        self.ybounds = None
        if mode is None:
            inp = input('Mode for {}: (1: n-g; 2: n-tot) >'.format(namestr))
            self.mode = {'1':'n-g', '2':'n-tot'}.get(inp,cf.default_smode)
            del inp
        else:
            self.mode = mode
        self.spectrum_tof = arrayin
        self.spectrum = self.arr_t2E(self.spectrum_tof)
        self.xmagnitude = 'ToF (us)'
        self.ymagnitude = 'Counts'
        self.ma_tof, self.mai_tof = func.maxima(arrayin, None, None, cf.itersmoothsamp)
        self.ma = self.arr_t2E(self.ma_tof)
        super().__init__(namestr,arrayin)
        func.sampprocess(self)

class Peak:
    kind = 'peak'
    def __init__(self,info):
        assert len(info) == peakattr.size, 'info tuple must have the same lenght as peakattr, but' + len(info) +' / '+ peakattr.size
        for i in range(peakattr.size):
            setattr(self,peakattr.get(i),info[i])

    pick = func.pick