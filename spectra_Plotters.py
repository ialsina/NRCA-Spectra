#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:08:15 2019

@author: ivan
"""

import os
import sys
from tqdm import tqdm
from datetime import datetime
import time
import shutil

import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as plt

from spectra_Basics import *
#from spectra_Settings import *
from spectra_Objects import *
from spectra_Finders import *
from spectra_InitSettings import cf


def Plotter(data,tit='',showlim=0,showma=0,tof=False,axlabsin=False,peaklabs=False,vlines=None,ax=None,**kwargs):
    """The main plotting function. It does not create a figure, just its content.
    inputs:
        - data: instance of Substance class or either dictionary of instances.
            data series whose spectrum either in energy or tof are to be plotted.
            If it is a dictionary, the keys are displayed in the legend.
        - tit: str
            Figure title
        - showlim: bool or int
            True or 1 / False or 0: Are bounds used for peak detection shown?
        - showma:  bool or int
            True or 1 / False or 0: If True, maxima are marked with red points.
        - tof: bool or int
            True or 1: a graph in Time of Flight is shown.
            False or 0: a graph in Energy is shown.
        - axlabsin: bool or int
            if True, the axes labels are displayed inside the plotting area.
        - peaklabs: bool or int
            if True, the peak labels are displayed.
        - vlines: tuple or None
            tuple of x coordinates where a vertical line is to be displayed
            if None, disabled.
        - ax: axis instance of matplotlib or None
            if provided, using that axis. If None, plot is done in plt.gca()
        - additional keyword arguments can be passed.
    ouptuts:
        - ax: the plotted axis.
    Notice: usually there is no need to store the output as a variable or a new instance.
    The usual procedure is to create a figure outside, call this function and show it, i.e.,
    plt.figure()
    Plotter(...)
    plt.show()"""
    islegend = False
    if type(data) not in [list,dict]:
        data = [data]
    elif type(data) == dict:
        islegend = True
        legend = list(data.keys())
        data = list(data.values())
    if tof:
        xscale = 'linear'
    else:
        xscale = 'log' if not data[0].kind == 'peak' else 'linear'
    yscale = 'linear' if data[0].kind in ['peak', 'sample'] else 'log'

    def units(arr,intof,mode):
        if tof^intof:
            if tof: return E2t(arr,mode)
            if not tof: return t2E(arr,mode)
        else:
            return arr

    ax = ax or plt.gca()
    ax.set_title(tit)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel('Time of flight (us)' if tof else 'Energy (eV)')
    ax.set_ylabel(data[0].ymagnitude)
    if axlabsin:
        ax.xaxis.set_label_coords(.5,.05)
        ax.yaxis.set_label_coords(.05,.5)
    for i in range(len(data)):
        series = data[i]
        label = legend[i] if islegend else None
        arr = series.pick('spectrum',tof)
        ax.plot(arr[0], arr[1], label=label, lw=0.8, color=kwargs.get('color'))
        
        if showlim:
            ax.axvline(x=units(series.xbounds[0],series.intof,series.mode), color='black', lw=0.5)
            ax.axvline(x=units(series.xbounds[1],series.intof,series.mode), color='black', lw=0.5)
            ax.axhline(y=series.ybounds,color='black',lw=0.5)
        
        if not vlines is None:
            for vline in vlines:
                ax.axvline(x=units(vline,series.mode),color='black',lw=0.5,ls='--')
        
        if showma:
            ma = series.pick('ma', tof)
            ax.scatter(ma[0], ma[1], color='red', marker='o', s=1)
            if len(data) == 1:
                plt.text(.95,1.02,'# Peaks: {}'.format(np.shape(series.ma)[1]),fontsize=8,ha='right',va='center',transform=ax.transAxes)
        if peaklabs:
            if series.kind != 'sample':
                for peak in list(series.peaks.values()):
                    coords = peak.pick('coords',tof)
                    ax.annotate(str(peak.num),(coords[0], coords[1]),fontsize=7)
            else:
                for peak in range(np.shape(series.ma_tof)[1]):
                    coords = series.pick('ma', tof)[0:2,peak]
                    ax.annotate(str(peak), (coords[0], coords[1]), fontsize=7)
    if islegend: ax.legend()

    return ax


def Plot(self,same=True):
    """This is an instance of Catalog.
    Asks what elements are to be shown together and plots them together.
    inputs:
        self: Category instance.
        same: boolean.
            True: what is there to be plotted will be showed together, in the same figure.
            False: an independent figure will be opened for each item.
    """
    Dict = self.Substances()
    querylist = Select(self.get_as_dict(),recursive=True)
    if querylist == []: return None
    querydict = dict()
    tof = not AskAxis()
    maxs = int(AskLim())
        
    for queryel in querylist:
        if queryel in Dict: querydict[queryel] = Dict[queryel]
    
    if not same:
        for queryel in querydict:
            plt.figure()
            Plotter(querydict[queryel],tit=queryel,tof=tof,showlim=maxs,showma=maxs,peaklabs=1)
            plt.show()
    
    elif same:
        plt.figure()
        Plotter(querydict,tit='Query results',tof=tof,showlim=maxs,showma=maxs,peaklabs=1)
        plt.show()


def plotone(self,num,side=None,title=None,prplot=cf.prplot,ax=None):
    """Function called by a method member of Data.
    Plots a single peak:
    inputs:
        - self: Data instance
        - num: int
            peak number to plot
        -side: str or None
            text to be placed on the left side of the plot
        - title: str or None
            text to be placed on top of the plot
        - prplot: int
            number of mesh points each side of the peak to plot
        - ax: plt axes instance or None
            if given, plots in that axes. If None, plots in plt.gca()"""
    ax = ax or plt.gca()
    if not self.peaks[num].xlims is None:
        indx = GetIndex(self.spectrum[0],self.peaks[num].xlims)
        x0 = indx[0] - prplot
        x1 = indx[1] + prplot +1
    else:
        x0 = self.mai[num] - 120
        x1 = self.mai[num] + 120
    if not title is None: plt.title(title)
    if not side is None: plt.text(.02,.5,side,fontsize=10,color='red',ha='left',va='center',transform=ax.transAxes)
    ax.plot(self.spectrum[0,x0:x1],self.spectrum[1,x0:x1])
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    if not self.peaks[num].xlims is None:
        for edge in self.peaks[num].xlims:
            ax.axvline(x=edge,color='black',lw=0.5)
    else:
        ax.axvline(x=ma[0,num],color='red',lw=0.5)
    plt.text(.95,.5,'{}-{}'.format(self.peaks[num].peakreason[0],self.peaks[num].peakreason[1]),fontsize=8,ha='right',va='center',transform=ax.transAxes)
    #ax.set_xlabel('Energy (eV)')
    #ax.set_ylabel('Cross section (b)')
    return ax

def plotpeaks(self,save=False,showlim=1,showma=1,peaklabs=True,**kwargs):
    """Function called by a method member of Data.
    Plots all the peaks in the Data instance in a fancy 4-column way.
    inputs:
        - self: Data instance
        - save: bool
            True: saves the dataplot
            False: shows the dataplot
        - the rest of the inputs are the ones documented in Plotter"""
    nrows = int((len(self.mai)-9)/4)+4 if len(self.mai) > 8 else (2 if len(self.mai) in [0,1,2,4] else 3)
    ncols = 4 if len(self.mai) > 5 or len(self.mai) == 4 else (2 if len(self.mai) == 0 else 3)
    
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=(4*ncols,3*nrows))
    gs = axes[1, 2].get_gridspec()
    for ax in axes[0:2,0:2].flatten():
        ax.remove()
    axes = axes.flatten()
    if ncols==3: axes = np.delete(axes,[0,1,3,4])
    if ncols==4: axes = np.delete(axes,[0,1,4,5])
    if len(self.mai) > 8:
        for ax in axes[len(self.mai):ncols*nrows]: ax.remove()
    for i in range(len(self.mai)):
        plotone(self,i,side=str(i),ax=axes[i])
    axbig = fig.add_subplot(gs[0:2,0:2],xscale='log',yscale='log')
    Plotter(self,self.fullname,showlim=1,showma=1,peaklabs=peaklabs,axlabsin=True,ax=axbig)
    fig.tight_layout()

    if not save: plt.show()
    if save: plt.savefig(os.path.join(os.getcwd(),'output','Peaks_'+self.fullname+'.png'))

def PlotBars(samp,isots=dict()):
    """Function called by a method member of Catalog.
    Given a sample and a dictionary of Data instances, plots a bar plot where the heights are proportional to the Data peak intensities.
    inputs:
        - samp: sample instance
        - isots: Data instances dictionary"""
    fig, ax = plt.subplots(1,1)
    lowest = np.min(samp.background_tof[1])
    if len(isots) != 0:
        supreme_s = np.max(samp.spectrum_tof[1])
        supreme_b = np.max(np.concatenate([isots[isot].get_from_peaks('integral') for isot in isots]))
        scaling = supreme_s/supreme_b*cf.bar_scaling
        for isotname in isots:
            isot = isots[isotname]
            ax.bar(isot.get_from_peaks('center_tof'), isot.get_from_peaks('integral')*scaling, width=cf.bar_width, bottom=lowest, label=isot.fullname)
        ax.axhline(y=lowest, lw=0.5, color='black')
    Plotter(samp,tit=samp.fullname,tof=True,showlim=False,showma=True,peaklabs=1, ax=ax, color='black')
    ax.plot(samp.background_tof[0], samp.background_tof[1], color='yellow')
    if not isots==dict(): ax.legend()
    fig.show()