#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:08:15 2019

@author: Ivan Alsina Ferrer
"""

#This is a program written by Ivan Alsina Ferrer (ivanalsinaferrer@gmail.com)
#during a summer internship at STFC (Oxfordshire, UK) on the period Jun-Aug 2019.
#It intends to create a spectrum database
#Please read the documentation in the spectra_Manual.docx file.

import os
import sys
from tqdm import tqdm
from datetime import datetime
import time
import shutil
import pickle

import numpy as np
import scipy.integrate as spint#
import matplotlib.pyplot as plt

from .spectra_Basics import isd, isf, isfx
from .spectra_Objects import Catalog
from .spectra_InitSettings import cf, err, paths
from .spectra_FileHandlers import pload

#This allows plots to be detached from the command line.
plt.ion()

#Required directories:


assert isd('data'), 'Missing directory: {}'.format(paths.data)

if not isd('output'):
    os.mkdir(paths.output)

if not isd('input'):
    os.mkdir(paths.input)

if not isd('load'):
    os.mkdir(paths.load)

#If 'spcat.pickle' pickle exists, ask.
#Either load it or create the insntace
if paths.isfx('spcat.pickle', 'load'):
    inp = input('Load catalog from file? ([y]/n) >')
    if not inp in ['n','no','q','quit']:
        spcat = pload()
    elif not inp in ['q','quit']:
        spcat = Catalog()
else:
    spcat = Catalog()

