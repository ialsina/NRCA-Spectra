#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:55:06 2020

@author: ivan
"""

from os import path as p
from os import getcwd

cwd = p.dirname(p.realpath(__file__))

paths_ = {
          'cwd'       :   cwd,
          'source'    :   p.join(cwd, 'src'),
          'data'      :   p.join(cwd, 'data'),
          'samples_n-tot':        p.join(cwd, 'samples_n-tot'),
          'samples_n-g'  :        p.join(cwd, 'samples_n-g'),
          'load'      :   p.join(cwd, 'load'),
          'input'     :   p.join(cwd, 'input'),
          'output'    :   p.join(cwd, 'output'),
          'history'   :   p.join(cwd, 'history'),
         }
