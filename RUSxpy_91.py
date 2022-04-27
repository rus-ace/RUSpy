#!/usr/bin/env python3.8
# Control program for the ACE Resonant Ultrasound Spectrometer
# Ref: Reviews of Scientific Instruments Vol.90 Issue 12 Dec 2019 pgs. 121401 ff.
# Copyright (C) 2022  Albert Migliori
# Based heavily on the Control program for the Red Pitaya vector network analyzer
# Copyright (C) 2021  Pavel Demin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. 

from asyncio import sleep
from asyncio.windows_events import NULL
from ctypes.wintypes import SC_HANDLE
from multiprocessing.sharedctypes import Value
from operator import pos
from datetime import datetime
#from os import wait
from pickle import NONE, TRUE
from shutil import register_unpack_format
from sqlite3 import Row
from statistics import mode
import sys
import struct
from tracemalloc import stop
from turtle import onrelease
from unittest.mock import NonCallableMock
import warnings
from xmlrpc.server import XMLRPCDocGenerator
from matplotlib import pyplot as plt
import os
from pathlib import Path
from pyparsing import lineEnd

import serial
import serial.tools.list_ports
import time

from functools import partial
from matplotlib.backend_bases import Event

import numpy as np
import math
import csv

import matplotlib
from numpy.core.arrayprint import str_format
from numpy.lib.type_check import imag, real
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import Formatter, FuncFormatter
from matplotlib.widgets import Cursor, MultiCursor
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
#import matplotlib.pyplot as plt

from PyQt5.uic import loadUiType
from PyQt5.QtCore import QRegExp, QTimer, QSettings, QDir, Qt
#from PyQt5.QtCore import Null, QRegExp, QTimer, QSettings, QDir, Qt
from PyQt5.QtGui import QRegExpValidator, QPalette, QColor, QBitmap, QPixmap
from PyQt5.QtWidgets import QApplication, QDoubleSpinBox, QMainWindow, QMessageBox, QDialog, QFileDialog, QPushButton, QLabel, QSpinBox
from PyQt5.QtNetwork import QAbstractSocket, QTcpSocket
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSplashScreen


Ui_rus, QMainWindow = loadUiType('rusx.ui')

class Measurement:
  def __init__(self, start, stop, size):
    self.freq = np.linspace(start, stop, size)
    self.data = np.zeros(size, np.complex64)
    #self.period = 62500  
    
class NavigationToolbar(NavigationToolbar):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('','Save')]   

class FigureTab:
  def __init__(self, layout, rus):
    event = 'motion_notify_event'
    # event = 'mouse_event'
    # create figure
    global comFlag
    self.figure = Figure()
    if sys.platform != 'win32':
      self.figure.set_facecolor('none')
    self.canvas = FigureCanvas(self.figure)
    self.canvas.mpl_connect('button_press_event', self.onpress)
    self.canvas.mpl_connect('button_release_event', self.onrelease)
    layout.addWidget(self.canvas)
    self.toolbar = NavigationToolbar(self.canvas, None, False)
    self.toolbar.layout().setSpacing(40)
    #remove subplots action
    actions = self.toolbar.actions()
    layout.addWidget(self.toolbar)
    self.tf = ''
    ##################################################################################
    #add set phase
    self.mode = 'reim'
    
    if comFlag == True:
      ser.write(b's\n')     # set point
      time.sleep(0.1)
      td = ser.readline()
      td = td.decode('UTF-8')
      self.tidF = td.strip('\n')
    else: self.tidF = '0.0'

    self.phase = QDoubleSpinBox()
    self.phase.setFixedHeight(40)
    self.phase.setSingleStep(5)
    self.phase.setMaximum(180)
    self.phase.setMinimum(-180)
    self.phase.setDecimals(1)
    self.phase.valueChanged.connect(self.phase_set)
    self.toolbar.addWidget(self.phase)
    #add value at crsor
    self.temp = QDoubleSpinBox()
    self.temp.setFixedHeight(50)
    self.temp.setReadOnly(True)
    self.temp.setButtonSymbols(2)
    self.temp.setDecimals(4)
    self.toolbar.addWidget(self.temp)
    #add crsor
    self.crsor = QDoubleSpinBox()
    self.crsor.setDecimals(3)
    self.crsor.setFixedHeight(50)
    self.crsor.setFixedWidth(125)
    self.crsor.setStyleSheet('background-color:pink')
    self.crsor.valueChanged.connect(self.plot_crsr)
    self.toolbar.addWidget(self.crsor)
    #add mark peak button
    self.peak = QPushButton('Mark')
    self.peak.setStyleSheet('background-color:orange')
    self.peak.setFixedWidth(80)
    self.peak.clicked.connect(self.mark)
    self.peak.clicked.connect(self.wFreq)
    self.toolbar.addWidget(self.peak)
    layout.addWidget(self.toolbar)
    #add unmark button
    self.unpeak = QPushButton('unMark')
    self.unpeak.setStyleSheet('background-color:red; color:rgb(255, 255, 255)')
    self.unpeak.setFixedWidth(80)
    self.unpeak.clicked.connect(self.unmark)
    self.unpeak.clicked.connect(self.wFreq)
    self.unpeak.clicked.connect(self.reset)
    self.toolbar.addWidget(self.unpeak)
    layout.addWidget(self.toolbar)
    #add rescale button
    self.plotButton = QPushButton('Replot')
    self.plotButton.setStyleSheet('background-color:black; color:rgb(255, 255, 255)')
    self.plotButton.setFixedWidth(80)
    self.plotButton.clicked.connect(self.TcontStop)
    self.plotButton.clicked.connect(self.reset)
    self.plotButton.clicked.connect(self.wFreq) #writes freq data
    self.toolbar.addWidget(self.plotButton)
    axes1 = self.figure.add_subplot(1,1,1)
    axes2 = self.figure.add_subplot(1,1,1)
    self.axes1 = axes1
    self.axes2 = axes2
    self.zt = 0.0
    self.rus = rus
    global ix0,ix1
    ix0 = 0
    ix1 = len(self.rus.reim.freq)
    self.ix0 = ix0
    self.ix1 = ix1
    self.izl = 0
    self.izr = self.ix1
    global rusf
    rusf = []
    self.TsweepTimer = QTimer()
    self.TsweepTimer.timeout.connect(self.Tcont)
    ###############################################################################
    self.rus.strT.clicked.connect(self.Tcont) #plots temperature
    self.rus.stpT.clicked.connect(self.TcontStop)
    self.rus.CM.clicked.connect(self.clearrusf)
    self.rus.CM.clicked.connect(self.wFreq)####################################
    self.rus.AP.clicked.connect(self.autoPeak)
    self.rus.AP.clicked.connect(self.reset)
    self.rus.ClearT.clicked.connect(self.tsReset)
  
  def tsReset(self): #resets time axis for temperature plot
    global ts, tout,ttout
    ts = time.time()
    self.axes1.cla()
    self.axes2.cla()
    self.axes1.xaxis.grid()
    self.axes1.yaxis.grid()
    tout = []
    ttout = []
    tout.append([])
    ttout.append([])

  def Tcont(self): #Starts the continuous plot of temperature
        global tout, ttout
        if self.rus.strT.isEnabled():
         self.rus.strT.setEnabled(False)
         self.rus.stpT.setEnabled(True)
         if (tout[0] != []):
            ymin = np.min(tout)
            ymax = np.max(tout)
            ylim = [ymin-.01,ymax+.01]
            self.axes1.tick_params(axis= 'x', colors = 'black', labelsize=16)
            self.axes1.tick_params(axis= 'y', colors = 'black', labelsize=16)
            self.axes1.set_xlabel('Time(s)', fontsize=20)
            self.axes1.set_ylabel('Temperature('+sx+')', fontsize=20)     
            self.axes1.set_xlim([0,ts1])
            self.axes1.set_ylim(ylim)
            self.axes1.plot(ttout,tout,color = 'black',  linestyle = 'solid', linewidth = 2, marker = 'o' , markersize =2)
            self.TsweepTimer.start(5000)
        else:
            self.plot_tout() #plots a temperature point each dTms

  def TcontStop(self): #stops continuous temperature plots
    self.TsweepTimer.stop()
    self.rus.strT.setEnabled(True)
    self.rus.stpT.setEnabled(False)
    #print(tout)
    #self.update(self.mode)

  def autoPeak(self):# Finds most peaks
    global aPf, rusf, sx, comFlag, dt, id
    mode = self.mode
    if aPf == False:
      aPf = True
      return
    aPf = True
    pdata = self.rus.reim.data
    pfreq = self.rus.reim.freq
    pdata = abs(pdata)
    r = len(pfreq)
    pmax = max(pdata) 
#find rms noise
    l = 16
    pfm = pdata.argmin()
    pfmin = pfreq[pfm]
    qbar = 0.0
    for q in range(1, l+1):
      qbar = qbar +pdata[q]
    qbar = qbar/l
    qrms = 0.0
    for q in (0, l):
      qrms = qrms + (pdata[q]-qbar)**2
    qrms = math.sqrt(qrms/l)
#test peak
    n = 6
    percent = self.rus.percent.value()
    for i in range(n, r-n-4):
      x1 = pfreq[i]; y1 = pdata[i]
      x2 = pfreq[i+1]; y2 = pdata[i+1]
      x3 = pfreq[i+2]; y3 = pdata[i+2]
      y1p = pdata[i-1]
      y3p = pdata[i+3]
      m = self.rus.noise.value()
      m = m*qrms
      if (y2 > 0.01*percent*pmax):
        if (y2 > y1):
            if (y2 >y3):
                if (y3 > y3p+m):
                    if (y1 > y1p+m):
                        #print(y1, y1p, y2, y3p, y3, 0.01*percent*pmax, x2, m)
                        p = []
                        p.append([])
                        p[-1].append(x2)
                        p[-1].append(y2)
                        p[-1].append(1.0)
                        dup = False
                        for k in range(0,len(rusf)):
                            if x2 == rusf[k][0]:
                                dup = True
                        if dup == False: rusf = rusf + p
                        aPf = False
                        self.wFreq()
                        self.update(mode)
  
  def plot_curves(self, izl, izr, xhgt, tf): #main plotting routine
    matplotlib.rcdefaults()
    matplotlib.rcParams['axes.formatter.use_mathtext'] = True
    self.figure.clf()
    self.figure.text(0.04, 0.005,r'       $\phi$              Z             f', fontsize = 18)
    freqP = self.rus.reim.freq
    sstep = (freqP[2]-freqP[1])
    dataP = self.rus.reim.data
    swF = self.rus.swF
    if swF == True: #ensures runs only once
        izl = 0
        #izr = self.rus.sweep_size
        izr = len(freqP)
        self.ix0 = izl
        self.ix1 = izr
        self.izl = 0
        self.izr = izr
        self.rus.swF = False
    if izl > 0:
         freqP = self.rus.reim.freq[izl:izr]
         dataP = self.rus.reim.data[izl:izr]
    data1 = np.real(dataP)
    data2 = np.imag(dataP)
    data3 = np.absolute(dataP)
    self.data3 = data3
    y = 1.07*np.max(data3)
    if xhgt[0] == 0: xhgt = [-y,y]
    h = [-y,y]
    lim = y
    if izl > 0: h = xhgt
    lim = xhgt
    ph = 0.0
    phD = 0.0
    co = 1.0
    si = 0.0
    #operate on phase
    if self.mode == 'reim':
        if tf == 'set':
            phD = self.phase.value()
            ph = math.radians(phD)
            co = math.cos(ph)
            si = math.sin(ph) 
    d1 = co*data1 - si*data2
    d2 = si*data1 + co*data2
    #find maximum
    i = np.argmax(data3)
    freqmax = float(freqP[i])
    vmax = data3[i]
    freq3 = [freqmax,freqmax]
    reim_Max = [-vmax,vmax]
    #add axes
    axes1 = self.figure.add_subplot(1,1,1)
    axes1.cla()
    axes1.xaxis.grid()
    axes1.yaxis.grid()
    axes1.set_xlabel('kHz',fontsize=20)
    axes1.set_ylabel('real',fontsize=20)
    #Setup crsor spinbox
    self.crsor.setMaximum(freqP[-1])
    self.crsor.setDecimals(3)
    self.crsor.setMinimum(freqP[0])
    self.crsor.setSingleStep(sstep)##########################3
    self.crsor.setValue(freqmax)
    #setup axes
    hgt = data3[i]
    xlim = [freqP[0],freqP[-1]]
    #axes1.set_xlim(xlim)
    axes1.tick_params('y', color = 'black', labelcolor = 'black', labelsize=16)
    axes1.tick_params('x', color = 'black', labelcolor = 'black', labelsize=16)
    axes1.yaxis.label.set_color('blue')
    axes1.set_ylabel('imag',fontsize=20)
    axes1.set_xlim(xlim)
    if self.mode == 'reim': 
      lim = h        
    if self.mode == 'absm': 
      lim = [z, h[1]]
    #axes3 = axes1.twinx()
    self.axes1 = axes1
    #lim = [-0.1,0.1]#############################################
    if lim is not None: 
        axes1.set_ylim(lim)
    # initiate plots
    if self.mode == 'reim': 
        self.temp.setValue(np.abs(data3[i]))
        axes1.set_ylabel('Volts',fontsize=20)
        axes1.yaxis.label.set_color('black')
        axes1.plot(freqP, d1, color = 'blue')
        axes1.plot(freqP, d2, color = 'red')  ##
    if self.mode == 'absm': 
        self.temp.setValue(np.abs(data3[i]))
        axes1.set_ylabel('Volts',fontsize=20)
        axes1.yaxis.label.set_color('black')
        axes1.plot(freqP, data3, color = 'red', label = 'magnitude')
    #plot line at max
    axes2 = axes1.twinx()
    self.axes2 = axes2
    axes2.tick_params(right = False)
    axes2.yaxis.set_ticklabels([])
    axes2.tick_params(right = False)
    #enable cursor       
    self.cur = Cursor(axes2, horizOn=True, vertOn=True, useblit=True, color = 'blue', linewidth = 1)
    if izl > 0:self.plot_crsr()
    self.wFreq()
    self.canvas.draw()
  
  def mark(self): #marks peaks in window
    global rusf
    x = []
    x.append([])
    #z = round(self.crsor.value(),3)
    z = self.crsor.value()
    for y in rusf:
        if z == y[0]: return 
    x[-1].append(z)
    x[-1].append(self.temp.value())
    x[-1].append(1.0)
    rusf = rusf + x

  def unmark(self): #unmarks peaks in window
    global rusf
    if rusf == []:return
    span = self.axes1.get_xlim()
    fzl = span[0]
    fzr = span[1]
    rows = len(rusf)
    rusa = rusf
    k = 0
    for i in range(0,rows):
      i = i - k
      if i > (len(rusa) - 1):break
      m = rusa[i][0]     
      if m >= fzl and m <= fzr:
        rusa = np.delete(rusa, i, 0)
        k = k + 1
    if k == 0:return
    rusfb = rusa.tolist()
    rusf = rusfb

  def clearrusf(self): #clears all peaks (stored in rusf)
    global rusf
    rusf = []
    self.wFreq()
    self.update(self.mode)

  def wFreq(self):#writes found peaks to file
    global rusf, sx, comFlag, dt, id, fn
    if self.rus.auto == True:return
    if rusf == []:
      return
    dat = dt
    id0 = self.rus.Tkc.value()
    tidF = str(int(id0+0.5))
    #tidF = str(self.tidF)
    #tid = str(self.rus.setPt.value())
    fn = os.path.join(path, tidF + '_' +sx + '_' + dat + '_F.dat')
    b = sorted(rusf)
    bx = [i[0] for i in b]
    by = [i[1] for i in b]
    b1 = [i[2] for i in b]
    if bx != []:
       with open(fn, 'w', newline='') as f:
          f.write('%8.2f' % id0)
          f.write(' ' +sx +'\n')
          for j in range(len(b)):
              f.write('%3.7f  %5.4f  %5.1f \n' % (0.001*bx[j], by[j], b1[j]))
       f.close()  
    self.axes1.plot(bx,by,color = 'black', marker = matplotlib.markers.CARETDOWNBASE , markersize =10, linestyle = 'None')
    self.canvas.draw()

  def update(self, mode):
    if self.rus.EnableCursor == False:
      self.rus.EnableCursor = True
      return
    getattr(self, 'update_%s' % mode)()

  def reset(self):#resets plot
    izl = 0
    izr = len(self.rus.reim.data)
    self.ix0 = 0
    self.ix1 = len(self.rus.reim.data)
    self.izl = 0
    self.izr = len(self.rus.reim.data)
    tf = ''
    h = np.max(abs(self.rus.reim.data[izl:izr]))
    xhgt = [-1.07*h,1.07*h]
    self.rus.swF = True
    self.plot_curves(izl, izr, xhgt, tf)

  def onpress(self,event): #starts window selction
    #gets cursor coordinates
     if event.inaxes:
      fx = event.xdata
      ix = self.find_dPoint(fx)[0]
      self.ix0 = ix

  def onrelease(self,event):#ends window selction#gets cursor coordinates
    if event.inaxes:
        fx = event.xdata
        ix = self.find_dPoint(fx)[0]
        self.ix1 = ix
        izl = self.ix0
        izr = self.ix1
        if izr <= izl: return
        #if izl > izr: return
        fl = self.rus.reim.freq[izl]
        fr = self.rus.reim.freq[izr]
        lim = [fl,fr]
        self.axes1.set_xlim(lim)
        xh = np.max(abs(self.rus.reim.data[izl:izr]))
        xhgt = [-1.07*xh, 1.07*xh]
        tf = self.tf
        self.izl = izl
        self.izr = izr
        self.il = izl
        self.ir = izr
        self.plot_curves(izl, izr, xhgt, tf)
        self.wFreq()

  def plot_crsr(self):#plots vertical line at peak position and can be moved
    if self.rus.EnableCursor == True:
        vc = self.crsor.value()
        axes2 = self.axes2
        hgt = self.axes1.get_ylim() #remember hgt before cla
        wdt = self.axes1.get_xlim()
        axes2.cla()
        axes2.tick_params(right = False)
        axes2.yaxis.set_ticklabels([])
        fre = self.find_dPoint(vc)[1]
        j = self.find_dPoint(vc)[0]
        posf = [fre,fre]
        h = np.abs(self.rus.reim.data[j])
        axes2.set_ylim(hgt) #reset y limit
        axes2.set_xlim(wdt)
        if self.mode == 'absm':val = [z, h]
        if self.mode == 'reim':val = [-h, h]
        self.temp.setValue(h)
        axes2.plot(posf, val, color = 'black', linewidth = 1, linestyle = 'solid')
        self.canvas.draw()
        self.wFreq()
        return
  
  def phase_set(self): #sets phase
    if self.mode == 'absm': return
    tf = 'set'
    self.tf = tf
    izl = self.ix0
    izr = self.ix1
    #print(izl,izr)
    h = np.max(abs(self.rus.reim.data[izl:izr]))
    xhgt = [-1.07*h,1.07*h]
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()

  def plot_reim(self):# plots real and imginary components
    global ix0, ix1
    self.ix0 = ix0
    self.ix1 = ix1
    izl = self.ix0
    izr = self.ix1
    #print(izl,izr)
    h = np.max(abs(self.rus.reim.data[izl:izr]))
    xhgt = [-1.07*h,1.07*h]
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()

  def update_reim(self):
      self.plot_reim()
      #getattr(self, 'plot_%s' % self.mode)()
  
  def plot_tout(self): #continuous plot of temperaure
    global tout, ttout, comFlag, ts, sx
    self.TsweepTimer.start(5000)
    self.axes2.yaxis.set_ticklabels([])
    self.axes2.tick_params(right = False)
    if comFlag == False: 
        return
    ser.write(b't\n')     # write a string to get temperature
    td = ser.readline()
    td = td.decode('UTF-8')
    id = td.strip('\r\n')
    if id != '':
        id = float(id)
        tout[-1].append(id)
        ts1 = time.time()-ts
        ttout[-1].append(ts1)
        ymin = np.min(tout)
        ymax = np.max(tout)
        ylim = [ymin-.01,ymax+.01]
        self.axes1.tick_params(axis= 'x', colors = 'black', labelsize=16)
        self.axes1.tick_params(axis= 'y', colors = 'black', labelsize=16)
        self.axes1.set_xlabel('Time(s)', fontsize=20)
        self.axes1.set_ylabel('Temperature('+sx+')', fontsize=20)     
        self.axes1.set_xlim([0,ts1])
        self.axes1.set_ylim(ylim)
        self.axes1.plot(ts1,id,color = 'black',  linestyle = 'solid', linewidth = 2, marker = 'o' , markersize =2)
        self.rus.Tkc.setValue(id)
        if sx == 'K': 
            ser.write(b'c\n')     # control block temperature for Thermoelectric stage
            bd = ser.readline()
            if bd == b'' :
                pass
            bd = bd.decode('UTF-8')
            bv = bd.strip('\n')
            bv = float(bv)
            self.rus.bT.setSuffix(sx)
            self.rus.bT.setValue(bv)
        self.canvas.draw()

  def plot_absm(self): #plots absolute value
    global ix0, ix1
    self.ix0 = ix0
    self.ix1 = ix1
    izl = self.ix0
    izr = self.ix1
    h = np.max(abs(self.rus.reim.data[izl:izr]))
    xhgt = [-1.07*h,1.07*h]
    self.mode = 'absm'
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()
    ##############################################################################################################

  def update_absm(self):
     self.plot_absm()

  def find_dPoint(self, vc):
      #finds the closest data point to the crsr
    fv = vc
    f = self.rus.reim.freq
    i0 = 0
    dl = abs(fv-f[0])
    for i in range(len(f)-1):
        i = i+1
        if abs(fv-f[i]) < dl: 
            dl = abs(fv-f[i])
            i0 = i
    fre = f[i0]
    j = i0
    return j, fre

class rus(QMainWindow, Ui_rus):
  graphs = ['reim', 'absm']
  def __init__(self):
    global comflag, sx
    super(rus, self).__init__()
    self.setupUi(self)
    global comFlag, comport
    settings = QSettings('rus.ini', QSettings.IniFormat)
    self.read_cfg_settings(settings)
    self.idle = True
    self.reading = False
    self.auto = False
    self.sweep_Hz = self.stepValue.value()
    self.sweep_start = self.startValue.value()
    self.sweep_stop = self.stopValue.value()
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/self.sweep_Hz)
    if(self.sweep_size)>32766:
      self.sweep_size = 32766
      self.sweep_Hz = 1000*(self.sweep_stop-self.sweep_start)/32766
    self.nSize.setValue(self.sweep_size)
    # buffer and offset for the incoming samples
    self.buffer = bytearray(16 * 32768)
    self.offset = 0
    self.data = np.frombuffer(self.buffer, np.complex64)
    # create measurements
    self.reim = Measurement(self.sweep_start, self.sweep_stop, self.sweep_size)
    self.mode = 'reim'  
    # create figures
    self.tabs = {}
    for i in range(len(self.graphs)):
      layout = getattr(self, '%sLayout' % self.graphs[i])
      self.tabs[i] = FigureTab(layout, self)
    # configure widgets
    #self.rateValue.addItems([ '1500','500', '150', '50', '15', '5', '2'])
    self.rateValue.addItems([ '1500','500', '150', '50'])
    #self.rateValue.addItems(['5000', '1000', '500', '100', '50', '10', '5', '1'])
    #rate = [10, 50, 100, 500, 1000, 5000, 10000, 50000][value]
    self.rateValue.lineEdit().setReadOnly(True)
    self.rateValue.lineEdit().setAlignment(Qt.AlignRight)
    for i in range(self.rateValue.count()):
      self.rateValue.setItemData(i, Qt.AlignRight, Qt.TextAlignmentRole)
    self.set_enabled(False)
    self.stopSweep.setEnabled(False)
    # create TCP socket
    self.socket = QTcpSocket(self)
    self.socket.connected.connect(self.connected)
    self.socket.readyRead.connect(self.read_data)
    self.socket.error.connect(self.display_error)
    # connect signals from widgets
    self.connectButton.clicked.connect(self.start)
    self.writeButton.clicked.connect(self.write_cfg)
    self.readButton.clicked.connect(self.read_cfg)
    
    self.freqX.clicked.connect(self.read_ini)
    self.freqX.clicked.connect(self.sweep_auto)
    #self.freqX.clicked.connect(self.trigResc)

    self.singleSweep.clicked.connect(self.trigResc)
    self.autoSweep.clicked.connect(self.sweep_auto)
    self.stopSweep.clicked.connect(self.cancel)
    self.datButton.clicked.connect(self.write_dat)
    self.readDat.clicked.connect(self.read_dat)
    self.readDat.clicked.connect(self.update_tab)
    self.startValue.valueChanged.connect(self.set_start)
    self.startValue.valueChanged.connect(self.set_size)
    self.stopValue.valueChanged.connect(self.set_stop)
    self.stopValue.valueChanged.connect(self.set_size)
    self.stepValue.valueChanged.connect(self.set_size)
    self.rateValue.currentIndexChanged.connect(self.set_rate)
    self.level1Value.valueChanged.connect(self.set_level1)
    self.tabWidget.currentChanged.connect(self.update_tab)
    # create timers
    self.startTimer = QTimer(self)
    self.startTimer.timeout.connect(self.timeout)
    self.sweepTimer = QTimer(self)
    self.sweepTimer.timeout.connect(self.sweep_timeout)
    self.TstepTimer = QTimer()
    self.TstepTimer.timeout.connect((self.TconTimer)) #contols temperature control time between updates

    self.swF = True #ensures that plot curves runs only once
    self.EnableCursor = True #control writing data file
    self.stopChangeSetPoint = False # prevents set point change during open loop arduino sweep
    #operate temperature controller
    self.strT.setEnabled(True)  
    self.stpT.setEnabled(False)
    self.StopC.setEnabled(False)
    self.StopC.clicked.connect(self.TconStop)
    self.comID.setText(comport)
    self.setPt.valueChanged.connect(self.set_Temp)
    self.updtT.clicked.connect(self.read_T)
    self.sSw.clicked.connect(self.TconTimer) #starts temperature control stuff
    self.Tstart.valueChanged.connect(self.final_T)
    self.Tstep.valueChanged.connect(self.final_T)
    self.nStep.valueChanged.connect(self.final_T)
    self.Stime.valueChanged.connect(self.final_T)
    self.Tkc.setValue(id)
    self.setPt.setSuffix(sx)
    self.setPt.setValue(int(id))
    self.Treset.clicked.connect(self.Reset_Tc)
    self.saveT.clicked.connect(self.writeTout)
    self.read_T()
    self.n = 0
    self.navg = 16 #length of averaging ARRAY
    self.nrds0 = 32 #64 number of reads for getting average T
    self.nrds = 0
    self.Tbar2 = 0.0
    self.xcount = 0 #number of temperature  steps executed
    self.fstring = ''
    self.xflag = False # determines if magnified or full scan
    self.contRun = False # indicates when controller is running
    self.xlen = 0
    self.Tbar = np.zeros(self.navg)
    if sx == 'C': 
      self.deltaS.setValue(180)
      self.deltaT.setValue(1)
      self.Tstart.setValue(60)
    if sx == 'K': 
        self.deltaS.setValue(240)
        self.deltaT.setValue(0.05)
        self.Tstart.setValue(293)
    if comFlag == True:
       self.set_Tenab(True)
       self.sSw.setEnabled(True)
       self.strT.setEnabled(True)
       self.setPt.setEnabled(True)
    else:
       self.set_Tenab(False)
       self.sSw.setEnabled(False)
       self.strT.setEnabled(False)
       self.setPt.setEnabled(False)

  def Reset_Tc(self): # rests the serial port
    global comport, ser
    ser.close()
    ser = serial.Serial(comport, 9600, timeout=2)  # open serial port 
    
  def read_T(self):
    global sx
    if comFlag == False: return
    ser.write(b't\n')     # specimen temperature
    time.sleep(0.1)
    td = ser.readline()
    td = td.decode('UTF-8')
    tds = td.strip('\n')
    if tds == '': return
    tds = float(tds)
    self.Tkc.setSuffix(sx)
    self.Tkc.setValue(tds)
    self.Tstart.setSuffix(sx)
    self.Tstep.setSuffix(sx)
    self.EndT.setSuffix(sx)
    ser.write(b'p\n')     # relative power
    time.sleep(0.1)
    td = ser.readline()
    td = td.decode('UTF-8')
    id = td.strip('\n')
    if id == '': return
    id = float(id)
    if sx == 'C':id = id/10.23
    self.rP.setSuffix(' %')
    self.rP.setValue(id)
    ser.write(b's\n')     # set point
    time.sleep(0.1)
    td = ser.readline()
    td = td.decode('UTF-8')
    id = td.strip('\n')
    if id == '': return
    id = float(id)
    self.setPt.setValue(id)
    if sx == 'K': # control block temperature for Thermoelectric stage
        ser.write(b'c\n')     
        time.sleep(0.1)
        td = ser.readline()
        td = td.decode('UTF-8')
        id = td.strip('\n')
        id = float(id)
        self.bT.setSuffix(sx)
        self.bT.setValue(id)

  def set_Temp(self):   #set point
    if comFlag == False: return
    if self.stopChangeSetPoint == True: return #in arduino mode prevents changing set point during sweep
    Tset = self.setPt.value()
    if Tset < 10.0: Tset = 10.0
    Tset = str(Tset) + '\n'
    Tset = Tset.encode('UTF-8')
    ser.write(Tset)

  def final_T(self):   #computes final temperature
      if comFlag == False: return
      Tstrt = self.Tstart.value()
      Tstp = self.Tstep.value()
      nStp = self.nStep.value()
      self.EndT.setValue(Tstrt + nStp*Tstp)

  def Tsweep(self):  #sets up arduino temperature sweep
      #self.strT.setEnabled(False)
      if comFlag == False: return
      Tstrt = self.Tstart.value()
      Tstp = self.Tstep.value()
      nStp = self.nStep.value()
      sTm = self.Stime.value()
      self.EndT.setValue(Tstrt + (nStp+1)*Tstp)
      self.final_T()
      sTstrt = str(Tstrt)
      sTstp = str(Tstp)
      snStp = str(int(nStp))
      ssTm = str(int(sTm))
      self.stopChangeSetPoint = True #blocks set point change
      self.xflag = False
      scan = 'x,'+ sTstrt +','+ sTstp +','+ snStp +','+ ssTm
      tscan = scan.strip()+'\n'
      tss= tscan.encode('UTF-8')
      ser.write(tss)
      #self.TconStop()

  def TconTimer(self): # Timer for temperature stuff
    self.sSw.setEnabled(False)
    self.StopC.setEnabled(True)
    self.deltaT.setEnabled(False)
    self.deltaS.setEnabled(False)
    self.TstepTimer.start(5000)

    Tstp = self.Tstep.value()
    nStp = int(self.nStep.value())
    sTm = self.Stime.value()
    self.EnableCursor = False #blocks cursor from plotting
    self.contRun = True
    self.TconSet(Tstp, nStp,sTm)
    ###########################################################
    if self.stopChangeSetPoint == False: self.read_T()
          
  def TconStop(self):#stops temperature control
    self.TstepTimer.stop()
    self.sSw.setEnabled(True)
    self.StopC.setEnabled(False)
    self.deltaT.setEnabled(True)
    self.deltaS.setEnabled(True)
    self.EnableCursor = True
    #self.stopChangeSetPoint = True
    self.contRun = False
    self.n = 0

  def writeTout(self):
      if comFlag == False: return
      ft = os.path.join(path, sx +'_' + dt + '_T.dat')
      fh = open(ft, 'w' )
      lT =len(tout[0])
      for k in  range(lT):
        fh.write('%12.2f %12.2f \n' %(ttout[0][k], tout[0][k]) )  
      fh.close()
      #print(tout)

  def TconSet(self, Tstp, nStp,sTm):#choses between periodic sweep and step and control
    global sx
    if comFlag == False: return
    if sx =='C':
      self.navg = 16 #16 = 3.5 minutes
      Terr = 0.0
    if sx =='K':
      self.navg = 16
      Terr = self.Terror
    if sTm > 0: #use arduino mode
        self.TconStop()
        self.stopChangeSetPoint = True   
        self.Tsweep()
        return
    else: #step and wait mode
        if self.EnableCursor == False:
            self.stopChangeSetPoint = False
            ser.write(b't\n')     # get specimen temperature
            td = ser.readline()
            td = td.decode('UTF-8')
            tid = td.strip('\n')
            if tid == '': return
            tid = float(tid)
            self.nrds0 = int(self.deltaS.value()/5)
            n = self.n
            if n > nStp: #end of process
                self.setPt.setValue(self.Tstart.value())
                self.Tbar = np.zeros(self.navg)
                self.writeTout()# saves temperature log
                self.n = 0
                self.nrds = 0
                self.Tbar2 = 0.0
                self.cancel()
                self.TconStop()
                return           
            Tstrt = self.Tstart.value() + n*Tstp
            self.setPt.setValue(Tstrt)
            self.Tbar = np.delete(self.Tbar,0)
            self.Tbar = np.append(self.Tbar,[tid])
            self.Tbar2 = np.mean(self.Tbar)
            if abs(self.Tbar2-Tstrt - Terr) < self.deltaT.value(): 
                self.nrds = self.nrds + 1
                print(self.nrds,self.deltaT.value(), self.Tbar2, Terr )
            if(self.nrds   > self.nrds0):
                self.TstepTimer.stop() #stops temperature stuff to give time to write last file
                if self.Xstep.value() == 0:  self.trigResc()# saves full sweep on reaching equilibrium
                if self.Xstep.value() > 0:  # saves xscans on reaching equilibrium
                    self.read_ini()
                    self.sweep_auto()
                n = n+1
                self.n = n
                self.Tbar = np.zeros(self.navg) #resets average
                self.nrds = 0
                self.Tbar2 = 0.0
      
  def trigResc(self): # starts single sweep and sets swF flag
    mode = self.mode
    self.swF = True
    # here will switch between fullscna and xscan
    self.sweep(mode)
    self.auto = False

  def set_enabled(self, enabled):
    widgets = [self.rateValue, self.level1Value, self.startValue, self.stopValue, self.stepValue, \
      self.singleSweep, self.autoSweep, self.swpD, self.comID, self.freqX, self.Xstep]  
    for entry in widgets:
      entry.setEnabled(enabled)

  def set_Tenab(self,enabled):
    Twidg = [self.updtT, self.Tkc, self.bT, self.rP, self.Tstart, self.Tstep, self.nStep, self.Stime, \
      #self.EndT, self.sSw,  self.strT,self.stpT, self.comID, self.Treset, self.StopC, self.ClearT]
      self.EndT, self.comID, self.Treset,  self.ClearT, self.saveT, self.deltaT, self.deltaS]
    for entry in Twidg:
      entry.setEnabled(enabled)

  def start(self):
    if self.idle:
      self.connectButton.setEnabled(False)
      self.socket.connectToHost(self.addrValue.text(), 1001)
      self.startTimer.start(2000)
    else:
      self.stop()

  def stop(self):
    self.idle = True
    self.cancel()
    self.socket.abort()
    self.connectButton.setText('Connect')
    self.connectButton.setEnabled(True)
    self.set_enabled(False)
    self.stopSweep.setEnabled(False)

  def timeout(self):
    self.display_error('timeout')
##########################################################################################################
  def connected(self):
    self.startTimer.stop()
    self.idle = False
    self.set_rate(self.rateValue.currentIndex())
    corr = 0
    self.set_corr(corr)
    phase1 = 0
    self.set_phase1(phase1)
    phase2 = 0
    self.set_phase2(phase2)
    self.set_level1(self.level1Value.value())
    level2 = 0
    self.set_level2(level2)
    self.set_gpio(1)
    self.connectButton.setText('Disconnect')
    self.connectButton.setEnabled(True)
    self.set_enabled(True)
    self.stopSweep.setEnabled(True)

  def read_data(self):
    while(self.socket.bytesAvailable() > 0):
      if not self.reading:
        self.socket.readAll()
        return
      size = self.socket.bytesAvailable()
      self.progressBar.setValue((self.offset + size) / 16)
      self.progressBar.setTextVisible(True)
      limit = 16 * self.sweep_size
      if self.offset + size < limit:
        self.buffer[self.offset:self.offset + size] = self.socket.read(size)
        self.offset += size
      else:
        self.buffer[self.offset:limit] = self.socket.read(limit - self.offset)
        adc1 = self.data[0::2]
        #adc2 = self.data[1::2]
        attr = getattr(self, 'reim')
        start = self.sweep_start
        stop = self.sweep_stop 
        size = self.sweep_size
        attr.freq = np.linspace(start, stop, size)
        attr.data = adc1[0:size].copy()
        attr.unity = np.full((size),1)
        if self.EnableCursor == True and self.auto == False: self.update_tab()  #May need to fix this
        if self.EnableCursor == False: 
            if self.xflag == False: self.write_dat()
            if self.stopChangeSetPoint == True: self.TstepTimer.start(5000) #ensures in arduino no reset
        self.reading = False
        if not self.auto:
            self.progressBar.setValue(0)
            self.set_enabled(True)
        if self.xflag == False and self.contRun == True: 
            if self.stopChangeSetPoint == False: 
                self.TstepTimer.start(5000)

  def display_error(self, socketError):
    self.startTimer.stop()
    if socketError == 'timeout':
      QMessageBox.information(self, 'rus', 'Error: connection timeout.')
    else:
      QMessageBox.information(self, 'rus', 'Error: %s.' % self.socket.errorString())
    self.stop()

  def set_start(self, value):
    self.sweep_start = value
    self.set_size

  def set_stop(self, value):
    if value < self.sweep_start + 1: 
      value = value+1
      #self.stopValue.setValue(value)
    self.sweep_stop = value
    self.nSize.setValue(self.sweep_size)

  def set_size(self):
    Hz = self.stepValue.value()
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/Hz)
    if self.sweep_size > 32766:
      self.sweep_size = 32766
      Hz = int((1000*(self.sweep_stop-self.sweep_start))/32766 + 0.5)
      if Hz <=1: Hz = 1
    self.stepValue.setValue(Hz)
    self.nSize.setValue(self.sweep_size)

  def set_rate(self, value):
    if self.idle: return
    rate = [30, 100, 300, 1000][value]
    #rate = [10, 50, 100, 500, 1000, 5000, 10000, 50000][value]
    #self.rateValue.addItems(['5000', '1000', '500', '100', '50', '10', '5', '1'])
    self.socket.write(struct.pack('<I', 3<<28 | int(rate)))

  def set_corr(self, value):
    if self.idle: return
    self.socket.write(struct.pack('<I', 4<<28 | int(value & 0xfffffff)))

  def set_phase1(self, value):
    if self.idle: return
    self.socket.write(struct.pack('<I', 5<<28 | int(value)))

  def set_phase2(self, value):
    if self.idle: return
    self.socket.write(struct.pack('<I', 6<<28 | int(value)))
##############################################################################################################################
  def set_level1(self, value):
    if self.idle: return
    data = 0 if (value <=0 or value >25) else int(32766*value/25)   
    self.socket.write(struct.pack('<I', 7<<28 | int(data)))

  def set_level2(self, value):
    if self.idle: return
    data = 0 if (value <=0 or value >25) else int(32766*value/25)
    self.socket.write(struct.pack('<I', 8<<28 | int(data)))

  def set_gpio(self, value):
    if self.idle: return
    self.socket.write(struct.pack('<I', 9<<28 | int(value)))

############################################################################################################
  def sweep(self, mode):
    global ix0,ix1
    ix0 = 0
    ix1 = self.sweep_size-1
    if self.idle: return
    self.set_enabled(False)
    self.mode = mode
    self.offset = 0
    self.reading = True
    self.stepValue.valueChanged.connect(self.set_size)
    #print(self.sweep_start, self.sweep_stop, 'sweep')
    self.socket.write(struct.pack('<I', 0<<28 | int(self.sweep_start * 1000)))
    self.socket.write(struct.pack('<I', 1<<28 | int(self.sweep_stop * 1000)))
    self.socket.write(struct.pack('<I', 2<<28 | int(self.sweep_size)))
    self.socket.write(struct.pack('<I', 10<<28))
    self.progressBar.setMinimum(0)
    self.progressBar.setMaximum(self.sweep_size)
    self.progressBar.setValue(0)
    if self.auto == True and self.xflag == False: self.write_dat() # writes data file at each temperature for full scan
    #print('sweep')

  def cancel(self):
    self.sweepTimer.stop()
    self.auto = False
    self.reading = False
    self.socket.write(struct.pack('<I', 11<<28))
    self.progressBar.setValue(0)
    self.set_enabled(True)

  def sweep_auto(self):
    self.auto = True
    if self.xflag == True: dtime = 500
    else: dtime = 1000*self.swpD.value()
    self.sweepTimer.start(dtime)

  def sweep_timeout(self):
    if not self.reading:
        if self.xflag == True and self.stopChangeSetPoint == False:  #captures regions around detected peaks

            if self.xcount >= self.xlen-1: #ends x scans
                self.xflag = False
                self.sweepTimer.stop
                f = self.reim.freq
                d = self.reim.data
                for i in range(len(f)):
                    self.fd.write('%12.5f  %12.10f  %12.10f\n' % (f[i], d.real[i], d.imag[i]))#saves last scan
                self.fd.write('%12.5f  %12.10f  %12.10f\n' % (self.inir, 0.0, 0.0))#adds point at end or original sweep
                self.startValue.setValue(self.inil)
                self.stopValue.setValue(self.inir)
                self.rateValue.setCurrentIndex(self.inirate)
                self.stepValue.setValue(self.iniHz)
                self.fxn.close()
                self.fd.close()
                self.cancel()
                self.xcount = 0
                if self.setPt.value() == self.EndT.value():
                  self.setPt.setValue(self.Tstart.value())
                  self.writeTout()
                  self.TconStop()
                  self.xflag = False
                else:self.TstepTimer.start(5000) ############################################################   

            elif self.Xstep.value() > 0 and self.stopChangeSetPoint == False: # executes x scans
                self.TstepTimer.stop()
                line = self.fxn.readline()
                self.xcount = self.xcount+1
                fields = line.split()
                fx = (float(fields[0]))
                xr = self.stopValue.value() #gets last scan end point
                indx = 3
                self.set_rate(indx)
                self.rateValue.setCurrentIndex(indx)
                Hz = self.Xstep.value()
                self.sweep_Hz = Hz
                strt = 1000*fx-0.9*Hz
                if strt <= self.inil: strt = self.inil
                stpt = strt + Hz
                if self.xcount >= 2:
                    if strt < xr: strt = xr + 0.001*Hz
                if stpt > self.inir: 
                    self.inir = stpt + Hz
                self.sweep_stop = stpt
                self.sweep_start = strt
                self.startValue.setValue(strt)
                self.stopValue.setValue(stpt)
                self.stepValue.setValue(Hz)
                self.sweep_size = 1000
                self.nSize.setValue(self.sweep_size)
                self.sweep('reim')
                f = self.reim.freq
                d = self.reim.data
                for i in range(len(f)): # appends each window
                    self.fd.write('%12.5f  %12.10f  %12.10f\n' % (f[i], d.real[i], d.imag[i]))              
        else: self.sweep('reim')

  def update_tab(self):
    global comFlag
    index = self.tabWidget.currentIndex()
    self.tabs[index].update(self.graphs[index])

  def write_cfg(self):
    dialog = QFileDialog(self, 'Write configuration settings', '.', '*.ini')
    dialog.setDefaultSuffix('ini')
    dialog.selectFile('rus.ini')
    dialog.setAcceptMode(QFileDialog.AcceptSave)
    dialog.setOptions(QFileDialog.DontConfirmOverwrite)
    if dialog.exec() == QDialog.Accepted:
      name = dialog.selectedFiles()
      settings = QSettings(name[0], QSettings.IniFormat)
      self.write_cfg_settings(settings)

  def read_cfg(self):
    dialog = QFileDialog(self, 'Read configuration settings', '.', '*.ini')
    dialog.setDefaultSuffix('ini')
    dialog.selectFile('rus.ini')
    dialog.setAcceptMode(QFileDialog.AcceptOpen)
    if dialog.exec() == QDialog.Accepted:
      name = dialog.selectedFiles()
      settings = QSettings(name[0], QSettings.IniFormat)
      self.read_cfg_settings(settings)
      window.update_tab()

  def read_ini(self): #reads frequencies to expand
        global dt, fn
        id0 = self.Tkc.value()
        tid = str(int(id0+0.5))
        tdat= time.time()-ts
        self.xcount = 0
        #self.fstring = os.path.join(path, tid + '_' +sx + '_' + dt + '_F.dat') #frequencies found
        self.fstring = fn
        ###############################################################################################
        self.xstring = os.path.join(path, tid + '_' +sx + '_' + dt + '_X.dat') # saves x scans
        self.inil = self.startValue.value()
        self.inir = self.stopValue.value()
        self.inirate = self.rateValue.currentIndex()
        self.iniHz = self.stepValue.value()
        fxd = self.xstring
        fd = open(fxd, 'w') # fd is x data file and is opened as new
        self.fd = fd
        Tsample = float(self.Tkc.value())
        Vsample = self.level1Value.value()
        fd.write('%8.2f' % (Tsample))
        fd.write(' ' +sx +'\n')
        fd.write('%8.2f' % (Vsample))
        fd.write(' V' +'\n')
        fd.write('%8.2f' % (tdat))
        fd.write(' s' +'\n')
        fd.write('      frequency    real         imag\n')
        fd.write('%12.5f  %12.10f  %12.10f\n' % (self.inil, 0.0, 0.0)) #writes start of original scan
        self.fd.close()
        fd = open(fxd, 'a') #now open x file in append mode
        self.fd = fd
        self.endF = self.reim.freq[-1]
        self.endR = self.reim.data.real[-1]
        self.endI = self.reim.data.imag[-1]
        self.reim.freq =  []
        self.reim.data = []
        if os.path.exists(self.fstring):
          fx = self.fstring
          fxn = open(fx, 'r') #fxn is list of peaks
          self.xflag = True
          self.xlen = int(len(fxn.readlines())) #gets file length
          fxn.seek(0) #resets read pointer
          self.fxn = fxn
          xT = fxn.readline() #reads temperature and tosses it
        else: return
     

  def write_cfg_settings(self, settings):
    settings.setValue('addr', self.addrValue.text())
    settings.setValue('rate', self.rateValue.currentIndex())
    settings.setValue('level_1', self.level1Value.value())
    settings.setValue('reim_start', int(self.reim.freq[0]))
    settings.setValue('reim_stop', int(self.reim.freq[-1]))
    #settings.setValue('reim_size', self.reim.freq.size)
    settings.setValue('step',self.sweep_Hz)
  
  def read_cfg_settings(self, settings):
    self.addrValue.setText(settings.value('addr', '192.168.1.100'))
    self.rateValue.setCurrentIndex(settings.value('rate', 0, type = int))
    self.level1Value.setValue(settings.value('level_1', 1, type = float))
    reim_start = settings.value('reim_start',200, type = int)
    reim_stop = settings.value('reim_stop', 450, type = int)
    self.sweep_Hz = settings.value('step',50,type = int)
    self.Terror = settings.value('Terr', -0.13, type = float)
    self.startValue.setValue(reim_start)
    self.stopValue.setValue(reim_stop)
    self.stepValue.setValue(self.sweep_Hz)

  def read_dat(self): #read data from file
    dialog = QFileDialog(self, 'Read dat file', '.', '*_*.dat')
    dialog.setDefaultSuffix('dat')
    dialog.setAcceptMode(QFileDialog.AcceptSave)
    dialog.setOptions(QFileDialog.DontConfirmOverwrite)
    if dialog.exec() == QDialog.Accepted: 
      name = dialog.selectedFiles()
      fh = open(name[0], 'r')
      id = fh.readline()
      id = id.strip('\n')
      id = id[0:-2]
      fid = float(id)
      vd = fh.readline()
      vd = vd.strip('\n')
      vd = vd[0:-2]
      vfd = float(vd)
      tdat = fh.readline()
      #print(vfd)
      self.level1Value.setValue(vfd)#################################
      self.Tkc.setSuffix(sx)
      self.Tkc.setValue(fid)
      hd = fh.readline()
      freq, real, imag, = [], [], []
      for line in fh.readlines():
        fields = line.split()
        freq.append(float(fields[0]))
        real.append(float(fields[1]))
        imag.append(float(fields[2]))
      l = len(freq)
      data = np.zeros(l, np.complex64)
      data.real = real
      data.imag = imag
      self.reim.freq = freq
      self.reim.data = data
      fh.close 

  def write_dat(self): #write data to file
    global sx, comFlag,dt,ts
    trds = '0.0'
    tflt = 0.0
    tdat= time.time()-ts
    if self.auto == True:
      t = datetime.now()
      dat = t.strftime('%d%m%Y%H%M%S')
      if comFlag == True:
        ser.write(b't\n')     # get specimen temperature
        tread = ser.readline()
        treads = tread.decode('UTF-8')
        trds = treads.strip('\r\n')
        tflt = float(trds)
      fd = os.path.join(path, trds + '_' +sx +'_' + dat + '_D.dat')
    else:
      dat = dt
      if comFlag == True:
        ser.write(b's\n')     # get set point temperature
        tread = ser.readline()
        treads = tread.decode('UTF-8')
        trds = treads.strip('\r\n')
        tflt = float(trds)
        tid = str(tflt)
        if self.Tbar2 != 0 :tid = str(self.Tbar2)
        fd = os.path.join(path, tid + '_' +sx +'_' + dat + '_D.dat')
      else:
        fd = os.path.join(path, trds + '_' +sx +'_' + dat + '_D.dat')
    f = self.reim.freq
    d = self.reim.data
    Tsample = float(self.Tkc.value())
    Vsample = self.level1Value.value()
    if (d.real[0]) and (d.imag[0]) != 0.0:
        fh = open(fd, 'w' )
        fh.write('%8.2f' % (Tsample))
        fh.write(' ' +sx +'\n')
        fh.write('%8.2f' % (Vsample))
        fh.write(' V' +'\n')
        fh.write('%8.2f' % (tdat))
        fh.write(' s' +'\n')
        fh.write('      frequency    real         imag\n')
        for i in range(len(f)):
            fh.write('%12.5f  %12.10f  %12.10f\n' % (f[i], d.real[i], d.imag[i]))     
        fh.close()
        if self.xflag == True: self.TstepTimer.start(5000)

#Main control code follows
warnings.filterwarnings('ignore')
app = QApplication(sys.argv)
splash_pix = QPixmap('Splash.png')
splash = QSplashScreen(splash_pix)
splash.show()
#time.sleep(1)
app.processEvents()

path = os.getcwd()
path = os.path.join(path, 'RUS_Data')
if os.path.exists(path) == False: os.makedirs(path)
sqr2 = 1/math.sqrt(2)
z = 0.0
ts = time.time()
t = datetime.now()
dt = t.strftime('%d%m%Y%H%M%S')
id = 0.0
comFlag = True
comport = ''
sx = ''
aPf = True
#tstest = 350
ports = serial.tools.list_ports.comports()
if ports ==[]:
    comFlag = False
if comFlag == True:
    p = ports[0]
    comport = p.device
    ser = serial.Serial(comport, 9600, bytesize=8, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, timeout=2, rtscts=0)  # open serial port 
    #ser = serial.Serial(comport, 9600, timeout=2)  # open serial port 
    time.sleep(3)
    ts = time.time()
    tout = []
    ttout = []
    tout.append([])
    ttout.append([])
    l = ser.readline()

    ser.write(b'*IDN?\n')     # get id
    l2 = ser.readline()
    lid = l2.decode('UTF-8')
    lid = lid.strip('\n')

    ser.write(b't\n')     # get temperature
    td = ser.readline()
    td0 = td.decode('UTF-8')
    id = td0.strip('\n')
    id = float(id)

    ser.write(b'm\n')     # get controller type
    md = ser.readline()
    md0 = md.decode('UTF-8')
    #sx = str(md0.strip('\n'))
    sx = str(md0.strip())
    tout[-1].append(id)
    ts1 = time.time()-ts
    ttout[-1].append(ts1)

window = rus()
window.update_tab()
window.show()
splash.finish(window)
sys.exit(app.exec_())
