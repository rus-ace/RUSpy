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
from ctypes.wintypes import SC_HANDLE
from operator import pos
from datetime import datetime
#from os import wait
from pickle import NONE, TRUE
from sqlite3 import Row
from statistics import mode
import sys
import struct
from tracemalloc import stop
from turtle import onrelease
from unittest.mock import NonCallableMock
import warnings
from matplotlib import pyplot as plt

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
from PyQt5.QtGui import QRegExpValidator, QPalette, QColor
from PyQt5.QtWidgets import QApplication, QDoubleSpinBox, QMainWindow, QMessageBox, QDialog, QFileDialog, QPushButton, QLabel, QSpinBox
from PyQt5.QtNetwork import QAbstractSocket, QTcpSocket
from PyQt5.QtWidgets import QWidget, QVBoxLayout

sqr2 = 1/math.sqrt(2)
z = 0.0
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
    #print('no ports')
    comFlag = False
if comFlag == True:
    p = ports[0]
    comport = p.device
    #ser = serial.Serial(comport, 9600, bytesize=8, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, timeout=1, rtscts=0)  # open serial port 
    ser = serial.Serial(comport, 9600, timeout=2)  # open serial port 
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

Ui_rus, QMainWindow = loadUiType('rus.ui')

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
    self.figure = Figure()
    if sys.platform != 'win32':
      self.figure.set_facecolor('none')
    self.canvas = FigureCanvas(self.figure)
    self.canvas.mpl_connect('button_press_event', self.onpress)
    self.canvas.mpl_connect('button_release_event', self.onrelease)
    layout.addWidget(self.canvas)
    self.toolbar = NavigationToolbar(self.canvas, None, False)
    self.toolbar.layout().setSpacing(6)
    #remove subplots action
    actions = self.toolbar.actions()
    layout.addWidget(self.toolbar)
    self.tf = ''
    ##################################################################################
    #add set phase
    self.mode = 'reim'
    self.phase = QDoubleSpinBox()
    self.phase.setSingleStep(5)
    self.phase.setMaximum(180)
    self.phase.setMinimum(-180)
    self.phase.setDecimals(1)
    self.phase.valueChanged.connect(self.phase_set)
    self.toolbar.addWidget(self.phase)
    #add  best phase button
    self.bPhase = QPushButton('Best Phase')
    self.bPhase.clicked.connect(self.phase_adj)
    self.toolbar.addWidget(self.bPhase)
    #add value at crsor
    self.temp = QDoubleSpinBox()
    self.temp.setReadOnly(True)
    self.temp.setButtonSymbols(2)
    self.temp.setDecimals(6)
    self.toolbar.addWidget(self.temp)
    #add crsor
    self.crsor = QDoubleSpinBox()
    self.crsor.setDecimals(3)
    self.crsor.setStyleSheet('background-color:pink')
    self.crsor.valueChanged.connect(self.plot_crsr)
    self.toolbar.addWidget(self.crsor)
    #add mark peak button
    self.peak = QPushButton('Mark')
    self.peak.setStyleSheet('background-color:orange')
    self.peak.clicked.connect(self.mark)
    self.peak.clicked.connect(self.wFreq)
    self.toolbar.addWidget(self.peak)
    layout.addWidget(self.toolbar)
    #add unmark button
    self.unpeak = QPushButton('unMark')
    self.unpeak.setStyleSheet('background-color:red; color:rgb(255, 255, 255)')
    self.unpeak.clicked.connect(self.unmark)
    #self.unpeak.clicked.connect(self.wFreq)
    self.toolbar.addWidget(self.unpeak)
    layout.addWidget(self.toolbar)
    #add rescale button
    self.plotButton = QPushButton('Replot')
    self.plotButton.setStyleSheet('background-color:black; color:rgb(255, 255, 255)')
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
    self.ix0 = 0
    self.ix1 = len(self.rus.reim.freq)
    self.izr = 0
    self.izl = self.ix1
    global rusf
    rusf = []
    self.TsweepTimer = QTimer()
    self.TsweepTimer.timeout.connect(self.Tcont)
    ###############################################################################
    self.rus.strT.clicked.connect(self.Tcont) #plots temperature
    self.rus.stpT.clicked.connect(self.TcontStop)
    self.rus.CM.clicked.connect(self.clearrusf)
    self.rus.AP.clicked.connect(self.autoPeak)
    self.rus.AP.clicked.connect(self.reset)
    self.rus.ClearT.clicked.connect(self.tsReset)
    #self.cur = Cursor(axes2, horizOn=True, vertOn=True, useblit=True, color = 'blue', linewidth = 1)
    #self.canvas.mpl_connect('motion_notify_event', self.onEnter)
  
  def tsReset(self): #resets time axis for temperature plot
    global ts, tout,ttout
    ts = time.time()
    self.axes2.cla()
    tout = []
    ttout = []
    tout.append([])
    ttout.append([])

  def unmark(self): #unmarks peaks in window
    global rusf
    if rusf == []:return
    izl = self.zbox()[0]
    izr = self.zbox()[1]
    fzl = self.rus.reim.freq[izl]
    fzr = self.rus.reim.freq[izr]
    rows = len(rusf)
    rusa = rusf
    rusf = []
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
    self.update(self.mode)
    return

  def Tcont(self): #Starts the continuous plot of temperature
        self.rus.strT.setEnabled(False)
        self.rus.stpT.setEnabled(True)
        if self.rus.Tsw == False: self.axes2.cla()
        self.TsweepTimer.start(1000)
        self.plot_tout()#plots a temperature point each 1000ms

  def TcontStop(self): #stops continuous temperature plots
    self.TsweepTimer.stop()
    self.rus.strT.setEnabled(True)
    #self.update(self.mode)

  def clearrusf(self): #clears all peaks (stired in rusf)
    global rusf
    rusf = []
    self.update(self.mode)

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
    for q in range(0, l):
      qbar = qbar +pdata[q]
    qbar = qbar/l
    qrms = 0.0
    for q in (0, l):
      qrms = qrms + (pdata[q]-qbar)**2
    qrms = math.sqrt(qrms/l)
    #print(qrms)
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
      if ((y2 > y1) and (y2 >y3) and (y3 > y3p+m) and (y1 > y1p+m)):
        if (y2 > 0.01*percent*pmax):
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
          dat = dt
          tid = str(self.rus.setPt.value())
          fn = tid + '_' +sx +'_' + dat + '_F.dat'
          b = sorted(rusf)
          bx = [i[0] for i in b]
          by = [i[1] for i in b]
          if bx != []:
            with open(fn, 'w', newline='') as f:
                f.write('%8.2f' % id)
                f.write(' ' +sx +'\n')
                writer = csv.writer(f, delimiter = ' ')
                writer.writerows(b)
                self.axes1.plot(bx,by,color = 'black', marker = matplotlib.markers.CARETDOWNBASE , markersize =10, linestyle = 'None')
                self.canvas.draw()
                self.update(mode)
  
  def plot_curves(self, izl, izr, xhgt, tf): #main plotting routine
    matplotlib.rcdefaults()
    matplotlib.rcParams['axes.formatter.use_mathtext'] = True
    self.figure.clf()
    self.figure.text(0.04, 0.005,r'$\phi$                    f          Z', fontsize = 18)
    freqP = self.rus.reim.freq
    sstep = (freqP[1]-freqP[0])
    dataP = self.rus.reim.data
    swF = self.rus.swF
    if swF == True:
        izl = 0
        izr = self.rus.sweep_size
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
    y = 1.05*np.max(data3)
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
        if tf == 'adj':
            ds1 = sum(data1)
            ds2 = sum(data2)
            dsq = math.sqrt(ds1*ds1+ds2*ds2)
            co = (ds1/dsq)
            si = abs(ds2/dsq)
            ztheta = math.acos(co)
            zt = math.degrees(ztheta)
            self.phase.setValue(zt)
            self.zt = zt
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
    axes1.set_xlabel('kHz')
    axes1.set_ylabel('real')
    #Setup crsor spinbox
    self.crsor.setMaximum(freqP[-1])
    self.crsor.setDecimals(3)
    self.crsor.setMinimum(freqP[0])
    self.crsor.setSingleStep(sstep)##########################3
    self.crsor.setValue(freqmax)
    #setup axes
    hgt = data3[i]
    xlim = [freqP[0],freqP[-1]]
    axes1.set_xlim(xlim)
    axes1.tick_params('y', color = 'black', labelcolor = 'black')
    axes1.yaxis.label.set_color('blue')
    axes1.set_ylabel('imag')
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
        axes1.set_ylabel('Volts')
        axes1.yaxis.label.set_color('black')
        axes1.plot(freqP, d1, color = 'blue')
        axes1.plot(freqP, d2, color = 'red')  ##
    if self.mode == 'absm': 
        self.temp.setValue(np.abs(data3[i]))
        axes1.set_ylabel('Volts')
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
    self.plot_crsr()
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

  def wFreq(self):#writes found peaks to file
    global rusf, sx, comFlag, dt, id
    dat = dt
    id0 = self.rus.setPt.value()
    tid = str(self.rus.setPt.value())
    fn = tid + '_' +sx +'_' + dat + '_F.dat'
    b = sorted(rusf)
    bx = [i[0] for i in b]
    by = [i[1] for i in b]
    if bx != []:
      with open(fn, 'w', newline='') as f:
          f.write('%8.2f' % id0)
          f.write(' ' +sx +'\n')
          writer = csv.writer(f, delimiter = ' ')
          writer.writerows(b)
          self.axes1.plot(bx,by,color = 'black', marker = matplotlib.markers.CARETDOWNBASE , markersize =10, linestyle = 'None')
          self.canvas.draw()

  def update(self, mode):
    getattr(self, 'update_%s' % mode)()

  def reset(self):#resets plot
    izl = 0
    izr = len(self.rus.reim.data)
    self.ix0 = 0
    self.ix1 = len(self.rus.reim.data)
    self.izl = 0
    self.izr = self.ix0
    tf = ''
    xhgt = self.zbox()[3]
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
        if izr == izl: return
        if izl > izr: return
        lim = [izl,izr]
        self.axes1.set_xlim(lim)
        xhgt = self.zbox()[3]
        tf = self.tf
        self.plot_curves(izl, izr, xhgt, tf)
        self.wFreq()
        return izl, izr, xhgt, tf

  def plot_crsr(self):#plots vertical line at peak position and can be moved
    if self.rus.Tsw == True:
        vc = self.crsor.value()
        #axes1 = self.axes1
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

  def phase_adj(self):#roll the phase
    if self.mode =='absm': return
    tf = 'adj'
    self.tf = tf
    getattr(self, 'plot_%s' % self.mode)()
  

  def phase_set(self): #determines best phase
      if self.mode == 'absm': return
      tf = 'set'
      self.tf = tf
      getattr(self, 'plot_%s' % self.mode)()

  def plot_reim(self):# plots real and imginary components
    izl = self.ix0
    izr = self.ix1
    xhgt = self.zbox()[3]
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()


  def update_reim(self):
      self.plot_reim()
      #getattr(self, 'plot_%s' % self.mode)()
  
  def plot_tout(self): #continuous plot of temperaure
    global tout, ttout, comFlag, ts, sx
    #self.axes2.cla()
    self.axes2.yaxis.set_ticklabels([])
    self.axes2.tick_params(right = False)
    time.sleep(0.2)
    if comFlag == False: 
        return
    ser.write(b't\n')     # write a string
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
        self.axes1.tick_params(axis= 'x', colors = 'black')
        self.axes1.set_xlabel('Time(s)')
        self.axes1.set_ylabel('Temperature('+sx+')')     
        self.axes1.set_xlim([0,ts1])
        self.axes1.set_ylim(ylim)
        #self.axes1.plot(ts1,id,color = 'black', marker = 'o' , markersize =2, linestyle = 'None')
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

  def plot_absm(self): #plots absolte value
    izl = self.ix0
    izr = self.ix1
    xhgt = self.zbox()[3]
    self.mode = 'absm'
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()
    ##############################################################################################################

  def update_absm(self):
     self.plot_absm()

  def zbox(self): #finds coordinates of freq range selected with cursor
    span = self.axes1.get_xlim()
    izl = self.ix0
    izr = self.ix1
    sz = izr-izl
    km = 1
    if span[0] == 0:
      return izl, izr, km, span
    f = self.rus.reim.freq
    d = np.absolute(self.rus.reim.data[izl:izr])
    xhgt = [-1.05*np.max(d),1.05*np.max(d)]
    f0 = span[0]
    for i in range(sz-1):
        if f[i] >= f0:
            izl = i
            break
    f1 = span[1]-1
    for j in range(sz-2):
        if f[j] >= f1:
            izr = j
            break
    km = np.argmax(d)+izl
    izr = izr + 1
    return izl, izr, km, xhgt

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
    super(rus, self).__init__()
    self.setupUi(self)
    global comFlag, comport
    # address validator
    #rx = QRegExp('^(([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\.){3}([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])|rp-[0-9A-Fa-f]{6}\.local$')
    #self.addrValue.setValidator(QRegExpValidator(rx, self.addrValue))
    # state variables
    self.idle = True
    self.reading = False
    self.auto = False
    # sweep parameters
    self.sweep_start = 300
    self.sweep_stop = 400
    self.sweep_Hz = 100
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/self.sweep_Hz)
    if(self.sweep_size)>32766:
      self.sweep_size = 32766
      self.sweep_Hz = 1000*(self.sweep_stop-self.sweep_start)/32766
    # buffer and offset for the incoming samples
    self.buffer = bytearray(16 * 32768)
    self.offset = 0
    self.data = np.frombuffer(self.buffer, np.complex64)
    # create measurements
    self.reim = Measurement(self.sweep_start, self.sweep_stop, self.sweep_size)
    ######################################################################################################
    self.mode = 'reim'  
    # create figures
    self.tabs = {}
    for i in range(len(self.graphs)):
      layout = getattr(self, '%sLayout' % self.graphs[i])
      self.tabs[i] = FigureTab(layout, self)
    # configure widgets
    self.rateValue.addItems([ '0.7','2', '6', '20', '60', '200', '600'])
    #self.rateValue.addItems(['5000', '1000', '500', '100', '50', '10', '5', '1'])
    #rate = [10, 50, 100, 500, 1000, 5000, 10000, 50000][value]
    self.rateValue.lineEdit().setReadOnly(True)
    self.rateValue.lineEdit().setAlignment(Qt.AlignRight)
    for i in range(self.rateValue.count()):
      self.rateValue.setItemData(i, Qt.AlignRight, Qt.TextAlignmentRole)
    self.set_enabled(False)
    self.stopSweep.setEnabled(False)
    # read settings
    settings = QSettings('rus.ini', QSettings.IniFormat)
    self.read_cfg_settings(settings)
    # create TCP socket
    self.socket = QTcpSocket(self)
    self.socket.connected.connect(self.connected)
    self.socket.readyRead.connect(self.read_data)
    self.socket.error.connect(self.display_error)
    # connect signals from widgets
    self.connectButton.clicked.connect(self.start)
    self.writeButton.clicked.connect(self.write_cfg)
    self.readButton.clicked.connect(self.read_cfg)
    #self.singleSweep.clicked.connect(self.sweep)
    self.singleSweep.clicked.connect(self.trigResc)
    self.autoSweep.clicked.connect(self.sweep_auto)######################################################################
    self.stopSweep.clicked.connect(self.cancel)
    self.datButton.clicked.connect(self.write_dat)
    self.readDat.clicked.connect(self.read_dat)
    self.readDat.clicked.connect(self.update_tab)
    self.startValue.valueChanged.connect(self.set_start)
    self.stopValue.valueChanged.connect(self.set_stop)
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
    self.TstepTimer.timeout.connect((self.TconTimer))
    self.swF = True #ensures that plot curves runs only once
    self.Tsw = True #blocks cursor from plotting during step T
    #operate temperature controller
    self.setPt.setEnabled(True)
    self.StopC.clicked.connect(self.TconStop)
    self.comID.setText(comport)
    self.setPt.valueChanged.connect(self.set_Temp)
    self.updtT.clicked.connect(self.read_T)
    #self.sSw.clicked.connect(self.tsReset)
    self.sSw.clicked.connect(self.TconTimer)
    #self.ClearT.clicked.connect(self.tsReset)
    self.Tstart.valueChanged.connect(self.final_T)
    self.Tstep.valueChanged.connect(self.final_T)
    self.nStep.valueChanged.connect(self.final_T)
    self.Stime.valueChanged.connect(self.final_T)
    self.Tkc.setValue(id)
    self.setPt.setSuffix(sx)
    self.setPt.setValue(int(id))
    self.Treset.clicked.connect(self.Reset_Tc)
    self.read_T()
    self.n = 0
    self.navg = 8
    self.Tbar = np.zeros(self.navg)
    if comFlag == True:
       self.set_Tenab(True)
    else:
       self.set_Tenab(False)

  def Reset_Tc(self): # rests the serial port
    global comport, ser
    ser.close()
    ser = serial.Serial(comport, 9600, timeout=2)  # open serial port 
    
  def read_T(self):
    global tout,sx, ts
    if comFlag == False: return
    ser.write(b't\n')     # specimen temperature
    #time.sleep(0.1)
    td = ser.readline()
    td = td.decode('UTF-8')
    id = td.strip('\n')
    if id == '': return
    id = float(id)
    self.Tkc.setSuffix(sx)
    self.Tkc.setValue(id)
    self.Tstart.setSuffix(sx)
    self.Tstep.setSuffix(sx)
    self.EndT.setSuffix(sx)
    ser.write(b'p\n')     # relative power
    #time.sleep(0.1)
    td = ser.readline()
    td = td.decode('UTF-8')
    id = td.strip('\n')
    if id == '': return
    id = float(id)
    self.rP.setSuffix(' %')
    self.rP.setValue(id)
    if sx == 'K': # control block temperature for Thermoelectric stage
        ser.write(b'c\n')     
        #time.sleep(0.1)
        td = ser.readline()
        td = td.decode('UTF-8')
        id = td.strip('\n')
        id = float(id)
        self.bT.setSuffix(sx)
        self.bT.setValue(id)

  def set_Temp(self):   #set point
    if comFlag == False: return
    Tset = self.setPt.value()
    if Tset < 10.0: Tset = 10.0
    Tset = str(Tset) + '\n'
    Tset = Tset.encode('UTF-8')
    ser.write(Tset)
    self.read_T()

  def final_T(self):   #computes final temperature
      if comFlag == False: return
      Tstrt = self.Tstart.value()
      Tstp = self.Tstep.value()
      nStp = self.nStep.value()
      self.EndT.setValue(Tstrt + nStp*Tstp)

  def Tsweep(self):  #sets up temperature sweep
      #self.strT.setEnabled(False)
      if comFlag == False: return
      ser.flushInput()
      time.sleep(0.2)
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
      scan = 'x,'+ sTstrt +','+ sTstp +','+ snStp +','+ ssTm
      tscan = scan.strip()+'\n'
      tss= tscan.encode('UTF-8')
      #print(tss)
      ser.write(tss)
################################################################################################################

  def TconTimer(self): # Timer for step control
    self.sSw.setEnabled(False)
    self.TstepTimer.start(5000)
    Tstp = self.Tstep.value()
    nStp = int(self.nStep.value())
    sTm = self.Stime.value()
    self.TconSet(Tstp, nStp,sTm)
    self.Tsw = False
  
  def TconStop(self):#stops temperature control
    self.TstepTimer.stop()
    self.sSw.setEnabled(True)
    self.Tsw = True
    #self.update(self.mode)

  def TconSet(self, Tstp, nStp,sTm):#choses between periodic sweep and step and control
    global sx
    if comFlag == False: return
    if sx =='C':self.navg = 16
    if sx =='K':self.navg = 32
    if sTm > 0: 
        self.TstepTimer.stop()
        self.Tsweep()
        return
    else:
        if self.Tsw == False:
            ser.write(b't\n')     # specimen temperature
            td = ser.readline()
            td = td.decode('UTF-8')
            tid = td.strip('\n')
            if tid == '': return
            tid = float(tid)
            n = self.n
            Tstrt = self.Tstart.value() +n*Tstp
            self.setPt.setValue(Tstrt)
            Tbar0 = np.sum(self.Tbar)
            self.Tbar = np.delete(self.Tbar,0)
            self.Tbar = np.append(self.Tbar,[tid])
            Tbar1 = np.sum(self.Tbar)
            Tbar2 = np.mean(self.Tbar)
            dT1 = 0.1
            dT2 = 0.3
            test1 = 999
            test2 = 376
            #print(Tbar2, Tstrt)
            if sx =='C':
                dT1 = 0.1
                test1  = abs(Tbar0-Tbar1)
            if sx =='K':
                dT2 =0.3
                test2 = abs(Tbar2 - Tstrt)
            if((test1 <= 0.02)and (test2<dT2)) or(test1 <= dT1):
                self.TstepTimer.stop()
                self.trigResc()
                n = n+1
                self.n = n
                self.Tbar = np.zeros(self.navg)
            self.Tkc.setSuffix(sx)
            self.Tkc.setValue(Tbar2)
            if n >= nStp:
                self.TstepTimer.stop() 
                self.n=n
                self.Tsw = True
 
########################################################################################################
      
  def trigResc(self): # starts single sweep and sets swF flag
    mode = self.mode
    self.swF = True
    self.sweep(mode)
    self.auto = False

  def set_enabled(self, enabled):
    widgets = [self.rateValue, self.level1Value, self.startValue, self.stopValue, self.stepValue, \
      self.singleSweep, self.autoSweep, self.swpD, self.comID]  
    for entry in widgets:
      entry.setEnabled(enabled)

  def set_Tenab(self,enabled):
    Twidg = [self.updtT, self.Tkc, self.bT, self.rP, self.Tstart, self.Tstep, self.nStep, self.Stime, \
      self.EndT, self.sSw,  self.strT,self.stpT, self.comID, self.Treset, self.StopC, self.ClearT]
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
        if self.Tsw == True: self.update_tab()
        if self.Tsw == False: 
            self.write_dat()
            self.TstepTimer.start()
        self.reading = False
        if not self.auto:
            self.progressBar.setValue(0)
            self.set_enabled(True)

  def display_error(self, socketError):
    self.startTimer.stop()
    if socketError == 'timeout':
      QMessageBox.information(self, 'rus', 'Error: connection timeout.')
    else:
      QMessageBox.information(self, 'rus', 'Error: %s.' % self.socket.errorString())
    self.stop()

  def set_start(self, value):
    self.sweep_start = value

  def set_stop(self, value):
    self.sweep_stop = value

  def set_size(self):
    Hz = self.stepValue.value()
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/Hz)
    if self.sweep_size > 32766:
      self.sweep_size = 32766
      Hz = int((1000*(self.sweep_stop-self.sweep_start))/32766 + 0.5)
      if Hz <=1: Hz = 1
      self.stepValue.setValue(Hz)

  def set_rate(self, value):
    if self.idle: return
    rate = [30, 100, 300, 1000, 3000, 10000, 30000][value]
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
    if self.idle: return
    self.set_enabled(False)
    self.mode = mode
    self.offset = 0
    self.reading = True
    self.stepValue.valueChanged.connect(self.set_size)
    self.socket.write(struct.pack('<I', 0<<28 | int(self.sweep_start * 1000)))
    self.socket.write(struct.pack('<I', 1<<28 | int(self.sweep_stop * 1000)))
    self.socket.write(struct.pack('<I', 2<<28 | int(self.sweep_size)))
    self.socket.write(struct.pack('<I', 10<<28))
    self.progressBar.setMinimum(0)
    self.progressBar.setMaximum(self.sweep_size)
    self.progressBar.setValue(0)
    if self.auto == True: self.write_dat() ################################################################

  def cancel(self):
    self.sweepTimer.stop()
    self.auto = False
    self.reading = False
    self.socket.write(struct.pack('<I', 11<<28))
    self.progressBar.setValue(0)
    self.set_enabled(True)

  def sweep_auto(self):
    self.auto = True
    dtime = 1000*self.swpD.value()
    self.sweepTimer.start(dtime)

  def sweep_timeout(self):
    if not self.reading:
      self.sweep('reim')

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
    self.rateValue.setCurrentIndex(settings.value('rate', 2, type = int))
    self.level1Value.setValue(settings.value('level_1', 1, type = int))
    reim_start = settings.value('reim_start',300, type = int)
    reim_stop = settings.value('reim_stop', 400, type = int)
    #reim_size = settings.value('reim_size', 1, type = int)
    self.sweep_Hz = settings.value('step',50,type = int)
    self.startValue.setValue(reim_start)
    self.stopValue.setValue(reim_stop)
    self.stepValue.setValue(self.sweep_Hz)

  def read_dat(self): #read data from file
    dialog = QFileDialog(self, 'Read dat file', '.', '*_D.dat')
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
      self.Tkc.setSuffix(sx)
      self.Tkc.setValue(fid)
      self.setPt.setSuffix(sx)
      self.setPt.setValue(fid)
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
    global sx, comFlag,dt
    trds = '0.0'
    tflt = 0.0
    if self.auto == True:
      t = datetime.now()
      dat = t.strftime('%d%m%Y%H%M%S')
      if comFlag == True:
        ser.write(b's\n')     # get set point temperature
        tread = ser.readline()
        treads = tread.decode('UTF-8')
        trds = treads.strip('\r\n')
        tflt = float(trds)
      fn = trds + '_' +sx +'_' + dat + '_D.dat'
    else:
      dat = dt
      if comFlag == True:
        ser.write(b's\n')     # get set point temperature
        tread = ser.readline()
        treads = tread.decode('UTF-8')
        trds = treads.strip('\r\n')
        tflt = float(trds)
        tid = str(tflt)
        fn = tid + '_' +sx +'_' + dat + '_D.dat'
      else:
        fn = trds + '_' +sx +'_' + dat + '_D.dat'
    f = self.reim.freq
    d = self.reim.data
    Tsample = float(self.Tkc.value())
    if (d.real[0]) and (d.imag[0]) != 0.0:
        fh = open(fn, 'w' )
        fh.write('%8.2f' % Tsample)
        fh.write(' ' +sx +'\n')
        fh.write('      frequency    real         imag\n')
        for i in range(f.size):
            fh.write('%12.2f  %12.7f  %12.7f\n' % (f[i], d.real[i], d.imag[i]))     
        fh.close()

warnings.filterwarnings('ignore')
app = QApplication(sys.argv)
window = rus()
window.update_tab()
window.show()
sys.exit(app.exec_())