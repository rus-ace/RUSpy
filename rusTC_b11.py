#RUS Alamo Creek Engineering Data Acquisition code ver. 3.0
# #!/usr/bin/env python3.11
# Control program for the ACE Resonant Ultrasound Spectrometer
# Ref: Reviews of Scientific Instruments Vol.90 Issue 12 Dec 2019 pgs. 121401 ff.
# Copyright (C) 2022  Albert Migliori
# Using pieces from the control program for the Red Pitaya vector network analyzer
# Copyright (C) 2021  Pavel Demin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

from asyncio import sleep
from operator import pos
from datetime import datetime
import sys
import struct
import warnings
import matplotlib
from matplotlib import pyplot as plt
import os
from pathlib import Path
import serial
import serial.tools.list_ports
import time

from functools import partial
from matplotlib.backend_bases import Event

import numpy as np
import math
import csv

from numpy.core.arrayprint import str_format
from numpy.lib.type_check import imag, real

from scipy.optimize import curve_fit

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import Formatter, FuncFormatter
from matplotlib.widgets import Cursor, MultiCursor
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

from PyQt5.uic import loadUiType
from PyQt5.QtCore import QRegExp, QTimer, QSettings, QDir, Qt
#from PyQt5.QtCore import Null, QRegExp, QTimer, QSettings, QDir, Qt
from PyQt5.QtGui import QRegExpValidator, QPalette, QColor, QBitmap, QPixmap
from PyQt5.QtWidgets import QApplication, QDoubleSpinBox, QMainWindow, QMessageBox, QDialog, QFileDialog, QPushButton, QLabel, QSpinBox
from PyQt5.QtNetwork import QAbstractSocket, QTcpSocket
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSplashScreen
import pyqtgraph as pg
pg.setConfigOption('background','w')
pg.setConfigOption('foreground', 'k')

Ui_rus, QMainWindow = loadUiType('rusTb02.ui')

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
    global comFlag, conTrunning, ix0, ix1, rusf
    self.rus = rus
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
    self.mode = 'reim'
    if comFlag == True:
      self.rus.read_T() #this reads and sets set point
    else: self.tidF = '0.0'
    #add set phase
    self.phase = QDoubleSpinBox()
    self.phase.setFixedHeight(40)
    self.phase.setSingleStep(5)
    self.phase.setMaximum(180)
    self.phase.setMinimum(-180)
    self.phase.setDecimals(1)
    self.phase.valueChanged.connect(self.phase_set)
    self.toolbar.addWidget(self.phase)
    #add height value at black crsor
    self.temp = QDoubleSpinBox()
    self.temp.setFixedHeight(50)
    self.temp.setReadOnly(True)
    self.temp.setButtonSymbols(2)
    self.temp.setDecimals(4)
    self.temp.setMinimum(-25)
    self.temp.setStyleSheet('background-color:gray;color:white')
    self.toolbar.addWidget(self.temp)
    #add frequency value of black crsor
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
    self.plotButton.clicked.connect(self.reset)
    self.plotButton.clicked.connect(self.wFreq) #writes freq data
    self.toolbar.addWidget(self.plotButton)
    axes1 = self.figure.add_subplot(1,1,1)
    axes2 = self.figure.add_subplot(1,1,1)
    self.axes1 = axes1
    self.axes2 = axes2
    self.zt = 0.0
    ix0 = 0
    ix1 = len(self.rus.reim.freq)
    rusf = []
    self.rus.CM.clicked.connect(self.clearrusf) #CM = clear marks
    self.rus.CM.clicked.connect(self.wFreq)####################################
    self.rus.AP.clicked.connect(self.autoPeak) #AP = auto peak finder
    self.rus.AP.clicked.connect(self.reset)
    self.rus.w0g.clicked.connect(self.w0_guess)
    self.rus.f1g.clicked.connect(self.f1_guess)
    popt = list()
    self.rus.save.clicked.connect(self.save_Lrnttz)
    self.Lresult = list()
    self.popt = list()
    self.niter = 0
    self.br = 0.0
    self.bi = 0.00
      
  def autoPeak(self):# Finds most peaks
    global aPf, rusf, sx, comFlag, dt, id, Fflag, conTrunning
    mode = self.mode
    Fflag = True
    if conTrunning == True: return
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
    #print(izl,izr, ix0, 'curve')
    self.figure.clf()
    self.figure.text(0.04, 0.005,r'       $\phi$              Z             f', fontsize = 18)
    freqP = self.rus.reim.freq
    sstep = (freqP[2]-freqP[1])
    dataP = self.rus.reim.data
    if izl > 0:
         freqP = self.rus.reim.freq[izl:izr]
         dataP = self.rus.reim.data[izl:izr]
    data1 = np.real(dataP)
    data2 = np.imag(dataP)
    data3 = np.absolute(dataP)
    self.data3 = data3
    self.Lm = 0
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
        if tf == 'set': #this is where phase is user adjustable
            phD = self.phase.value()
            ph = math.radians(phD)
            co = math.cos(ph)
            si = math.sin(ph) 
    d1 = co*data1 - si*data2 #These are the curves to be plotted
    d2 = si*data1 + co*data2
    #find maximum
    i = np.argmax(data3)
    freqmax = float(freqP[i])
    #Lorentzian stuff *********************************************************************************************
    if self.rus.xFit == True:
      if izl > 1:
        zL = list() #inmitialize lists
        fL = list()
        mL = list()
        err = 0
        O0 = 90 #the guessed phase in degrees
        O1 = 90
        #now compute the fitted Lorentzian and the errot squared
        for self.niter in range(3):
            M0 = self.rus.h0.value() #the maximum data point for initial guess
            f0 = self.rus.f0.value() #the frequency at guessed maximum
            g0 = self.rus.w0.value() #the (guessed?) width
            f1 = self.rus.f1.value() #the frequency at guessed maximum
            g1 = self.rus.w1.value()
            M1 = self.rus.h1.value()
            O0 = math.radians(O0) # converted to radians
            O1 = math.radians(O1)
            zL = list()
            fL = list()
            mL = list()
            if M1 != 0:
                guess = [M0, f0, g0, O0, M1, f1, g1, O1, self.br, self.bi]
                try:
                    popt,pcov = curve_fit(self.cLrntz, freqP, data1, p0 = guess) #finds real background
                    self.rus.result.setText('success')
                    self.popt = popt
                    self.rus.h0.setValue(popt[0])
                    self.rus.f0.setValue(popt[1])
                    self.rus.w0.setValue(popt[2])
                    self.rus.h1.setValue(popt[4])
                    self.rus.f1.setValue(popt[5])
                    self.rus.w1.setValue(popt[6])
                    O0 = popt[3]
                    O1 = popt[7]
                    self.br = popt[8]
                    #self.bi = popt[9]
                    guess = popt
                    popt,pcov = curve_fit(self.dLrntz, freqP, data2, p0 = guess) #finds imag background
                    self.rus.result.setText('success')
                    self.bi = popt[9]
                    for j in range(izl, izr):
                      indx = j-izl
                      fj = freqP[indx]  #the measured frequency at current data point
                      out = self.cLrntz( fj, *popt) #compute Lorentzian
                      a0L = math.sqrt(self.Lm)
                      zL.insert(indx, out)
                      mL.insert(indx, a0L)
                      fL.insert(indx, fj)

                  #Lz = [f0, g0, M0]
                except:
                   self.rus.result.setText('retry')
            else:
                guess = [M0, f0, g0, O0, self.br ,self.bi]
                try:
                    popt,pcov = curve_fit(self.cLrntz, freqP, data1, p0 = guess) #finds real background
                    self.rus.result.setText('success')
                    #print(popt[1], popt[2], popt[0])
                    self.popt = popt
                    self.rus.h0.setValue(popt[0])
                    self.rus.f0.setValue(popt[1])
                    self.rus.w0.setValue(popt[2])
                    O0 = popt[3] 
                    self.br = popt[4]
                    #self.bi = popt[5]
                    guess = popt
                    popt,pcov = curve_fit(self.dLrntz, freqP, data2, p0 = guess) #finds imag background
                    self.rus.result.setText('success')
                    self.bi = popt[5]
                    #print(popt[1], popt[2], popt[0])
                    for j in range(izl, izr):
                      indx = j-izl
                      fj = freqP[indx]  #the measured frequency at current data point
                      out = self.cLrntz( fj, *popt) #compute Lorentzian
                      a0L = math.sqrt(self.Lm)
                      zL.insert(indx, out)
                      mL.insert(indx, a0L)
                      fL.insert(indx, fj)
                except:
                   self.rus.result.setText('retry')
    #add axes
    axes1 = self.figure.add_subplot(1,1,1)
    axes1.cla()
    axes1.xaxis.grid()
    axes1.yaxis.grid()
    if self.mode =='absm' or self.mode == 'reim':
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
            axes1.set_ylabel('Volts',fontsize=20)
            axes1.yaxis.label.set_color('black')
            #axes1.plot(freqP, d2, color = 'red') #adds other phase to plot
            axes1.plot(freqP, d1, color = 'blue')  ##
            if ((self.rus.xFit == True) and (self.niter ==2)):
                  axes1.plot(fL, zL, color = 'black', marker = 'o', markersize =1, linewidth = 0) # Lorentzian fit plot
                  self.phase.setValue(0.0)
        if self.mode == 'absm': 
            axes1.set_ylabel('Volts',fontsize=20)
            axes1.yaxis.label.set_color('black')
            axes1.plot(freqP, data3, color = 'red', label = 'magnitude')
            if ((self.rus.xFit == True) and (self.niter == 2)):
                  #print(self.br, self.bi)
                  axes1.plot(fL, mL, color = 'black', marker = 'o', markersize =1, linewidth = 0) # Lorentzian fit plot
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
    self.rus.xFit = False  #this blocks plots of data

  def cLrntz(self,fj, *gs): #computes in-phase Lorentzian
      if len(gs) == 6:
        M0 = gs[0]
        f0 = gs[1]
        g0 = gs[2]
        O0 = gs[3]
        br = gs[4]
        bi = gs[5]
        fsq0 = f0-fj
      if len(gs)== 10:
         M0 = gs[0]
         f0 = gs[1]
         g0 = gs[2]
         O0 = gs[3]
         M1 = gs[4]
         f1 = gs[5]
         g1 = gs[6]
         O1 = gs[7]
         br = gs[8]
         bi = gs[9]
         fsq0 = f0-fj
         fsq1 = f1-fj
      x0L = M0*g0*g0/(fsq0**2 + g0**2)
      y0L = M0*g0*fsq0/(fsq0**2 + g0**2)
      z0Lr = x0L*math.cos(O0) - math.sin(O0)*y0L
      z0Li = x0L*math.sin(O0) + math.cos(O0)*y0L
      if len(gs) == 10:
        x1L = M1*g1*g1/(fsq1**2 + g1**2)
        y1L = M1*g1*fsq1/(fsq1**2 + g1**2)
        z1Lr = x1L*math.cos(O1) - math.sin(O1)*y1L
        z1Li = x1L*math.sin(O1) + math.cos(O1)*y1L
      else:
         x1L = 0
         y1L = 0
         z1Lr = 0
         z1Li = 0
      zL0 = z0Lr + z1Lr +br
      zL1 = z0Li + z1Li +bi
      r0 = zL0
      r1 = zL1
      self.Lm = r1*r1 + r0*r0
      return zL0
  
  def dLrntz(self,fj, *gs):#computes quadrature Lorentzian
      if len(gs) == 6:
        M0 = gs[0]
        f0 = gs[1]
        g0 = gs[2]
        O0 = gs[3]
        br = gs[4]
        bi = gs[5]
        fsq0 = f0-fj
      if len(gs)== 10:
         M0 = gs[0]
         f0 = gs[1]
         g0 = gs[2]
         O0 = gs[3]
         M1 = gs[4]
         f1 = gs[5]
         g1 = gs[6]
         O1 = gs[7]
         br = gs[8]
         bi = gs[9]
         fsq0 = f0-fj
         fsq1 = f1-fj
      x0L = M0*g0*g0/(fsq0**2 + g0**2)
      y0L = M0*g0*fsq0/(fsq0**2 + g0**2)
      z0Lr = x0L*math.cos(O0) - math.sin(O0)*y0L
      z0Li = x0L*math.sin(O0) + math.cos(O0)*y0L
      if len(gs) == 10:
        x1L = M1*g1*g1/(fsq1**2 + g1**2)
        y1L = M1*g1*fsq1/(fsq1**2 + g1**2)
        z1Lr = x1L*math.cos(O1) - math.sin(O1)*y1L
        z1Li = x1L*math.sin(O1) + math.cos(O1)*y1L
      else:
         x1L = 0
         y1L = 0
         z1Lr = 0
         z1Li = 0
      zL0 = z0Lr + z1Lr +br
      zL1 = z0Li + z1Li +bi
      r0 = zL0
      r1 = zL1
      self.Lm = r1*r1 + r0*r0
      return zL1
  
  def save_Lrnttz(self):
    id0 = self.rus.Tkc.value()
    #self.Tlor.append(id0)
    popt = self.popt
    if popt == []: return
    l4 = popt[0:4]
    l4 = np.append(l4,id0)
    self.Lresult.append(l4)
    self.br = popt[4]
    self.bi = popt[5]
    if len(popt) == 10: 
       l8 = popt[4:8]
       l8 = np.append(l8,id0)
       self.Lresult.append(l8)
       self.br = popt[8]
       self.bi = popt[9]
    self.write_Lrntz_out(sorted(self.Lresult, key = lambda fq: fq[1]))
    self.rus.result.setText('')
  
  def write_Lrntz_out(self, lv):
    global dt
    # M0, f0, g0, O0, M1, f1, g1, O1
    print(self.rus.Tkc.value())
    if sx == '': 
       fd = os.path.join(path,'Lfit.dat')
    else: fd = os.path.join(path, + sx +'Lfit.dat')
    fd = os.path.join(path, sx + 'Lfit.dat')
    fh = open(fd, 'w' )
    for j in range (0, len(lv)):
            fh.write(str(f"{lv[j][4]:.4f}") + ' ' + str(f"{lv[j][1]:.4f}") + '  '+ str(f"{lv[j][2]:.3f}") + '  ' + str(f"{math.degrees(lv[j][3]):.3f}") + '  ' + str(f"{lv[j][0]:.2f}") + '\n')
    fh.close()

  def w0_guess(self):
     global ix0, ix1
     if ix0 > 0:
          freqP = self.rus.reim.freq[ix0:ix1]
          dataP = self.rus.reim.data[ix0:ix1]
          data3 = np.absolute(self.rus.reim.data[ix0:ix1])
          i = np.argmax(np.absolute(self.rus.reim.data[ix0:ix1])) + ix0
          j = np.argmax(data3)+ix0
          fz = float(self.rus.reim.freq[i])
          df = self.crsor.value()
          step = freqP[1]-freqP[0]
          if df > 0:
            if abs(fz-df) > step:
              self.rus.w0.setValue(abs(df-fz))
              self.rus.f0.setValue(self.rus.reim.freq[i])
              if (df-fz) != 0: self.rus.w1.setValue(abs(df-fz))
              self.rus.f0.setValue(freqP[i-ix0])

  def f1_guess(self):#fix this
     if self.crsor.value() >0:
        vc = self.crsor.value()
        pos = self.crsor.value()
        self.rus.f1.setValue(pos)
        g = self.rus.w0.value()
        self.rus.w1.setValue(g)
        j = self.find_dPoint(vc)[0]
        #hgt = self.temp.value()   #this is wrong 
        hgt = np.imag(self.rus.reim.data[j])
        self.rus.h1.setValue(hgt) 
  
  def mark(self): #marks peaks in window
    global rusf, conTrunning
    if conTrunning == True:return
    x = []
    x.append([])
    z = self.crsor.value()
    for y in rusf:
        if z == y[0]: return 
    x[-1].append(z)
    x[-1].append(self.temp.value())
    x[-1].append(1.0)
    rusf = rusf + x

  def unmark(self): #unmarks peaks in window
    global rusf, conTrunning
    if conTrunning == True:return
    if rusf == []:return
    self.wFreq()
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
    global rusf, c
    rusf = []
    #c = rusf
    self.wFreq()
    self.update(self.mode)

  def wFreq(self):#writes found peaks to file
    global rusf, c, b, sx, comFlag, dt, id, fn, Fflag, conTrunning
    if self.rus.auto == True:return
    if conTrunning == True:return
    #if rusf == []:return
    dat = dt
    id0 = self.rus.Tkc.value()
    tidF = str(int(id0+0.5))
    fn = os.path.join(path, tidF + '_' +sx + '_' + dat + '_F.dat')
    b = sorted(rusf)
    bx = [i[0] for i in b]
    by = [i[1] for i in b]
    b1 = [i[2] for i in b]
    if Fflag == True:
      if b != c:     
          with open(fn, 'w', newline='') as f:
              f.write('%8.2f' % id0)
              f.write(' ' +sx +'\n')
              for j in range(len(b)):
                  f.write('%3.7f  %5.4f  %5.1f \n' % (0.001*bx[j], by[j], b1[j]))
          f.close()  
    self.axes1.plot(bx,by,color = 'black', marker = matplotlib.markers.CARETDOWNBASE , markersize =10, linestyle = 'None')
    self.canvas.draw()
    c = b

  def update(self, mode):
    getattr(self, 'update_%s' % mode)()

  def reset(self):#resets plot
    global ix0, ix1
    ix0 = 0
    ix1 = len(self.rus.reim.data)
    izl = ix0
    izr = ix1
    tf = ''
    #print(izl,izr,ix0,'reset')
    h = np.max(abs(self.rus.reim.data[ix0:ix1]))
    xhgt = [-1.07*h,1.07*h]
    self.plot_curves(izl, izr, xhgt, tf)
    self.rus.set_Lenab(False)
    #self.wFreq()

  def onpress(self,event): #starts window selction
    global conTrunning, ix0, ix1
    #gets cursor coordinates
    if event.inaxes:
      fx = event.xdata
      if conTrunning == True:
        self.crsor.setValue(fx)
        return
      ix0 = self.find_dPoint(fx)[0]

  def onrelease(self,event):#ends window selection gets window coordinates into plot array
    global ix0, ix1
    if conTrunning == True:
        return
    if event.inaxes:
        fx = event.xdata
        ix1 = self.find_dPoint(fx)[0]
        izl = ix0
        izr = ix1
        if ix0 >= ix1: return
        fl = self.rus.reim.freq[ix0]
        fr = self.rus.reim.freq[ix1]
        lim = [fl,fr]
        self.axes1.set_xlim(lim)
        xh = np.max(abs(self.rus.reim.data[ix0:ix1]))
        xhgt = [-1.07*xh, 1.07*xh]
        tf = self.tf
        i = np.argmax(abs(self.rus.reim.data[ix0:ix1]))
        i = i +ix0
        #print(i)
        freqmax = self.rus.reim.freq[i]
        self.plot_curves(ix0, ix1, xhgt, tf)
        #self.rus.h0.setValue(xh)
        self.rus.h0.setValue(np.imag(self.rus.reim.data[i]))
        self.rus.f0.setValue(freqmax)
        self.wFreq()
        self.rus.set_Lenab(True)
  
  def plot_crsr(self):#plots vertical line at peak position and can be moved
    if self.mode =='tplot' : return
    if self.rus.EnableCursor == True:
        vc = self.crsor.value()
        axes2 = self.axes2
        hgt = self.axes1.get_ylim() #remember hgt before cla
        wdt = self.axes1.get_xlim()
        axes2.cla()
        axes2.tick_params(right = False)
        axes2.yaxis.set_ticklabels([])
        fre = self.find_dPoint(vc)[1] #findex of frequency of black cursor(find_dPoint returns j, fre)
        j = self.find_dPoint(vc)[0]
        posf = [fre,fre]
        h = np.abs(self.rus.reim.data[j])
        axes2.set_ylim(hgt) #reset y limit
        axes2.set_xlim(wdt)
        if self.mode == 'absm':val = [z, h]
        if self.mode == 'reim':val = [-h, h]
        self.temp.setValue(h) #puts the height of the cursor into spinbox
        #self.rus.h0.setValue(h)
        #self.rus.f0.setValue(fre)
        #print(self.temp.value(), self.crsor.value(), j)
        axes2.plot(posf, val, color = 'black', linewidth = 1, linestyle = 'solid')
        self.canvas.draw()
        self.wFreq()
        return
  
  def phase_set(self): #sets phase
    global ix0, ix1
    if self.mode == 'absm': return
    tf = 'set'
    self.tf = tf
    izl = ix0
    izr = ix1
    h = np.max(abs(self.rus.reim.data[izl:izr]))
    xhgt = [-1.07*h,1.07*h]
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()

  def update_reim(self): #this connects to re-im tab
    self.plot_reim()

  def plot_reim(self):# plots real and imginary components
    global ix0, ix1
    if len(self.rus.reim.data) < ix1: # this fixes if next data set is shorter than present one
       ix0 = 0
       ix1 = len(self.rus.reim.data)
    izl = ix0
    izr = ix1
    h = np.max(abs(self.rus.reim.data[izl:izr]))
    xhgt = [-1.07*h,1.07*h]
    self.mode = 'reim'
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()

  def update_absm(self):#this connects to abs tab
     self.plot_absm()     

  def plot_absm(self): #plots absolute value
    global ix0, ix1
    if len(self.rus.reim.data) < ix1: # this fixes if next data set is shorter than present one
       ix0 = 0
       ix1 = len(self.rus.reim.data)
    izl = ix0
    izr = ix1
    h = np.max(abs(self.rus.reim.data[izl:izr]))
    xhgt = [0,1.07*h]
    self.mode = 'absm'
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    self.wFreq()
    
  def update_tplot(self):# this connects to Tout tab
    if comFlag == False: self.plot_absm()
    else: self.plot_Tout()
     
  def plot_Tout(self): # plots temperature log
    self.mode = 'tplot'
    global tout, ttout, Fflag
    if(tout[0][0]) == '':return
    self.axes2.yaxis.set_ticklabels([])
    self.axes2.tick_params(right = False)
    tl = ix0
    tr = ix1
    th = np.max(tout)
    xhgt = [-1.07*th,1.07*th]
    #self.mode = 'tplot'
    tf = self.tf
    Fflag = False
    self.plot_curves(tl, tr, xhgt, tf)
    
  def find_dPoint(self, vc): #finds the closest data point to the crsr
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
    if j == 1: j = j + 1
    #print(j, fre)
    return j, fre

class rus(QMainWindow, Ui_rus):
  graphs = ['reim', 'absm']
  def __init__(self):
    global comFlag, sx
    super(rus, self).__init__()
    self.setupUi(self)
    global comFlag, comport, comport1
    settings = QSettings('rus.ini', QSettings.IniFormat)
    #self.rateValue.addItems([ '1500','500', '150', '50', '15', '5', '2'])
    self.rateValue.addItems([ '1500','500', '150', '50'])
    #self.rateValue.addItems(['5000', '1000', '500', '100', '50', '10', '5', '1'])
    #rate = [10, 50, 100, 500, 1000, 5000, 10000, 50000][value]
    self.rateValue.lineEdit().setReadOnly(True)
    self.rateValue.lineEdit().setAlignment(Qt.AlignRight)
    for i in range(self.rateValue.count()):
      self.rateValue.setItemData(i, Qt.AlignRight, Qt.TextAlignmentRole)
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
    self.f0.setValue(self.sweep_start)
    self.f1.setValue(self.sweep_start)
    #self.h0.setValue(0.5)
    self.h1.setValue(0)
    self.w0.setValue(0.01)
    self.mode = 'reim'  
    # create figures
    #Tab setup ###################################################################################################
    self.tabs = {}
    for i in range(len(self.graphs)):
      layout = getattr(self, '%sLayout' % self.graphs[i])
      self.mode = self.graphs[i]
      self.tabs[i] = FigureTab(layout, self)
    self.set_Lenab(False)
    self.stopSweep.setEnabled(False)
    # create TCP socket
    self.socket = QTcpSocket(self)
    self.socket.connected.connect(self.connected)
    self.socket.readyRead.connect(self.read_data)
    self.socket.error.connect(self.display_error)
    # connect signals from widgets
    self.connectButton.clicked.connect(self.start)
    #self.connectButton.clicked.connect(self.read_T)
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
    self.startTimer.timeout.connect(self.timeout) #freq sweep timer
    self.sweepTimer = QTimer(self)
    self.sweepTimer.timeout.connect(self.sweep_timeout)
    self.TstepTimer = QTimer()
    self.TstepTimer.timeout.connect((self.TconTimer)) #contols temperature control time between updates
    self.TsweepTimer = QTimer()
    self.TsweepTimer.timeout.connect(self.log_tout)
    self.EnableCursor = True #control writing data file
    self.stopChangeSetPoint = False # prevents set point change during open loop arduino sweep
    #operate temperature controller
    self.strT.setStyleSheet('QPushButton {background-color: yellow; color: black;}')
    self.strT.clicked.connect(self.log_tout) #starts log temperature to file
    self.stpT.clicked.connect(self.Tlog_Stop) #stops log temperature to file
    self.StopC.setEnabled(False)
    self.StopC.setStyleSheet('QPushButton {background-color: none; color: black;}')
    self.StopC.clicked.connect(self.TconStop) #stops temperature control
    self.Go.clicked.connect(self.FitResp)
    if sx == 'C' or sx == 'K': self.comID.setText(comdesig[0])
    self.setPt.valueChanged.connect(self.set_Temp)
    self.updtT.clicked.connect(self.read_T)
    self.sSw.clicked.connect(self.saveTstart)
    self.sSw.setStyleSheet('QPushButton {background-color: green; color: white;}')
    self.sSw.clicked.connect(self.TconTimer) #starts temperature controller
    self.sSw.clicked.connect(self.log_tout)
    self.Tstart.valueChanged.connect(self.final_T)
    self.Tstep.valueChanged.connect(self.final_T)
    self.nStep.valueChanged.connect(self.final_T)
    self.Stime.valueChanged.connect(self.final_T)
    #self.Tkc.setValue(id)
    self.Tstart0 = self.Tstart.value()
    self.setPt.setSuffix(sx)
    #self.setPt.setValue(int(id))
    self.saveT.setStyleSheet('QPushButton {background-color: green; color: white;}')
    self.saveT.clicked.connect(self.write_Tlog) #write temperature log to file
    self.delts = self.noise.value()
    self.dflag = True
    self.set_equil()
    self.read_T() #this reads and enters values in TC spinboxes
    if sx != '':
        tout[-1].append(float(self.Tkc.value())) #initializes temperature log list
        ttout[-1].append(ts1)
    self.n = 0
    self.Tsple = 0
    self.nrds = 0
    self.Tmean = 0.0
    self.ncycle = 0
    self.xcount = 0 #number of temperature  steps executed
    self.fstring = ''
    self.xflag = False # determines if magnified or full scan
    self.contRun = False # indicates when controller is running
    self.xlen = 0
    self.Tbar = np.zeros(16)
    self.Loop_on  = 0
    self.Loop.stateChanged.connect(self.LoopResp)
    self.xFit = False
    
    if sx == 'C' or sx == 'K': #enables or disables temperature controller inputs
       self.set_Tenab(True)
       self.set_Cenab(True)
       self.comID.setVisible(True)
       self.sSw.setEnabled(True)
       self.setPt.setEnabled(True)
    else:
       self.set_Tenab(False)
       self.set_Cenab(False)
       self.comID.setVisible(False)
       self.sSw.setEnabled(False)
       self.setPt.setEnabled(False)

    styles = {'color':'b', 'font-size':'20px'}
    self.Tlabel = 'Temperature(' + sx +')'
    self.graphWidget.setLabel('bottom', 'time(s)', **styles)
    self.graphWidget.setLabel('left', self.Tlabel, **styles)
    #self.graphWidget.setLabel('top', '', **styles)
    self.graphWidget.setLabel('top')
    self.graphWidget.setLabel('right')
    self.graphWidget.showGrid(x=True, y=True)
    if comFlag  == True: 
        self.Tplot(ttout[-1], tout[-1]) 

  def Tplot(self, time, temperature):
        #self.graphWidget.setSymbol('o')
        self.graphWidget.plot(time, temperature , border = 'True',  pen = None, symbol='o', symbolSize = '5')

  def FitResp(self):
       self.xFit = True
       self.update_tab()
       return

  def LoopResp(self):
     self.Loop_on  = self.Loop.checkState() 

  def saveTstart(self):
     self.Tstart0 = self.Tstart.value() 
 
  def set_equil(self): #sets equilibration test parameters and startup Ts to temperature on startup
    if self.dflag == False:return
    #self.read_T()
    if sx == 'C': 
        self.Tstart.setValue(round(float(self.Tkc.value()))) #sets starting T to what was measured on program start
        self.setPt.setValue(round(float(self.Tkc.value())))
        self.deltaS.setValue(120)
        self.deltaT.setValue(0.3)
        self.Tstep.setValue(5)
        self.dflag = False
    if sx == 'K': 
        self.Tstart.setValue(float(self.Tkc.value()))
        self.setPt.setValue(round(float(self.Tkc.value())))
        self.deltaS.setValue(360)
        self.deltaT.setValue(0.02)
        self.Tstep.setValue(1)
        self.dflag = False
    
  def log_tout(self): #continuous log of temperature
        global tout, ttout, comFlag, ts, sx
        self.strT.setStyleSheet('QPushButton {background-color: none; color: black;}')
        self.stpT.setStyleSheet('QPushButton {background-color: yellow; color: black;}')
        self.strT.setEnabled(False)
        self.stpT.setEnabled(True)
        self.TsweepTimer.start(5000)
        if comFlag == False: 
            return
        self.read_T()
        id = float(self.Tkc.value())
        tout[-1].append(id)
        ts1 = time.time()-ts #this is where log T array is generated to store
        ttout[-1].append(ts1)
        self.Tplot(ttout[-1], tout[-1]) 

  def slope(self): #sets slope for peak detection
      self.deltS = self.sweep_Hz/10
      self.noise.setValue(self.delts)
    
  def read_T(self): #primary temperature controller comm function
    global sx
    if comFlag == False: return
    
    if sx == 'C' or sx =='K':
      ser0.write(b't\n')     # specimen temperature for ACE units
      time.sleep(0.1)
      td = ser0.readline()
      td = td.decode('UTF-8')
      tds = td.strip('\n')
      if tds == '': return
      tds = float(tds)
      self.Tkc.setSuffix(sx)
      self.Tkc.setValue(tds)
      self.Tstart.setSuffix(sx)
      self.Tstep.setSuffix(sx)
      self.EndT.setSuffix(sx)

    if sx == 'k':
      ser0.write(b't\n')     # specimen temperature for user TC
      time.sleep(0.1)
      td = ser0.readline()
      td = td.decode('UTF-8')
      tds = td.strip('\n')
      if tds == '': return
      tds = float(tds)
      self.Tkc.setSuffix(sx)
      self.Tkc.setValue(tds)
      self.Tstart.setSuffix(sx)
      self.Tstep.setSuffix(sx)
      self.EndT.setSuffix(sx)

    if sx == 'C' or sx =='K':
      ser0.write(b's\n')     # get set point for ACE TC
      time.sleep(0.1)
      td = ser0.readline()
      td = td.decode('UTF-8')
      id = td.strip('\n')
      if id == '': return
      id = float(id)
      self.setPt.setValue(id)

    if sx == 'k':
      ser0.write(b's\n')     # get set point for user TC
      time.sleep(0.1)
      td = ser0.readline()
      td = td.decode('UTF-8')
      id = td.strip('\n')
      if id == '': return
      id = float(id)
      self.setPt.setValue(id)

    if sx == 'K': 
        ser0.write(b'c\n')     # control block temperature for ACE Thermoelectric stage
        time.sleep(0.1)
        td = ser0.readline()
        td = td.decode('UTF-8')
        id = td.strip('\n')
        id = float(id)
        self.bT.setSuffix(sx)
        self.bT.setValue(id)

    if sx == 'k': 
      ser0.write(b'c\n')     # control block temperature for user TC
      time.sleep(0.1)
      td = ser0.readline()
      td = td.decode('UTF-8')
      id = td.strip('\n')
      id = float(id)
      self.bT.setSuffix(sx)
      self.bT.setValue(id)

    
    if sx == 'C' or sx =='K':
        ser0.write(b'p\n')     # relative power
        time.sleep(0.1)
        td = ser0.readline()
        td = td.decode('UTF-8')
        id = td.strip('\n')
        id = float(id)
        if sx == 'C':id = id/10.23 #power correction factor for oven
        self.rP.setSuffix(' %')
        self.rP.setValue(id)

  def wSet_point(self,setT): #set point write to hardware
    global sx
    if (sx =='C' or sx =="K"): # for ACE TC
        ser0.write(setT)
    if sx =='k': # for user TC
        ser0.write(setT)

  def set_Temp(self):   
    if comFlag == False: return
    if self.stopChangeSetPoint == True: return #in arduino mode prevents changing set point during sweep
    if sx == 'K' or sx == 'C':  #for ACE TC
        Tset = self.setPt.value()
        if Tset < 10.0: Tset = 10.0
        Tset = str(Tset) + '\n'
        Tset = Tset.encode('UTF-8')
        self.wSet_point(Tset) # writes to TC
    if sx == 'k':      # for user TC
        Tset = self.setPt.value()
        if Tset < 10.0: Tset = 10.0
        Tset = str(Tset) + '\n'
        Tset = Tset.encode('UTF-8')
        self.wSet_point(Tset) # writes to TC

  def final_T(self):   #computes final temperature
      if comFlag == False: return
      Tstrt = self.Tstart.value()
      Tstp = self.Tstep.value()
      nStp = self.nStep.value()
      self.EndT.setValue(Tstrt + nStp*Tstp)

  def Tsweep(self):  #sets up arduino unmonitored temperature sweep
      self.strT.setEnabled(False)
      global conTrunning
      conTrunning = True
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
      self.stopChangeSetPoint = True #blocks set point change--code must be restarted to re enable this
      self.xflag = False
      scan = 'x,'+ sTstrt +','+ sTstp +','+ snStp +','+ ssTm
      tscan = scan.strip()+'\n'
      tss= tscan.encode('UTF-8')
      print(tss)
      ser0.write(tss) #arduino mode start T, T step, number of steps, time between steps
      conTrunning = False
      #self.TconStop()

  def TconTimer(self): # Timer for temperature stuff
    global conTrunning
    conTrunning = True
    self.sSw.setEnabled(False)
    self.StopC.setEnabled(True)
    self.deltaT.setEnabled(False)
    self.deltaS.setEnabled(False)
    self.Loop_on = self.Loop.checkState()
    self.TstepTimer.start(5000)
    Tstp = self.Tstep.value()
    nStp = int(self.nStep.value())
    sTm = self.Stime.value()
    self.EnableCursor = False #blocks cursor from plotting during T control
    self.contRun = True
    self.TconSet(Tstp, nStp,sTm) #starts the temperature control function
    if self.stopChangeSetPoint == False: self.read_T()
          
  def TconStop(self):#stops temperature control
    global conTrunning, Fflag
    self.TstepTimer.stop()
    self.StopC.setStyleSheet('QPushButton {background-color: none; color: black;}')
    self.sSw.setStyleSheet('QPushButton {background-color: green; color: white;}')
    self.sSw.setEnabled(True)
    self.StopC.setEnabled(False)
    self.CM.setEnabled(True)
    self.AP.setEnabled(True)
    self.set_Cenab(True)
    self.EnableCursor = True
    #self.stopChangeSetPoint = True
    self.contRun = False
    conTrunning = False
    Fflag = True
    self.n = 0

  def write_Tlog(self): #This writes the temperature log
      if comFlag == False: return
      ft = os.path.join(path, sx +'_' + dt + '_T.dat')
      fh = open(ft, 'w' )
      lT =len(tout[0])
      for k in  range(lT):
        fh.write('%12.2f %12.2f \n' %(ttout[0][k], tout[0][k]) )  
      fh.close()

  def Tlog_Stop(self): #stops continuous temperature log
    #global Fflag
    self.TsweepTimer.stop()
    self.strT.setEnabled(True)
    self.stpT.setEnabled(False)
    self.strT.setStyleSheet('QPushButton {background-color: yellow; color: black;}')
    self.stpT.setStyleSheet('QPushButton {background-color: none; color: black;}')      

  def TconSet(self, Tstp, nStp,sTm):#choses between periodic sweep and step and control and is main temp control function
    global sx, tout, Fflag
    self.EnableCursor = False
    self.StopC.setStyleSheet('QPushButton {background-color: red; color: black;}') #yellow is active
    self.sSw.setStyleSheet('QPushButton {background-color: none; color: black;}')
    if comFlag == False: return
    self.CM.setEnabled(False) #disable a bunch of controls
    self.AP.setEnabled(False)
    self.set_Cenab(False) 
    Fflag = False
    self.delS = self.deltaS.value()
    self.delT = self.deltaT.value()
    self.ncycle = self.ncycle + 1
    nStp0 = nStp
    if sx =='C':
        navg = 16 #16 = 3.5 minutes
    if sx =='K':
        navg = 16
    if sTm > 0: #use arduino mode
        self.TconStop()
        self.stopChangeSetPoint = True   
        self.Tsweep()
        return
    else: #step and wait mode
            self.stopChangeSetPoint = False
            #Tm = 0
            nlist = len(tout[0]) #length of the temperature log file
            Teq = tout[0][-1] #current specimen temperature
            nrds0 = int(self.delS/5) #of reads for averaging because 5 seconds between reads
            n = self.n #temperature step counter
            if n > nStp0: #end of process of  (first) sweep
                if self.Loop_on == 2:
                   self.Loop.setCheckState(False)
                   #self.Loop_on = 0
                   Tstrt = self.setPt.value()-Tstp #setpoint to last temperature of loop
                   self.n = 0
                   Tstp = -Tstp #Change sign of steps
                   nStp0 = nStp-1
                   self.nStep.setValue(nStp0)# is this the issue?
                   self.Tstep.setValue(Tstp)
                   self.Tstart.setValue(Tstrt) #new starting point 
                   self.EndT.setValue(Tstrt + nStp0*Tstp) #sets end point to original start point
                   #self.Tstart0 = self.EndT.value() #reset end temperature later
                   self.nrds = 0
                   self.Tmean = 0.0
                   self.ncycle = 0
                else:
                  self.setPt.setValue(self.Tstart0) 
                  self.Tbar = np.zeros(navg)
                  self.write_Tlog()# saves temperature log
                  self.n = 0
                  self.nrds = 0
                  self.Tmean = 0.0
                  self.ncycle = 0
                  self.cancel()
                  self.TconStop()
                  self.EnableCursor = True
                  return   
            #print(self.Loop_on)        
            Tstrt = self.Tstart.value() + n*Tstp #this is current set point ************************************************
            self.setPt.setValue(Tstrt) #This actually implements new setpoint *****************************************************
            if(nlist > navg): self.Tbar = tout[0][nlist-navg:nlist] #take last temperaure points to average
            if ((nlist >= nrds0) and (self.ncycle  > navg)): #nrds = 0 when T near set point
                self.Tmean = np.mean(self.Tbar) #get average temperature for most recent navg
            #self.delT is the delta T that defines near equilibrium
            if sx == 'C':
                 if abs(self.Tmean-Tstrt) < self.delT: #increments the counter for determining stability
                    self.Tsple = self.Tmean
                    self.nrds = self.nrds + 1 
                 #else: self.nrds = 0   #resets counter if any T is off

            if sx == 'K': #checks for difference between T nrds0 ago and now
                if abs(self.Tmean-Teq) <= self.delT:
                    self.Tsple = self.Tmean
                    self.nrds = self.nrds + 1 #increments counter when T close to equilibrium only
                #else:self.nrds = 0   #resets counter if any T is off

            if(self.nrds > nrds0): #this is where the number of good components of the average is added
                self.TstepTimer.stop() #stops temperature stuff to give time to write last freq file
                if self.Xstep.value() == 0:  self.trigResc()# saves full sweep on reaching equilibrium
                if self.Xstep.value() > 0:  # saves xscans on reaching equilibrium
                    self.read_ini()
                    self.sweep_auto()
                n = n+1
                self.n = n #step counter is incremented here
                self.Tbar = np.zeros(navg) #resets averages
                self.nrds = 0
                self.Tmean = 0.0
                self.Tsample = 0
                self.ncycle = 0
     
  def trigResc(self): # starts single sweep
    mode = self.mode
    # here will switch between fullscan and xscan
    self.sweep(mode) #this is where full scan is written to file
    self.auto = False

  def set_enabled(self, enabled):
    widgets = [self.rateValue, self.level1Value, self.startValue, self.stopValue, self.stepValue, \
     self.singleSweep, self.autoSweep, self.swpD, self.comID, self.freqX, self.Xstep]
    for entry in widgets:
      entry.setEnabled(enabled)

  def set_Tenab(self, enabled): # enables controls if T controller is connected
    Twidg = [self.updtT, self.Tkc, self.bT, self.rP, self.EndT, self.comID, self.saveT,self.strT, self.stpT]
    for entry in Twidg:
      entry.setEnabled(enabled)

  def set_Cenab(self, enabled): #enables or disables many controls while T controller is running
    Cwidg = [self.deltaT, self.deltaS, self.Tstart, self.Tstep, self.nStep, self.Stime]
    for entry in Cwidg:
      entry.setEnabled(enabled)

  def set_Lenab(self,enabled): #enable fitting inputs
     Lwidg = [self.f0, self.w0, self.h0, self.f1, self.w1, self.h1, self.w0g, self.f1g , self.result, self.save, self.Go]
     for entry in Lwidg:
        entry.setEnabled(enabled)

  def set_Venab(self, enabled):
    Vwidg = [self.DVM]
    for entry in Vwidg:
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
      self.progressBar.setValue(int((self.offset + size) / 16))
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
    self.sweep_stop = value
    self.nSize.setValue(self.sweep_size)

  def set_size(self):
    global ix0,ix1
    Hz = self.stepValue.value()
    #ix0 = int(1000*self.sweep_start/Hz)
    ix0 = 0
    #ix1 = int(1000*self.sweep_stop/Hz)
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/Hz)
    if self.sweep_size > 32766:
      self.sweep_size = 32766
      Hz = int((1000*(self.sweep_stop-self.sweep_start))/32766 + 0.5)
      if Hz <=1: Hz = 1
    self.stepValue.setValue(Hz)
    ix1 = self.sweep_size
    self.nSize.setValue(self.sweep_size)
    self.noise.setValue(0.2*Hz +1) #guess for slope parameter

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
    #global ix0, ix1
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
    if self.auto == True and self.xflag == False: self.write_dat() # writes data file at each temperature for full scan

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
    global ix0, ix1
    if not self.reading:
        if self.xflag == True and self.stopChangeSetPoint == False:  #captures regions around detected peaks
            if self.xcount >= self.xlen-1: #ends x scans
                self.xflag = False
                self.sweepTimer.stop
                f = self.reim.freq
                d = self.reim.data
                for i in range(len(f)):
                    self.fd.write('%12.5f  %12.10f  %12.10f\n' % (f[i], d.real[i], d.imag[i]))#saves last scan
                self.fd.write('%12.5f  %12.10f  %12.10f\n' % (self.inir, 0.0, 0.0))#adds point at end of original sweep
                self.startValue.setValue(self.inil)
                self.stopValue.setValue(self.inir)
                self.rateValue.setCurrentIndex(self.inirate)
                self.stepValue.setValue(self.iniHz)
                self.fxn.close()
                self.fd.close()
                self.cancel()
                self.xcount = 0
                self.update_tab()
                if self.setPt.value() == self.EndT.value():
                  self.setPt.setValue(self.Tstart.value())
                  self.write_Tlog()
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
                self.set_rate(indx) #sets to slowest scan rate for best line shape
                self.rateValue.setCurrentIndex(indx)
                Hz = self.Xstep.value()
                self.stepValue.setValue(Hz)
                self.nSize.setValue(1000) #size  = 1000
                self.sweep_Hz = Hz
                strt = 1000*fx-0.5*Hz #this sets peak at center of window
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
    self.mode = self.graphs[index]
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
    settings.setValue('step',self.stepValue.value())
  
  def read_cfg_settings(self, settings):
    self.addrValue.setText(settings.value('addr', '192.168.1.100'))
    self.rateValue.setCurrentIndex(settings.value('rate', 0, type = int))
    self.level1Value.setValue(settings.value('level_1', 1, type = float))
    reim_start = settings.value('reim_start',200, type = int)
    reim_stop = settings.value('reim_stop', 450, type = int)
    self.sweep_Hz = settings.value('step',50,type = int)
    self.startValue.setValue(reim_start)
    self.stopValue.setValue(reim_stop)
    self.stepValue.setValue(self.sweep_Hz)

  def read_dat(self): #read data from file
    global ix0, ix1
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
      self.level1Value.setValue(vfd)
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
      Fstep = int(self.reim.freq[3] - self.reim.freq[2])
      Fstep = 1000*Fstep
      self.noise.setValue(0.2*Fstep +1) #guess for slope parameter
      fh.close 

  def write_dat(self): #write data to file
    global sx, comFlag,dt,ts
    trds = '0.0'
    tdat= time.time()-ts
    t = datetime.now()
    dat = t.strftime('%d%m%Y%H%M%S') #get date time string to identify files
    if comFlag == True: #Temperature controller connected and auto
      if self.auto == True:
        self.read_T()
        trds = str(self.Tkc.value())
        fd = os.path.join(path, trds + '_' +sx +'_' + dat + '_D.dat')
      else: #Temperature controller connected but not auto uses set point to label file but puts actual temperature in file
        #dat = dt
        tid = str(self.setPt.value())
        fd = os.path.join(path, tid + '_' +sx +'_' + dat + '_D.dat')
    else: #no temperature controller connected just uses 0.0 for temperature
        fd = os.path.join(path, trds + '_' +sx +'_' + dat + '_D.dat')
    f = self.reim.freq
    d = self.reim.data
    self.read_T()
    Tsample = float(self.Tkc.value())
    if self.Tsple != 0:
            Tsample = self.Tsple
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
os.environ['QT_AUTO_SCREEN_SCALE_FACTOR'] = '1'
warnings.filterwarnings('ignore')
app = QApplication(sys.argv)
splash_pix = QPixmap('Splash.png')
splash = QSplashScreen(splash_pix)
splash.show()
#time.sleep(1)
app.processEvents()
path = os.getcwd() #get current working directory
path = os.path.join(path, 'RUS_Data') #append RUS_Data to path
if os.path.exists(path) == False: os.makedirs(path)
sqr2 = 1/math.sqrt(2)
z = 0.0
ts = time.time()
t = datetime.now()
Fflag = True #Flag to block writes of T log file while T controller is running
conTrunning = False #Flag to indicate if T controller is running
dt = t.strftime('%d%m%Y%H%M%S') #time that program starts
id = 0.0
comFlag = False #below detects if ACE controller connected
comport = ''
comport1 = ''
sx = ''
dx = ''
ex = ''
sx1 = ''
aPf = True
comdesig = []
c = []
ix0 = 0
ix1 = 2
ports = serial.tools.list_ports.comports()
nports = len(ports)
if nports > 0:
    for iports in range(0, nports):
        pq = ports[iports]
        comport = pq.device #this gets each comport designator
        print(str(pq)[7:])
        #if (str(pq)[7:23] == 'USB-SERIAL CH340') or (str(pq)[8:24] == 'USB-SERIAL CH340') or (str(pq)[8:15] == 'Arduino') or  (str(pq)[7:14] == 'Arduino'): 
        if (str(pq)[7:12] != 'Intel') and (str(pq)[8:13] != 'Intel') :
            unkP = serial.Serial(comport, 9600, bytesize=8, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, timeout=2, rtscts=0)  # open serial port 
            l2 = unkP.readline()
            time.sleep(2) #this delay is critical to enable ports to open, 2 or 3 s required
            unkP.write(b'*IDN?\n')     # get id for controllers that sends an ID
            time.sleep(2.0)
            l2 = unkP.readline()
            lid = l2.decode('UTF-8')
            ly = lid
            lidy = str(ly.strip())
            #print(comdesig)

            if (lidy == 'ACETC'):
              comFlag = True
              ser0 = unkP
              comdesig.append(str(pq)[7:])
              ser0.write(b'm\n')     # get TC controller type (K, C)
              md = ser0.readline()
              md0 = md.decode('UTF-8')
              dx = str(md0.strip())
              if dx == 'K' or dx == 'C': sx = dx
              break
  
    ts = time.time()
    tout = []
    ttout = []
    tout.append([]) #setup of temperature log list
    ttout.append([])
    ts1 = time.time()-ts

window = rus()
window.update_tab()
window.show()
splash.finish(window)
sys.exit(app.exec_())
