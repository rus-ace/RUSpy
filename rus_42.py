#!/usr/bin/env python3

# Control program for the Red Pitaya vector network analyzer
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

from operator import pos
import sys
import struct
import warnings

from functools import partial
from matplotlib.backend_bases import Event

import numpy as np
import math

import matplotlib
from numpy.core.arrayprint import str_format
from numpy.lib.type_check import imag, real
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import Formatter, FuncFormatter
from matplotlib.widgets import Cursor
#import matplotlib.pyplot as plt

from PyQt5.uic import loadUiType
from PyQt5.QtCore import QRegExp, QTimer, QSettings, QDir, Qt
from PyQt5.QtGui import QRegExpValidator, QPalette, QColor
from PyQt5.QtWidgets import QApplication, QDoubleSpinBox, QMainWindow, QMessageBox, QDialog, QFileDialog, QPushButton, QLabel, QSpinBox
from PyQt5.QtNetwork import QAbstractSocket, QTcpSocket

sqr2 = 1/math.sqrt(2)
z = 0.0


Ui_rus, QMainWindow = loadUiType('rus.ui')

class Measurement:
  def __init__(self, start, stop, size):
    self.freq = np.linspace(start, stop, size)
    self.freq2 = self.freq
    self.data = np.zeros(size, np.complex64)
    self.period = 62500
    

class NavigationToolbar(NavigationToolbar):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ( 'Home', 'Pan', 'Zoom', 'Save')]   

class FigureTab:
  #event = 'motion_notify_event'
  #event = 'xlim_changed'
  def __init__(self, layout, rus):
    # create figure
    self.figure = Figure()
    if sys.platform != 'win32':
      self.figure.set_facecolor('none')
    self.canvas = FigureCanvas(self.figure)
    self.canvas.mpl_connect('button_press_event', self.onclick)
    layout.addWidget(self.canvas)
    self.toolbar = NavigationToolbar(self.canvas, None, False)
    self.toolbar.layout().setSpacing(6)
    layout.addWidget(self.toolbar)
    tf = ''
    self.tf = tf
    ##################################################################################
    #add crsor
    self.crsor = QDoubleSpinBox()
    self.crsor.setDecimals(3)
    self.crsor.valueChanged.connect(self.plot_crsr)
    self.toolbar.addWidget(self.crsor)
    #add value at crsor
    self.temp = QDoubleSpinBox()
    self.temp.setReadOnly(True)
    self.temp.setButtonSymbols(2)
    self.temp.setDecimals(6)
    self.toolbar.addWidget(self.temp)
    #add set phase
    self.mode = 'reim'
    self.phase = QSpinBox()
    self.phase.setSingleStep(5)
    self.phase.setMaximum(180)
    self.phase.setMinimum(-180)
    if self.mode == 'reim': self.phase.valueChanged.connect(self.phase_set)
    self.toolbar.addWidget(self.phase)
    #add  best phase button
    self.bPhase = QPushButton('Best Phase')
    self.bPhase.setStyleSheet('background-color:orange')
    self.bPhase.clicked.connect(self.phase_adj)
    self.toolbar.addWidget(self.bPhase)
    #read best phase
    self.best = QSpinBox()
    self.best.setReadOnly(True)
    self.best.setButtonSymbols(2)
    self.best.setMaximum(180)
    self.best.setMinimum(-180)
    self.toolbar.addWidget(self.best)
    #add rescale button
    self.plotButton = QPushButton('Rescale')
    self.plotButton.setStyleSheet('background-color:blue; color:rgb(255, 255, 255)')
    self.plotButton.clicked.connect(self.plot_reim)
    self.toolbar.addWidget(self.plotButton)
    #add mark peak button
    self.peak = QPushButton('Mark')
    self.peak.clicked.connect(self.mark)
    self.toolbar.addWidget(self.peak)
    layout.addWidget(self.toolbar)
    #self.figure.subplots_adjust(left = 0.09, bottom = 0.1, right = 0.91, top = 0.96)
    axes1 = self.figure.add_subplot(1,1,1)
    axes2 = self.figure.add_subplot(1,1,1)
    axes3 = self.figure.add_subplot(1,1,1)
    self.axes1 = axes1
    self.axes2 = axes2
    self.axes3 = axes3
    self.mode = 'reim'
    self.rus = rus

####################################################################################
  
  def plot_curves(self, izl, izr, xhgt, tf):
    matplotlib.rcdefaults()
    matplotlib.rcParams['axes.formatter.use_mathtext'] = True
    self.figure.clf()
    self.figure.text(0.12, 0.005,r'  f          z       $\phi$       $\phi$$_0$', fontsize = 18)
    freq = self.rus.reim.freq
    data = self.rus.reim.data
    print(izl,izr)
    if izl > 0:
         freq = self.rus.reim.freq[izl:izr]
         data = self.rus.reim.data[izl:izr]
    data1 = np.real(data)
    data2 = np.imag(data)
    data3 = np.absolute(data)
    self.data3 = data3
    y = data3.max()
    h = [-y,y]
    lim = y
    if izl > 0: h = xhgt
    lim = xhgt
    ph = 0
    phD = 0
    co = 1
    si = 0
    #operate on phase
    if tf!='':
      if self.mode == 'reim':
        if tf == 'set':
            phD = self.phase.value()
            ph = math.radians(phD)
            tf = ''
        if tf == 'adj':
            d1 = sum(data1)
            d2 = sum(data2)
            if d2 == z: d1 = 1.0
            ph = (math.atan2(-d2, d1))
            phD = math.degrees(ph)
            self.best.setValue(phD)
            tf = ''
    co = math.cos(ph)
    si = math.sin(ph)
    d1 = co*data1 - si*data2
    d2 = si*data1 + co*data2
    self.phase.setValue(phD)
    #find maximum
    i = np.argmax(data3)
    freqmax = float(freq[i])
    vmax = data3[i]
    freq3 = [freqmax,freqmax]
    reim_Max = [-vmax,vmax]
    #add axes
    axes1 = self.figure.add_subplot(1,1,1)
    axes1.cla()
    axes1.xaxis.grid()
    axes1.set_xlabel('kHz')
    axes1.set_ylabel('real')
    #Setup crsor spinbox
    self.crsor.setMaximum(freq[-1])
    self.crsor.setMinimum(freq[0])
    sstep = (freq[1]-freq[0])
    self.crsor.setSingleStep(sstep)
    self.crsor.setValue(freqmax)
    #setup axes
    hgt = data3[i]
    xlim = [freq[0],freq[-1]]
    axes1.set_xlim(xlim)
    axes1.tick_params('y', color = 'black', labelcolor = 'black')
    axes1.yaxis.label.set_color('blue')
    axes2 = axes1.twinx()
    axes2.set_ylabel('imag')
    axes2.set_xlim(xlim)
    if self.mode == 'reim': 
      lim = h        
    if self.mode == 'absm': 
      lim = [z, h[1]]
    axes3 = axes1.twinx()
    self.axes3 = axes3
    if lim is not None: 
        axes1.set_ylim(lim)
        axes2.set_ylim(lim)
        axes3.set_ylim(lim)
    # initiate plots
    if self.mode == 'reim': 
        self.temp.setValue(np.abs(data3[i]))
        axes2.yaxis.label.set_color('red')
        self.curve1, = axes1.plot(freq, d1, color = 'blue', label = 'real')
        self.curve2, = axes1.plot(freq, d2, color = 'red', label = 'imag')  ##
    if self.mode == 'absm': 
        self.temp.setValue(np.abs(data1[i]))
        #self.temp.setValue(np.abs(sqr2*data1[i]))
        axes1.set_ylabel('magnitude')
        axes1.yaxis.label.set_color('red')
        axes2.tick_params('y',color = 'black', labelcolor = 'black')
        axes2.set_ylabel('')
        #hgt = sqr2*hgt
        self.curve2, = axes1.plot(freq, data3, color = 'red', label = 'magnitude')
    #plot line at max
    self.curve3, = axes2.plot(freq3, reim_Max, color = 'blue', linewidth = 1, linestyle = '--')
    #enable cursor
    self.cur = Cursor(axes3, horizOn=True, vertOn=True, useblit=True, color = 'blue', linewidth = 1)
    self.canvas.draw()
    #enable zoom box
    axes3.callbacks.connect('xlim_changed',self.on_xlims_change)
  
  def mark(self):
    pass
  def update(self, mode):
    getattr(self, 'update_%s' % mode)()

  def onclick(self,event):
    #gets cursor coordinates
    ix, iy = event.xdata, event.ydata
    coords = [ix, iy]
    self.coords = coords
 
  def on_xlims_change(self, axes3):
    #sets value of max into spin boxes
    izl = self.zbox()[0]
    izr = self.zbox()[1]
    km = self.zbox()[2]
    xhgt = self.zbox()[3]
    d = np.abs(self.rus.reim.data)
    f = self.rus.reim.freq
    self.temp.setValue(d[km])
    self.crsor.setValue(f[km])
    return izl, izr, km, xhgt

  def plot_crsr(self):
    #plots vertical line at peak position and can be moved
    axes3 = self.axes3
    hgt = self.axes3.get_ylim() #remember hgt before cla
    wdt = self.axes3.get_xlim()
    axes3.cla()
    # axes3.tick_params(right = False)
    # axes3.yaxis.set_ticklabels([])
    fre = self.find_dPoint()[1]
    j = self.find_dPoint()[0]
    posf = [fre,fre]
    h = np.abs(self.rus.reim.data[j])
    axes3.set_ylim(hgt) #reset y limit
    axes3.set_xlim(wdt)
    offset = (wdt[1]-wdt[0])/20
    if self.mode == 'absm':val = [z, h]
    if self.mode == 'reim':val = [-h, h]
    self.temp.setValue(h)
    axes3.plot(posf, val, color = 'black', linewidth = 1, linestyle = 'solid')
    axes3.text(fre+offset, h, round(fre,3), bbox=dict(facecolor = 'yellow', alpha = 0.5))
    self.canvas.draw()
    return

  def phase_adj(self, tf):
    #trigger to roll the phase
    if self.mode =='absm': return
    tf = 'adj'
    self.tf = tf
    getattr(self, 'plot_%s' % self.mode)()
  

  def phase_set(self, tf):
    #trigger to set best phase
      if self.mode == 'absm': return
      tf = 'set'
      self.tf = tf
      getattr(self, 'plot_%s' % self.mode)()

  def plot_reim(self):
    izl = self.zbox()[0]
    izr = self.zbox()[1]
    km = self.zbox()[2]
    xhgt = self.zbox()[3]
    c = self.axes3.get_xlim()
    tf = self.tf
    #print(izl, izr, km, fl)
    self.plot_curves(izl, izr, xhgt, tf)


  def update_reim(self):
      self.plot_reim()
      #getattr(self, 'plot_%s' % self.mode)()

  def plot_absm(self):
    izl = self.zbox()[0]
    izr = self.zbox()[1]
    km = self.zbox()[2]
    xhgt = self.zbox()[3]
    tf = self.tf
    self.mode = 'absm'
    tf = self.tf
    self.plot_curves(izl, izr, xhgt, tf)
    ##############################################################################################################

  def update_absm(self):
     self.plot_absm()

  def zbox(self):
    sz = self.rus.sweep_size
    span = self.axes3.get_xlim()
    xhgt = self.axes3.get_ylim()
    izl = 0
    izr = sz
    km = 1
    if span[0] == 0:
      #print('zero')
      return izl, izr, km, span
    f = self.rus.reim.freq
    d = np.absolute(self.rus.reim.data)
    f0 = span[0]
    for i in range(sz):
        if f[i] >= f0:
            izl = i
            break
    f1 = span[1]
    for j in range(sz):
        if f[j] >= f1:
            izr = j
            break
    km = np.argmax(d[izl:izr])
    km = km +izl-1
    izr = izr + 1
    #print(izl,izr)
    return izl, izr, km, xhgt

  def find_dPoint(self):
      #finds the closest data point to the crsr
    fv = self.crsor.value()
    f = self.rus.reim.freq
    i = 0
    for fre in f:
        i = i+1
        if fre >= fv: 
            break
    fre = f[i]
    j = i
    return j, fre

class rus(QMainWindow, Ui_rus):
  graphs = ['reim', 'absm']
  def __init__(self):
    super(rus, self).__init__()
    self.setupUi(self)
    # address validator
    rx = QRegExp('^(([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\.){3}([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])|rp-[0-9A-Fa-f]{6}\.local$')
    self.addrValue.setValidator(QRegExpValidator(rx, self.addrValue))
    # state variables
    self.idle = True
    self.reading = False
    self.auto = False
    # sweep parameters
    self.sweep_start = 300
    self.sweep_stop = 400
    self.sweep_Hz = 40
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/self.sweep_Hz)
    #sstep = self.sweep_Hz
    if(self.sweep_size)>32766:
      self.sweep_size = 32766
      self.sweep_Hz = (self.sweep_stop-self.sweep_start)/32766
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
    #########################################################################################################
    self.singleSweep.clicked.connect(self.sweep)
    self.autoSweep.clicked.connect(self.sweep_auto)
    self.stopSweep.clicked.connect(self.cancel)
    self.datButton.clicked.connect(self.write_dat)
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

  def set_enabled(self, enabled):
    widgets = [self.rateValue, self.level1Value, self.startValue, self.stopValue, self.stepValue, self.singleSweep, self.autoSweep]
    
    for entry in widgets:
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
        self.update_tab()
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

  def set_size(self, value):
    Hz = self.stepValue.value()
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/Hz)
    if self.sweep_size > 32766:
      self.sweep_size = 32766
      Hz = int(1000*(self.sweep_stop-self.sweep_start))/32766
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

  def cancel(self):
    self.sweepTimer.stop()
    self.auto = False
    self.reading = False
    self.socket.write(struct.pack('<I', 11<<28))
    self.progressBar.setValue(0)
    self.set_enabled(True)

  def sweep_auto(self):
    self.auto = True
    self.sweepTimer.start(1000)

  def sweep_timeout(self):
    if not self.reading:
      self.sweep('reim')

  def update_tab(self):
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
    self.rateValue.setCurrentIndex(settings.value('rate', 0, type = int))
    self.level1Value.setValue(settings.value('level_1', 1, type = int))
    reim_start = settings.value('reim_start', 300, type = int)
    reim_stop = settings.value('reim_stop', 400, type = int)
    #reim_size = settings.value('reim_size', 1, type = int)
    self.sweep_Hz = settings.value('step',40,type = int)
    self.startValue.setValue(reim_start)
    self.stopValue.setValue(reim_stop)
    #self.stepValue.setValue(reim_size)
    self.stepValue.setValue(self.sweep_Hz)
    
  def write_dat(self):
    dialog = QFileDialog(self, 'Write dat file', '.', '*.dat')
    dialog.setDefaultSuffix('dat')
    dialog.setAcceptMode(QFileDialog.AcceptSave)
    dialog.setOptions(QFileDialog.DontConfirmOverwrite)
    if dialog.exec() == QDialog.Accepted:   
      name = dialog.selectedFiles()
      fh = open(name[0], 'w')
      f = self.reim.freq
      d = self.reim.data
      fh.write('frequency    real         imag\n')
      for i in range(f.size):
        fh.write('%12.2f  %12.7f  %12.7f\n' % (f[i] * 1000, d.real[i], d.imag[i]))
      fh.close()

warnings.filterwarnings('ignore')
app = QApplication(sys.argv)
window = rus()
window.update_tab()
window.show()
sys.exit(app.exec_())
