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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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


Ui_rus, QMainWindow = loadUiType('rus.ui')

class Measurement:
  def __init__(self, start, stop, size):
    self.freq = np.linspace(start, stop, size)
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
    # create navigation toolbar
    self.toolbar = NavigationToolbar(self.canvas, None, False)
    self.toolbar.layout().setSpacing(6)
    # remove subplots action
    #actions = self.toolbar.actions()
    layout.addWidget(self.toolbar)
    tf = 'adj'
    self.tf = tf
    ##################################################################################
    #add crsor
    self.crsor = QDoubleSpinBox()
    self.crsor.setDecimals(3)
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
    self.plotButton.clicked.connect(self.plot)
    self.toolbar.addWidget(self.plotButton)
    #add mark peak button
    self.peak = QPushButton('Mark')
    self.peak.clicked.connect(self.mark)
    self.toolbar.addWidget(self.peak)
##################################################################################
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

  def xlim(self, freq):
    start = freq[0]
    stop = freq[-1]
    min = np.minimum(start, stop)
    max = np.maximum(start, stop)
    margin = (max - min) / 50
    return (min - margin, max + margin)

  def plot(self):
    #self.crsor.valueChanged.connect(self.plot_crsor)
    getattr(self, 'plot_%s' % self.mode)()


  def update(self, mode):
    #self.crsor.valueChanged.connect(self.plot_crsor)
    getattr(self, 'update_%s' % mode)()

####################################################################################
####################################################################################
  
  def plot_curves(self, freq, data1, label1, limit1, data2, label2, limit2, tf):
    matplotlib.rcdefaults()
    matplotlib.rcParams['axes.formatter.use_mathtext'] = True
    self.figure.clf()
    self.figure.text(0.12, 0.005,r'  f                   z               $\phi$                            best $\phi$')
    data4 = np.abs(data1 +1j*data2)
    self.data4 = data4
    #limit2 = limit1
    ph = 0
    phD = 0
    co = 1
    si = 0
    #self.tf = ''
    if tf!='':
      if self.mode == 'reim':
          if tf == 'set':
                  phD = self.phase.value()
                  ph = math.radians(phD)
                  tf = ''
          if tf == 'adj':
                  d1 = sum(data1)
                  d2 = sum(data2)
                  if d2 == 0.0: d1 = 1.0
                  ph = (math.atan2(-d2, d1))
                  phD = math.degrees(ph)
                  self.best.setValue(phD)
                  tf = ''
          co = math.cos(ph)
          si = math.sin(ph)
    self.phase.setValue(phD)
    d1 = co*data1 - si*data2
    d2 = si*data1 + co*data2
    self.freq = freq
    self.limit2 = limit2
    #find maximum
    i = np.argmax(data4)
    freqmax = float(freq[i])
    vmax = data4[i]
    freq3 = [freqmax,freqmax]
    if self.mode == 'absm': data3 = [0.0, sqr2*vmax]
    data3 = [-vmax,vmax]
    #add axes
    #self.figure.subplots_adjust(left = 0.09, bottom = 0.1, right = 0.91, top = 0.96)
    axes1 = self.figure.add_subplot(1,1,1)
    axes1.cla()
    axes1.xaxis.grid()
    axes1.set_xlabel('kHz')
    axes1.set_ylabel(label1)
    xlim = self.xlim(freq)
    axes1.set_xlim(xlim)
    ###################################################################################
    #Setup crsor spinbox
    self.crsor.setMaximum(freq[-1])
    self.crsor.setMinimum(freq[0])
    sstep = (freq[1]-freq[0])
    self.crsor.setSingleStep(sstep)
    self.crsor.setValue(freqmax)
    hgt = data4[i]
    if self.mode == 'reim': self.temp.setValue(np.abs(sqr2*data4[i]))
    if self.mode == 'absm':  
        self.temp.setValue(np.abs(sqr2*data1[i]))
        hgt = sqr2*hgt
    if limit1 is not None: axes1.set_ylim(limit1)
    ##############################################################################
    self.curve1, = axes1.plot(freq, d1, color = 'blue', label = label1)
    if data2 is None:
      self.canvas.draw()
      return
    axes1.tick_params('y', color = 'blue', labelcolor = 'blue')
    axes1.yaxis.label.set_color('blue')
    axes2 = axes1.twinx()
    axes2.spines['left'].set_color('blue')
    axes2.spines['right'].set_color('red')
    axes2.set_ylabel(label2)
    axes2.set_xlim(xlim)
    if limit2 is not None: axes2.set_ylim(limit2)
    axes3 = axes1.twinx()
    self.axes3 = axes3
    if limit2 is not None: axes3.set_ylim(limit2)
    axes2.tick_params('y', color = 'red', labelcolor = 'red')
    axes2.yaxis.label.set_color('red')
    self.curve2, = axes2.plot(freq, d2, color = 'red', label = label2)
    self.curve3, = axes3.plot(freq3, data3, color = 'blue', linewidth = 1, linestyle = '--')
    self.curse = Cursor(axes3, horizOn=True, vertOn=True, useblit=True, color = 'blue', linewidth = 1)
    span = [0,0]
    self.span = span
    #self.crsor.valueChanged.connect(self.plot_crsr)
    #if tf == '': self.canvas.draw()
    self.canvas.draw()
    axes3.callbacks.connect('xlim_changed',self.on_xlims_change)
    self.crsor.valueChanged.connect(self.plot_crsr)
  
  def mark(self):
    self.canvas.draw

  def onclick(self,event):
    #gets cursor coordinates
    ix, iy = event.xdata, event.ydata
    #print ( ix, iy)
    coords = [ix, iy]
    self.coords = coords
 
  def on_xlims_change(self, axes3):
    #gets zoom box coordinates
    #print('xlims_changed')
    k = self.cbox()[2]
    d = self.cbox()[3]
    f = self.cbox()[4]
    self.temp.setValue(sqr2*d[k])
    self.crsor.setValue(f[k])

  def plot_crsr(self):
    #plots vertical line at peak position and can be moved
    #print('plot_crsr')
    # axes1 = self.axes1
    #axes2 = self.axes2
    axes3 = self.axes3
    data4 = self.data4
    if self.cbox() == None: 
      print(self.cbox())
      return
    d = self.cbox()[3]
    hgt = self.cbox()[6]
    print(hgt)
    max = np.max(d)
    # axes1.cla()
    axes3.cla()
    #axes3.tick_params(right = False)
    #axes3.yaxis.set_ticklabels([])
    limi2 = self.limit2
    axes3.set_ylim([-max,max])    
    fre = self.find_dPoint()[1]
    j = self.find_dPoint()[0]
    posf = [fre,fre]
    h = (data4[j])
    z = 0.0
    if self.mode == 'absm':
        #hgt = sqr2*hgt
        axes3.set_ylim(z,hgt[1]) 
        val = [0.0, h]
    if self.mode == 'reim':
        axes3.set_ylim(-hgt[0],hgt[1])
        val = [-h, h]
    self.temp.setValue(hgt[1])
    offset = self.freq[12]-self.freq[0]
    x = axes3.get_ylim()
    axes3.plot(posf, val, color = 'black', linewidth = 1, linestyle = 'solid')
    axes3.text(fre+offset, hgt[1], round(fre,3), bbox=dict(facecolor = 'yellow', alpha = 0.5))
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
      if self.mode == 'absm': 
        self.phase.setValue(0)
        return
      tf = 'set'
      self.tf = tf
      getattr(self, 'plot_%s' % self.mode)()

  def plot_reim(self):
        self.mode = 'reim'
        freq = self.rus.reim.freq
        data = self.rus.reim.data
        d1 = np.real(data)
        d2 = np.imag(data)
        d3 = np.absolute(data)
        max = np.fmax(0.001, d3.max())
        tf = self.tf
        self.plot_curves(freq, d1, 'real', (-1.03 * max, 1.03 * max), d2, 'imag', (-1.03*max, 1.03*max), tf)

  def update_reim(self):
    self.mode = 'reim'
    getattr(self, 'plot_%s' % self.mode)()
    #self.plot_reim()

  def plot_absm(self):
    self.mode = 'absm'
    freq = self.rus.reim.freq
    d3 = np.absolute(self.rus.reim.data)
    max = np.fmax(0.001,d3.max())
    tf = self.tf
    self.plot_curves(freq, d3, ' ', (0, 1.03*max), d3, 'Magnitude', (0, 1.03*max), tf)
    ##############################################################################################################

  def update_absm(self):
    self.mode = 'absm'
    getattr(self, 'plot_%s' % self.mode)()
    #self.plot_reim()

  def cbox(self):
    span = self.axes3.get_xlim()
    hight = self.axes3.get_ylim()
    #print(span, hight)
    if span[0] == 0.0: return
    #print(span)
    d = self.data4
    f = self.freq
    f0 = span[0]
    i = 0
    for fl in self.freq:
        i = i + 1
        if fl >= f0:
                i = i
                break
    f0 = span[1]
    j = 0
    for fr in self.freq:
        j = j + 1
        if fr >= f0:
                  break
    datas = d[i:j]
    freqs = f[i:j]
    k = np.argmax(datas)
    km = k+i
    return fl, fr, k, datas, freqs, span, hight

  def find_dPoint(self):
    f = self.crsor.value()
    i = 0
    for fre in self.freq:
        i = i+1
        if fre >= f: 
            break
    fre = self.freq[i-2]
    j = i-2
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
    self.sweep_Hz = 20
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
    self.sweep_Hz = settings.value('step',20,type = int)
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
