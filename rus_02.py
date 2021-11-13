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
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtWidgets import QApplication, QDoubleSpinBox, QMainWindow, QMessageBox, QDialog, QFileDialog, QPushButton, QLabel, QSpinBox
from PyQt5.QtNetwork import QAbstractSocket, QTcpSocket

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
  def __init__(self, layout, rus):
    # create figure
    self.figure = Figure()
    if sys.platform != 'win32':
      self.figure.set_facecolor('none')
    self.canvas = FigureCanvas(self.figure)
    layout.addWidget(self.canvas)
    # create navigation toolbar
    self.toolbar = NavigationToolbar(self.canvas, None, False)
    self.toolbar.layout().setSpacing(6)
    # remove subplots action
    #actions = self.toolbar.actions()
    layout.addWidget(self.toolbar)
    self.tf = ''
    ##################################################################################
    #add cursor
    self.cursor = QDoubleSpinBox()
    self.cursor.setDecimals(3)
    self.toolbar.addWidget(self.cursor)
    #add value at cursor
    self.temp = QDoubleSpinBox()
    self.temp.setReadOnly(True)
    self.temp.setButtonSymbols(2)
    self.temp.setDecimals(6)
    self.toolbar.addWidget(self.temp)
    #add set phase
    self.mode = 'dut'
    self.phase = QSpinBox()
    self.phase.setSingleStep(5)
    self.phase.setMaximum(180)
    self.phase.setMinimum(-180)
    if self.mode == 'dut': self.phase.valueChanged.connect(self.phase_set)
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
    self.plotButton.setStyleSheet('background-color:blue')
    self.plotButton.clicked.connect(self.plot)
    self.toolbar.addWidget(self.plotButton)
##################################################################################
    layout.addWidget(self.toolbar)
    #self.figure.subplots_adjust(left = 0.09, bottom = 0.1, right = 0.91, top = 0.96)
    axes1 = self.figure.add_subplot(1,1,1)
    axes3 = self.figure.add_subplot(1,1,1)
    self.axes1 = axes1
    self.axes3 = axes3
    self.mode = 'dut'
    self.rus = rus

  def xlim(self, freq):
    start = freq[0]
    stop = freq[-1]
    min = np.minimum(start, stop)
    max = np.maximum(start, stop)
    margin = (max - min) / 50
    return (min - margin, max + margin)

  def plot(self):
    self.cursor.valueChanged.connect(self.plot_cursor)
    getattr(self, 'plot_%s' % self.mode)()


  def update(self, mode):
    self.cursor.valueChanged.connect(self.plot_cursor)
    getattr(self, 'update_%s' % mode)()

####################################################################################
####################################################################################
  
  def plot_curves(self, freq, data1, label1, limit1, data2, label2, limit2, tf):
    matplotlib.rcdefaults()
    matplotlib.rcParams['axes.formatter.use_mathtext'] = True
    self.figure.clf()
    data4 = np.abs(data1 +1j*data2)
    ph = 0
    phD = 0
    co = 1
    si = 0
    #self.tf = ''
    if tf!='':
      if self.mode == 'dut':
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
    d1 = co*data1 - si*data2
    d2 = si*data1 + co*data2
    self.freq = freq
    self.limit2 = limit2
    self.data4 = data4
    #find maximum
    i = np.argmax(data4)
    freqmax = float(freq[i])
    self.vmax = float(data1[i])
    freq3 = [freqmax,freqmax] #############
    data3 = [-data4[i], data4[i]]
    #add axes
    #self.figure.subplots_adjust(left = 0.09, bottom = 0.1, right = 0.91, top = 0.96)
    axes1 = self.figure.add_subplot(1,1,1)
    #axes1 = self.axes1
    axes1.cla()
    axes1.xaxis.grid()
    axes1.set_xlabel('kHz')
    axes1.set_ylabel(label1)
    xlim = self.xlim(freq)
    axes1.set_xlim(xlim)
    ###################################################################################
    self.data3 = data3
    #Setup cursor spinbox
    self.cursor.setMaximum(freq[-1])
    self.cursor.setMinimum(freq[0])
    sstep = (freq[1]-freq[0])
    self.cursor.setSingleStep(sstep)
    self.cursor.setValue(freqmax)
    self.temp.setValue(data4[i])
    #############################################################################
    #self.cursor.valueChanged.connect(self.plot_cursor)
    #########################################################################
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
    #axes2.spines['left'].set_linewidth(1)
    #axes2.spines['right'].set_linewidth(1)
    axes2.set_ylabel(label2)
    axes2.set_xlim(xlim)
    axes3 = axes1.twinx()
    axes3.tick_params(right = False)
    axes3.yaxis.set_ticklabels([])
    self.axes3 = axes3
    if limit2 is not None: axes2.set_ylim(limit2)
    axes2.tick_params('y', color = 'red', labelcolor = 'red')
    axes2.yaxis.label.set_color('red')
    self.curve2, = axes2.plot(freq, d2, color = 'red', label = label2)
    self.curve3, = axes3.plot(freq3, data3, color = 'green', linewidth = 2, linestyle = 'dashed')
    if tf == '': self.canvas.draw()
  
  
  def plot_cursor(self):
    axes3 = self.axes3
    axes3.cla()
    axes3.tick_params(right = False)
    axes3.yaxis.set_ticklabels([])
    limi2 = self.limit2
    if limi2 is not None: axes3.set_ylim(limi2)
    f = self.cursor.value()
    i = 1
    for fre in self.freq:
        i = i + 1
        if fre >= f:
                i = i - 1
                break
    posf = [fre,fre]
    hgt = self.data4[i]
    val = [-hgt, hgt]
    offset = self.freq[12]-self.freq[0]
    self.temp.setValue(hgt)
    #self.curve3, = self.axes3.plot(posf, val, color = 'black', linewidth = 2, linestyle = 'solid')
    axes3.plot(posf, val, color = 'black', linewidth = 2, linestyle = 'solid')
    axes3.text(fre+offset, hgt/2, round(fre,3), bbox=dict(facecolor = 'yellow', alpha = 0.5))

    self.plot_absm
    self.canvas.draw()

  def phase_adj(self, tf):
    #self.mode = 'dut'
    self.tf = 'adj'
    getattr(self, 'plot_%s' % self.mode)()
  

  def phase_set(self, tf):
      if self.mode == 'absm': 
        self.phase.setValue(0)
        return
      self.tf = 'set'
      getattr(self, 'plot_%s' % self.mode)()

  def plot_re_im(self, freq, data, label, mode):
        #self.mode = 'dut'
        max = np.fmax(0.001, np.abs(data).max())
        label1 = 'real'
        label2 = 'imag'
        tf = self.tf
        #self.curve3 = self.plot_curves(freq, np.real(data), label1, (-1.03 * max, 1.03 * max), np.imag(data), label2, (-1.03*max, 1.03*max), tf)
        self.plot_curves(freq, np.real(data), label1, (-1.03 * max, 1.03 * max), np.imag(data), label2, (-1.03*max, 1.03*max), tf)

  def update_re_im(self, freq, data, label, mode):
     self.plot_re_im(freq, data, label, mode)

  def plot_dut(self):
    self.plot_re_im(self.rus.dut.freq, self.rus.dut.data, 'dut', 'dut')

  def update_dut(self):
    self.update_re_im(self.rus.dut.freq, self.rus.dut.data, 'dut', 'dut')

  def plot_absm(self):
    self.mode = 'absm'
    freq = self.rus.dut.freq
    d3 = np.absolute(self.rus.dut.data)
    max = np.fmax(0.001,d3.max())
    tf = self.tf
    self.plot_curves(freq, d3, ' ', (0, 1.03*max), d3, 'Magnitude', (0, 1.03*max), tf)
    ##############################################################################################################

  def update_absm(self):
    #self.mode = 'absm'
    self.plot_absm()

class rus(QMainWindow, Ui_rus):
  graphs = ['dut', 'absm']

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
    self.sweep_start = 336
    self.sweep_stop = 340
    self.sweep_Hz = 2
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/self.sweep_Hz)
    if(self.sweep_size)>32766:
      self.sweep_size = 32766
    # buffer and offset for the incoming samples
    self.buffer = bytearray(16 * 32768)
    self.offset = 0
    self.data = np.frombuffer(self.buffer, np.complex64)
    # create measurements
    self.dut = Measurement(self.sweep_start, self.sweep_stop, self.sweep_size)
    ######################################################################################################
    self.mode = 'dut'   
    # create figures
    self.tabs = {}
    for i in range(len(self.graphs)):
      layout = getattr(self, '%sLayout' % self.graphs[i])
      self.tabs[i] = FigureTab(layout, self)
    # configure widgets
    self.rateValue.addItems(['5000', '1000', '500', '100', '50', '10', '5', '1'])
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
    self.singleSweep.clicked.connect(partial(self.sweep, 'dut'))
    self.autoSweep.clicked.connect(self.sweep_auto)
    self.stopSweep.clicked.connect(self.cancel)
    self.datButton.clicked.connect(self.write_dat)
    self.startValue.valueChanged.connect(self.set_start)
    self.stopValue.valueChanged.connect(self.set_stop)
    self.sizeValue.valueChanged.connect(self.set_size)
    self.rateValue.currentIndexChanged.connect(self.set_rate)
    self.level1Value.valueChanged.connect(self.set_level1)
    self.tabWidget.currentChanged.connect(self.update_tab)
    # create timers
    self.startTimer = QTimer(self)
    self.startTimer.timeout.connect(self.timeout)
    self.sweepTimer = QTimer(self)
    self.sweepTimer.timeout.connect(self.sweep_timeout)

  def set_enabled(self, enabled):
    widgets = [self.rateValue, self.level1Value, self.startValue, self.stopValue, self.sizeValue, self.singleSweep, self.autoSweep]
    
    for entry in widgets:
      entry.setEnabled(enabled)

  def start(self):
    if self.idle:
      self.connectButton.setEnabled(False)
      self.socket.connectToHost(self.addrValue.text(), 1001)
      self.startTimer.start(5000)
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
      limit = 16 * self.sweep_size
      if self.offset + size < limit:
        self.buffer[self.offset:self.offset + size] = self.socket.read(size)
        self.offset += size
      else:
        self.buffer[self.offset:limit] = self.socket.read(limit - self.offset)
        adc1 = self.data[0::2]
        #adc2 = self.data[1::2]
        attr = getattr(self, 'dut')
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
    self.sweep_size = int(1000*(self.sweep_stop-self.sweep_start)/value)
    if self.sweep_size > 32766:
      self.sweep_size = 32766
    #print(self.sweep_size)

  def set_rate(self, value):
    if self.idle: return
    rate = [10, 50, 100, 500, 1000, 5000, 10000, 50000][value]
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
      self.sweep('dut')

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
    settings.setValue('dut_start', int(self.dut.freq[0]))
    settings.setValue('dut_stop', int(self.dut.freq[-1]))
    settings.setValue('dut_size', self.dut.freq.size)
    settings.setValue('step',self.sweep_Hz)
  
  def read_cfg_settings(self, settings):
    self.addrValue.setText(settings.value('addr', '192.168.1.100'))
    self.rateValue.setCurrentIndex(settings.value('rate', 0, type = int))
    self.level1Value.setValue(settings.value('level_1', 2, type = int))
    dut_start = settings.value('dut_start', 336, type = int)
    dut_stop = settings.value('dut_stop', 340, type = int)
    dut_size = settings.value('dut_size', 1, type = int)
    self.sweep_Hz = settings.value('step',2,type = int)
    self.startValue.setValue(dut_start)
    self.stopValue.setValue(dut_stop)
    self.sizeValue.setValue(dut_size)
    self.sizeValue.setValue(self.sweep_Hz)
    
  def write_dat(self):
    dialog = QFileDialog(self, 'Write dat file', '.', '*.dat')
    dialog.setDefaultSuffix('dat')
    dialog.setAcceptMode(QFileDialog.AcceptSave)
    dialog.setOptions(QFileDialog.DontConfirmOverwrite)
    if dialog.exec() == QDialog.Accepted:   
      name = dialog.selectedFiles()
      fh = open(name[0], 'w')
      f = self.dut.freq
      d = self.dut.data
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
