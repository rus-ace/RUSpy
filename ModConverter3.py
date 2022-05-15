from math import sqrt
import sys
import os
from pathlib import Path
from PyQt5.uic import loadUiType
from PyQt5.QtCore import QRegExp, QTimer, QSettings, QDir, Qt
from PyQt5.QtWidgets import QApplication, QDoubleSpinBox, QMainWindow, QMessageBox, QDialog, QFileDialog, QPushButton, QLabel, QSpinBox
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSplashScreen

Ui_rus, QMainWindow = loadUiType('ModulusConverter.ui')
class mod(QMainWindow, Ui_rus):
  def __init__(self):
    super(mod, self).__init__()
    self.setupUi(self)
    self.ShearV = 2.0
    #settings = QSettings('ModulusConverter.ini', QSettings.IniFormat)
    self.checkM.stateChanged.connect(self.check_all)
    self.checkB.stateChanged.connect(self.check_all)
    self.checkY.stateChanged.connect(self.check_all)
    self.checkC.stateChanged.connect(self.check_all)
    self.checkS.stateChanged.connect(self.check_all)
    self.checkL.stateChanged.connect(self.check_all)
    self.C44.valueChanged.connect(self.check_all)
    self.Bm.valueChanged.connect(self.check_all)
    self.Ym.valueChanged.connect(self.check_all)
    self.C11.valueChanged.connect(self.check_all)
    self.sg.valueChanged.connect(self.check_all)
    self.la.valueChanged.connect(self.check_all)
    self.rho.valueChanged.connect(self.check_all)

    self.clear.clicked.connect(self.clear_checks)
    self.Reset.clicked.connect(self.zeroBox)
    self.save.clicked.connect(self.write_moduli)

  def clear_checks(self):
    self.checkY.setEnabled(True)
    self.checkC.setEnabled(True)
    self.checkS.setEnabled(True)
    self.checkL.setEnabled(True)
    self.checkB.setEnabled(True)
    self.checkM.setEnabled(True)

    self.checkY.setChecked(False)
    self.checkC.setChecked(False)
    self.checkS.setChecked(False)
    self.checkM.setChecked(False)
    self.checkB.setChecked(False)
    self.checkL.setChecked(False)

  def zeroBox(self):  
    self.C44.setValue(0)
    self.Bm.setValue(0)
    self.Ym.setValue(0)
    self.C11.setValue(0)
    self.sg.setValue(0)
    self.la.setValue(0)
    

  def check_all(self):
    if self.checkM.isChecked() and self.checkB.isChecked(): self.mu_Bm()
    if self.checkM.isChecked() and self.checkY.isChecked(): self.mu_Ym()
    if self.checkM.isChecked() and self.checkC.isChecked(): self.mu_C11()
    if self.checkM.isChecked() and self.checkS.isChecked(): self.mu_sg()
    if self.checkM.isChecked() and self.checkL.isChecked(): self.mu_la()

    if self.checkB.isChecked() and self.checkY.isChecked(): self.Bm_Ym()
    if self.checkB.isChecked() and self.checkC.isChecked(): self.Bm_C11()
    if self.checkB.isChecked() and self.checkS.isChecked(): self.Bm_sg()
    if self.checkB.isChecked() and self.checkL.isChecked(): self.Bm_la()

    if self.checkY.isChecked() and self.checkC.isChecked(): self.Ym_C11()
    if self.checkY.isChecked() and self.checkS.isChecked(): self.Ym_sg()
    if self.checkY.isChecked() and self.checkL.isChecked(): self.Ym_la()

    if self.checkC.isChecked() and self.checkS.isChecked(): self.C11_sg()
    if self.checkC.isChecked() and self.checkL.isChecked(): self.C11_la()

    if self.checkS.isChecked() and self.checkL.isChecked(): self.sg_la()


  def mu_Bm(self):#GK
    mu = self.C44.value()
    Bm = self.Bm.value()
    if(Bm and mu != 0.0):
      rho = self.rho.value()
      la = Bm-2*mu/3
      Ym = mu*(3*la+2*mu)/(la+mu)
      c11 = la + 2*mu
      sg = la/(2*(la+mu))
      self.Ym.setValue(Ym)
      self.sg.setValue(sg)
      self.la.setValue(la)
      self.C11.setValue(c11)
      if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
      self.checkY.setEnabled(False)
      self.checkC.setEnabled(False)
      self.checkS.setEnabled(False)
      self.checkL.setEnabled(False)
    else: return

  def mu_Ym(self):#GE
    mu = self.C44.value()
    Ym = self.Ym.value()
    if(mu and Ym != 0.0):
        rho = self.rho.value()
        la = mu*(Ym-2*mu)/(3*mu-Ym)
        c11 = la + 2*mu
        sg = la/(2*(la+mu))
        Bm = la + 2*mu/3
        self.la.setValue(la)
        self.sg.setValue(sg)
        self.C11.setValue(c11)
        self.Bm.setValue(Bm)
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        #self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        #self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)

        

  def mu_C11(self):#GM
    mu = self.C44.value()
    c11 = self.C11.value()
    if(mu and c11 != 0.0):
        rho = self.rho.value()
        la = c11-2*mu
        Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        sg = la/(2*(la+mu))
        #c11 = la + 2*mu
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        #self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        #self.c11.setValue(c11)
        self.la.setValue(la)
        self.sg.setValue(sg)

        #self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        #self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)


  def mu_sg(self):#Gv
    mu = self.C44.value()
    sg = self.sg.value()
    la = 2*mu*sg/(1-2*sg)
    if(mu != 0.0):
        rho = self.rho.value()
        Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        #self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        #self.C11.setValue(c11)
        self.la.setValue(la)
        self.sg.setValue(sg)
        #self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        #self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)

  def mu_la(self):#Gl
    mu = self.C44.value()
    la = self.la.value()
    if(mu and la != 0.0):
        rho = self.rho.value()
        Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        #self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        self.C11.setValue(c11)
        #self.la.setValue(la)
        self.sg.setValue(sg)
        #self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        #self.checkL.setEnabled(False)

  def Bm_Ym(self):#KE
    Bm = self.Bm.value()
    Ym = self.Ym.value()
    mu = 3*Bm*Ym/(9*Bm-Ym)
    la = 3*Bm*(3*Bm-Ym)/(9*Bm-Ym)
    if(Ym and Bm != 0.0):
        rho = self.rho.value()
        #Bm = la + 2*mu/3
        #Ym = mu*(3*la+2*mu)/(la+mu)
        c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        #self.Bm.setValue(Bm)
        #self.Ym.setValue(Ym)
        self.C11.setValue(c11)
        self.la.setValue(la)
        self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        #self.checkB.setEnabled(False)
        #self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)



  def Bm_C11(self):#KM
    Bm = self.Bm.value()
    c11 = self.C11.value()   
    mu = 3*(c11-Bm)/4
    la = (3*Bm-c11)/2
    if(Bm and c11 != 0.0):
        rho = self.rho.value()
        #Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        #c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        #self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        #self.C11.setValue(c11)
        self.la.setValue(la)
        self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        #self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        #self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)

  def Bm_sg(self):
    Bm = self.Bm.value()
    sg = self.sg.value()
    mu = 3*Bm*(1-2*sg)/(2*(1+sg))
    la = 3*Bm*sg/(1+sg)
    if(Bm != 0.0):
        rho = self.rho.value()
        #Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        c11 = la + 2*mu
        #sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        #self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        self.C11.setValue(c11)
        self.la.setValue(la)
        #self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        #self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        #self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)

  def Bm_la(self):#Kla
    Bm = self.Bm.value()
    la = self.la.value()
    mu = 3*(Bm-la)/2
    if(Bm and la != 0.0):
        rho = self.rho.value()
        #Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        #self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        self.C11.setValue(c11)
        self.la.setValue(la)
        #self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        #self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        #self.checkL.setEnabled(False)


  def Ym_C11(self):#EM
    Ym = self.Ym.value()
    c11 = self.C11.value()
    s = sqrt(Ym*Ym + 9*c11*c11 - 10*Ym*c11)
    mu = (3*c11+Ym-s)/8
    la = (c11-Ym+s)/4
    if(Ym and c11 != 0.0):
        rho = self.rho.value()
        Bm = la + 2*mu/3
        #Ym = mu*(3*la+2*mu)/(la+mu)
        #c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        #self.Ym.setValue(Ym)
        #self.C11.setValue(c11)
        self.la.setValue(la)
        self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        #self.checkY.setEnabled(False)
        #self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)


  def Ym_sg(self):#Ev
    Ym = self.Ym.value()
    sg = self.sg.value()
    mu = Ym/(2*(1+sg))
    la = Ym*sg/((1+sg)*(1-2*sg))
    if(Ym != 0.0):
        rho = self.rho.value()
        Bm = la + 2*mu/3
        #Ym = mu*(3*la+2*mu)/(la+mu)
        c11 = la + 2*mu
        #sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        #self.Ym.setValue(Ym)
        self.C11.setValue(c11)
        self.la.setValue(la)
        self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        #self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        #self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)

  def Ym_la(self):#Ela
    Ym = self.Ym.value()
    la = self.la.value()
    Rs = sqrt(Ym*Ym + 9*la*la+2*Ym*la)
    mu = (Ym-3*la+Rs)/4
    if(Ym and la != 0.0):
        rho = self.rho.value()
        Bm = la + 2*mu/3
        #Ym = mu*(3*la+2*mu)/(la+mu)
        c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        #self.Ym.setValue(Ym)
        self.C11.setValue(c11)
        #self.la.setValue(la)
        self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        #self.checkY.setEnabled(False)
        self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        #self.checkL.setEnabled(False)


  def C11_sg(self):#Mv
    c11 = self.C11.value()
    sg = self.sg.value()
    mu = c11*(1-2*sg)/(2*(1-sg))
    la = c11*sg/(1-sg)
    if(c11 != 0.0):
        rho = self.rho.value()
        Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        #c11 = la + 2*mu
        #sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        #self.C11.setValue(c11)
        self.la.setValue(la)
        #self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        #self.checkC.setEnabled(False)
        #self.checkS.setEnabled(False)
        self.checkL.setEnabled(False)

  def C11_la(self):#Mla
    c11 = self.C11.value()
    la = self.la.value()
    mu = (c11-la)/2
    if(c11 and la != 0.0):
        rho = self.rho.value()
        Bm = la + 2*mu/3
        Ym = mu*(3*la+2*mu)/(la+mu)
        #c11 = la + 2*mu
        sg = la/(2*(la+mu))
        if rho != 0.0:
            self.vl.setValue(sqrt(c11/rho))
            self.vs.setValue(sqrt(mu/rho))
        self.C44.setValue(mu)
        self.Bm.setValue(Bm)
        self.Ym.setValue(Ym)
        #self.C11.setValue(c11)
        #self.la.setValue(la)
        self.sg.setValue(sg)
        self.checkM.setEnabled(False)
        self.checkB.setEnabled(False)
        self.checkY.setEnabled(False)
        #self.checkC.setEnabled(False)
        self.checkS.setEnabled(False)
        #self.checkL.setEnabled(False)

  def sg_la(self):#Gv
    sg = self.sg.value()
    if la !=0 and mu!= 0:
      la = self.la.value()
      mu = la*((1-sg)/(2*sg))
    else: return
    rho = self.rho.value()
    Bm = la + 2*mu/3
    Ym = mu*(3*la+2*mu)/(la+mu)
    c11 = la + 2*mu
    #sg = la/(2*(la+mu))
    if rho != 0.0:
        self.vl.setValue(sqrt(c11/rho))
        self.vs.setValue(sqrt(mu/rho))
    self.C44.setValue(mu)
    self.Bm.setValue(Bm)
    self.Ym.setValue(Ym)
    self.C11.setValue(c11)
    #self.la.setValue(la)
    #self.sg.setValue(sg)
    self.checkM.setEnabled(False)
    self.checkB.setEnabled(False)
    self.checkY.setEnabled(False)
    self.checkC.setEnabled(False)
    #self.checkS.setEnabled(False)
    #self.checkL.setEnabled(False)

  def write_moduli(self):
    c44 = str('Shear modulus(GPa) ' + '\u03bc' + ' = ' + str(round(self.C44.value(),2))+'\n')
    Ym = str('Youngs modulus(GPa) = ' + str(round(self.Ym.value(),2))+'\n')
    Bm = str('Bulk modulus(GPa) = ' + str(round(self.Bm.value(),2))+'\n')
    c11 = str('C11(Gpa) = ' + str(round(self.C11.value(),2))+'\n')
    sg = str('Poisons Ratio = ' + str(round(self.sg.value(),4))+'\n')
    la = str('lame constant(GPa) ' + '\u03bb' + ' = '+ str(round(self.la.value(),2))+'\n')
    dialog = QFileDialog(self, 'Save results', '.', '*.dat')
    dialog.setDefaultSuffix('ini')
    #dialog.selectFile('rus.ini')
    dialog.setAcceptMode(QFileDialog.AcceptSave)
    dialog.setOptions(QFileDialog.DontConfirmOverwrite)
    if dialog.exec() == QDialog.Accepted:
        name = dialog.selectedFiles()
        fd = os.path.join(path,name[0])
        fh = open(fd, 'w',encoding="utf-8" )
        fh.write(c44)
        fh.write(Bm)
        fh.write(Ym)
        fh.write(c11)
        fh.write(sg)
        fh.write(la)
        fh.close()

path = os.getcwd()
app = QApplication(sys.argv)
app.processEvents()
window = mod()
window.show()
#splash.finish(window)
sys.exit(app.exec_())