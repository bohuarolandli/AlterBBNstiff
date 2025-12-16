# This is a python interface for alter_BBN for use with stiff matter + GW cosmologies

import os
import numpy as np

class BBN_stiff:
    
    def __init__(self, namestring):
        self.fileBase = namestring
        self.DataDir = '/Users/bohuali/Documents/bbn_codes/alterbbn_v2.2/data' 
        self.FileNameBase = self.fileBase
        self.alterBBN_ExecDir = '/Users/bohuali/Documents/bbn_codes/alterbbn_v2.2'
        self.filename = self.DataDir + '/abundances_' + self.FileNameBase + '.txt'
    
    def calculateAbundances(self, dd0, dd0_rad, ombh2, ndd = 6., 
                            tau = 878.4,   # Neutron lifetime (PDG2022)
                            failsafe = 22, fast = True,
                            tau_err = 0.5, # Neutron lifetime uncertainty (PDG2022)
                           ):

        #ndd = 6.       #Stiff matter energy density exponent
        sd0 = 0.*2.*dd0   #Stiff matter entropy ratio, assuming w=1
        nsd = 0.*5.       #Stiff matter entropy density exponent
        Nnu = 3.044  #Standard model N_eff
        dNnu = 0.    #Extra neutrino species (e.g. sterile neutrinos)
        eta = 273.3036e-10*ombh2  #Baryon to photon ratio from arXiv:1505.01076
        
        if not os.path.exists(self.DataDir):
            os.makedirs(self.DataDir)

        #print(failsafe, fast)
        if fast:
            execString = self.alterBBN_ExecDir + '/alter_stiff.x'
        else: 
            execString = self.alterBBN_ExecDir + '/alter_stiff_mc.x'
        
        execString += ' ' + str(dd0)
        execString += ' ' + str(ndd)
        execString += ' ' + str(sd0)
        execString += ' ' + str(nsd)
        execString += ' ' + str(Nnu)
        execString += ' ' + str(dNnu)
        execString += ' ' + str(dd0_rad)
        execString += ' ' + str(eta)
        execString += ' ' + str(tau)
        execString += ' ' + self.DataDir
        execString += ' ' + self.FileNameBase
        execString += ' 0' # Turn off coupling to plasma (coupd = 0), not affecting entropy but only the expansion rate
        execString += ' 0 0' # Stiff matter energy and entropy cutoff temperatures set to 0, i.e., Tdend = Tsend = 0
        execString += ' ' + str(failsafe) # Switch ODE solver method
        execString += ' ' + str(tau_err)
        #print('*** running: ', execString)
        os.system(execString)
        
        #self.filename = DataDir + '/abundances_' + FileNameBase + '.txt'
        #print(self.filename)
        self.abundances = np.loadtxt(self.filename)
        
        #return abundances
        
    def cleanup(self):
        os.system('rm '+self.filename)
        
