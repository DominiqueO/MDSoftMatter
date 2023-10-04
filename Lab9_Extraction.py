# Molecular Simulations lab9 (Free Energy) Extraction Program
# Author: Dominique Ostermayer

from __future__ import print_function, division
import matplotlib.pyplot as plt
import re
import numpy as np
from numpy import trapz
import matplotlib.axes as ax
import math 
from scipy.integrate import simps
from scipy import stats
import os 
import fnmatch as fn
import pandas as pd
import statsmodels.api as sm


directory = "" # working directory format "x:/y/z/"

sig = 'dab2' # signature/name of file
sig2 = 'ab2'

format = '.dat' # format of output files to be analysed

filelist = []
TI_EPRTOT = []
TI_FOR = []
TP_EPRTOT = []
TP_FOR = []

for filename in os.listdir(directory):
    # search for the output files with prefix ab2 or dab2 and ending with .dat
    if filename.startswith(sig) and filename.endswith(format) and not filename.endswith('inter.dat') or filename.startswith(sig2) and filename.endswith(format) and not filename.endswith('inter.dat'):
        path_filename = os.path.join(directory, filename)
        filelist.append(filename)
        in_file = open(path_filename, "rt") # open file .dat (output file) 
        contents = in_file.read()         # read the entire file into a string variable
        in_file.close()  # close the file
        LAMBD = re.findall("LAMBDA= (.*)\D", contents) # extract value of LAMBDA
        LAMBD = LAMBD[1:] # remove 0 from LAMBD
        LAMBD = [float(x) for x in LAMBD]
        match =  re.findall(r'PERTURBATION> TI.*EPRTOT= (.*)EFORWARD=', contents)
        match = [float(x) for x in match]
        TI_EPRTOT.append([filename, LAMBD, match])
        match = re.findall(r'PERTURBATION> TI.*EFORWARD=\s*(-*\d*.\d+)\D+', contents)
        match = [float(x) for x in match]
        TI_FOR.append([filename, LAMBD,  match])
        match =  re.findall(r'PERTURBATION> TP.*EPRTOT= (.*)EFORWARD=', contents)
        match = [float(x) for x in match]
        TP_EPRTOT.append([filename, LAMBD, match])
        match = re.findall(r'PERTURBATION> TP.*EFORWARD=\s*(-*\d*.\d+)\D+.*EBACKWARD=\s*(-*\d*.\d+)\D+', contents)
        sumlist = []
        for i in range(len(match)):
            x = float(match[i][0])-float(match[i][1])
            sumlist.append(x)
        TP_FOR.append([filename, LAMBD, sumlist])
        
        
# Forward Free Energy Plots 
for i in range(len(TI_FOR)):
    filename = 'dab2ava.dat'
    if TI_FOR[i][0] == filename:
        ti = np.array(TI_FOR[i][2])
        tp = np.array(TP_FOR[i][2])
        tplam = np.array(TI_FOR[i][1])
        tilam = np.array(TP_FOR[i][1])
        tinorminc = ti
        tpnorminc = tp
        lamnorminc = tplam

err_for_incr = (tp[-1]-ti[-1])/tp[-1]

                
plt.figure('Incremental Free Energy')   
plt.title('Incremental Free Energy')  
plt.plot(tilam, ti, label='Thermodynamic Integration (forward)')
plt.plot(tplam, tp, label='Exponential Formula (forward)')
plt.plot(tilam, np.cumsum(ti), label='Thermodynamic Integration (forward), (Cumulative Sum)')
plt.plot(tplam, np.cumsum(tp), label='Exponential Formula (forward), (Cumulative Sum)')

        
for i in range(len(TI_FOR)):
    filename = 'dab2ava_b.dat'
    if TI_FOR[i][0] == filename:
        ti = np.array(TI_FOR[i][2])
        tp = np.array(TP_FOR[i][2])
        tplam = np.array(TI_FOR[i][1])
        tilam = np.array(TP_FOR[i][1])
        
err_back_incr = (tp[-1]-ti[-1])/tp[-1]


plt.plot(tilam, ti, label='Thermodynamic Integration (backward)')
plt.plot(tplam, tp, label='Exponential Formula (backward)')
plt.plot(tilam, np.cumsum(ti), label='Thermodynamic Integration (backward), (Cumulative Sum)')
plt.plot(tplam, np.cumsum(tp), label='Exponential Formula (backward), (Cumulative Sum)')
plt.xlabel(r'$\lambda$')
plt.ylabel('Free Energy (kcal/mol)')
plt.legend()
plt.show()    

    
# Accumulated Free Energy Plots Plots
for i in range(len(TI_EPRTOT)):
    filename = 'dab2ava.dat'
    if TI_EPRTOT[i][0] == filename:
        ti = np.array(TI_EPRTOT[i][2])
        tp = np.array(TP_EPRTOT[i][2])
        tplam = np.array(TI_EPRTOT[i][1])
        tilam = np.array(TP_EPRTOT[i][1])
        tinormacc = ti
        tpnormacc = tp
        
err_for_eprtot = (tp[-1]-ti[-1])/tp[-1]
        
plt.figure('Accumulated Free Energy')   
plt.title('Accumulated Free Energy')  
plt.plot(tilam, ti, label='Thermodynamic Integration (forward)')
plt.plot(tplam, tp, label='Exponential Formula (forward)')
# plt.plot(tilam, np.cumsum(ti), label='Thermodynamic Integration (forward), (Cumulative Sum)')
# plt.plot(tplam, np.cumsum(tp), label='Exponential Formula (forward), (Cumulative Sum)')


        
for i in range(len(TI_EPRTOT)):
    filename = 'dab2ava_b.dat'
    if TI_EPRTOT[i][0] == filename:
        ti = np.array(TI_EPRTOT[i][2])
        tp = np.array(TP_EPRTOT[i][2])
        tplam = np.array(TI_EPRTOT[i][1])
        tilam = np.array(TP_EPRTOT[i][1])

err_back_eprtot = (tp[-1]-ti[-1])/tp[-1]


plt.plot(tilam, ti, label='Thermodynamic Integration (backward)')
plt.plot(tplam, tp, label='Exponential Formula (backwards)')
# plt.plot(tilam, np.cumsum(ti), label='Thermodynamic Integration (backwards), (Cumulative Sum)')
# plt.plot(tplam, np.cumsum(tp), label='Exponential Formula (backwards), (Cumulative Sum)')
plt.xlabel(r'$\lambda$')
plt.ylabel('Free Energy (kcal/mol)')
plt.legend()
plt.show()        

# Stage 6 comparison of different seed numbers
seedlist = []
for i in range(len(filelist)):
    filename = 'dab2ava.+\d' # name of file, matches all files, beginning by dab2ava followed by an arbitrary number: example dab2ava_3.dat or dab2ava_seed3.dat
    if re.match(filename, TI_EPRTOT[i][0]):
        seedno = re.match(r'dab2ava_.*(\d+)', TI_EPRTOT[i][0])
        seedno = int(seedno.group(1))
        ti = TI_EPRTOT[i][2]
        tp = TP_EPRTOT[i][2]
        error = (tp[-1]-ti[-1])/tp[-1]
        seedlist.append([seedno, tp, ti, error]) # results are safed to list named seedlist.
        
# stage 7 different integration paths
# pathways Incremental Free Energy Plots Plots
for i in range(len(TI_FOR)):
    filename = 'dab2ab.dat'
    if TI_FOR[i][0] == filename:
        ti = np.array(TI_FOR[i][2])
        tp = np.array(TP_FOR[i][2])
        tplam = np.array(TI_FOR[i][1])
        tilam = np.array(TP_FOR[i][1])

err_for_incr = (tp[-1]-ti[-1])/tp[-1]

                
plt.figure('Incremental Free Energy (Pathways)')   
plt.title('Incremental Free Energy')  
plt.plot(tilam, ti, label='Thermodynamic Integration (dab2ab)')
plt.plot(tplam, tp, label='Exponential Formula (dab2ab)')
plt.plot(tilam, np.cumsum(ti), label='Thermodynamic Integration (dab2ab), (Cumulative Sum)')
plt.plot(tplam, np.cumsum(tp), label='Exponential Formula (dab2ab), (Cumulative Sum)')
plt.xlabel('LAMBDA')
plt.ylabel('Incremental Free Energy (kcal/mol)')
x = lamnorminc

plt.plot(x, tinorminc, label='Thermodynamic Integration (dab2ava)')
plt.plot(x, tpnorminc, label='Exponential Formula (dab2ava)')
plt.plot(x, np.cumsum(tinorminc), label='Thermodynamic Integration (dab2ava), (Cumulative Sum)')
plt.plot(x, np.cumsum(tpnorminc), label='Exponential Formula (dab2ava), (Cumulative Sum)')

plt.legend()
plt.show()
        
for i in range(len(TI_FOR)):
    filename = 'ab2ava.dat'
    if TI_FOR[i][0] == filename:
        ti = np.array(TI_FOR[i][2])
        tp = np.array(TP_FOR[i][2])
        tplam = np.array(TI_FOR[i][1])
        tilam = np.array(TP_FOR[i][1])
        
err_back_incr = (tp[-1]-ti[-1])/tp[-1]


plt.plot(tilam, ti, label='Thermodynamic Integration (ab2ava)')
plt.plot(tplam, tp, label='Exponential Formula (ab2ava)')
plt.plot(tilam, np.cumsum(ti), label='Thermodynamic Integration (backwards), (ab2ava)')
plt.plot(tplam, np.cumsum(tp), label='Exponential Formula (backwards), (ab2ava)')
plt.xlabel(r'$\lambda$')
plt.ylabel('Free Energy (kcal/mol)')
plt.legend()
plt.show()    

    
# Stage 7 pathways Accumulated Free Energy Plots Plots 
for i in range(len(TI_FOR)):
    filename = 'dab2ab.dat'
    if TI_EPRTOT[i][0] == filename:
        ti = np.array(TI_EPRTOT[i][2])
        tp = np.array(TP_EPRTOT[i][2])
        tplam = np.array(TI_EPRTOT[i][1])
        tilam = np.array(TP_EPRTOT[i][1])
        
err_for_eprtot = (tp[-1]-ti[-1])/tp[-1]
        
plt.figure('Accumulated Free Energy (Pathways)')   
plt.title('Accumulated Free Energy')  
plt.plot(tilam, ti, label='Thermodynamic Integration (dab2ab)')
plt.plot(tplam, tp, label='Exponential Formula (dab2ab)')
plt.xlabel(r'$\lambda$')
plt.ylabel('Free Energy (kcal/mol)')
plt.legend()
plt.show()
        
for i in range(len(TI_FOR)):
    filename = 'ab2ava.dat'
    if TI_EPRTOT[i][0] == filename:
        ti = np.array(TI_EPRTOT[i][2])
        tp = np.array(TP_EPRTOT[i][2])
        tplam = np.array(TI_EPRTOT[i][1])
        tilam = np.array(TP_EPRTOT[i][1])

err_back_eprtot = (tp[-1]-ti[-1])/tp[-1]


plt.plot(tilam, ti, label='Thermodynamic Integration (ab2ava)')
plt.plot(tplam, tp, label='Exponential Formula (ab2ava)')


plt.plot(x, tinormacc, label='Thermodynamic Integration (dab2ava)')
plt.plot(x, tpnormacc, label='Exponential Formula (dab2ava)')
plt.plot(x, np.cumsum(tinormacc), label='Thermodynamic Integration (dab2ava), (Cumulative Sum)')
plt.plot(x, np.cumsum(tpnormacc), label='Exponential Formula (dab2ava), (Cumulative Sum)')



plt.xlabel('LAMBDA')
plt.ylabel('Accumulated Free Energy (kcal/mol)')
plt.legend()
plt.show()        

 
  
    
        

