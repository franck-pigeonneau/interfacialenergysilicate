#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 09:18:15 2023

This program determines the surface tension between two immiscible oxide liquids.
The method is based on the development of Kaptay [1]. The chemical potentials of
solute and solvant are defined according to Kim et al. [2].

References:

[1] G. Kaptay (2012) Acta Mater., 60:6804-6813.
[2] S. S. Kim, J. Y. Park & T. H. Sanders Jr (2001) J. Alloys Compd., 321:84-90. 

@author: fpigeonneau
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import zero_Celsius
from scipy import optimize
import pandas as pd
from binarysystem import BinarySystem

def fgamma(TsTc,gamma0,n):
    return gamma0*(1.-TsTc)**n
#end fgamma
    
def fgammavsdx(x,A,n):
    return A*x**n
#end fgammavsdx

# ---------------------------------
# Loading of the database of oxides
# ---------------------------------

solvant='SiO2'
solute='Na2O'
EXPE=True
savefig=True
thermobinary=BinarySystem(solvant,solute)

# ----------------------
# 2. Critical conditions
# ----------------------

dbcritical=pd.read_csv('critcondbinarysystem.csv',index_col=0)
isys=np.argwhere(dbcritical.index==thermobinary.system)[0][0]
Tc=dbcritical['Tc'].values[isys]
xc=dbcritical['xc'].values[isys]
thermobinary.criticalpoint(xc,Tc)
print('Critical conditions: xc=',thermobinary.xc,' at Tc=',thermobinary.Tc,' K')

# --------------------------------
# 3. Determination monotectic line
# --------------------------------

if (thermobinary.nature=='monotectic'):
    dbmono=pd.read_csv('monotecticbinarysystem.csv',index_col=0)
    isys=np.argwhere(dbmono.index==thermobinary.system)[0][0]
    x1mono=dbmono['x1'].values[isys]
    x2mono=dbmono['x2'].values[isys]
    Tmono=dbmono['Tmono'].values[isys]
    thermobinary.monotectic(x1mono,x2mono,Tmono)
    print('x1mono=',thermobinary.x1mono,' x2mono=',thermobinary.x2mono,' Tmono=',thermobinary.Tmono,'K')
else:
    # Here, monotectic conditions do not exist, we determine the limit of the binodale for
    # a temperature below to the critical condition. Nevertheless these values are stored
    # in the same variables that in the previous situation.
    thermobinary.Tmono=thermobinary.Tc-300
    x0=np.array([0.01,0.3])
    thermobinary.x1mono,thermobinary.x2mono=thermobinary.binodale(x0,thermobinary.Tmono)
    print('x1mono=',thermobinary.x1mono,' x2mono=',thermobinary.x2mono,' Tmono=',thermobinary.Tmono-zero_Celsius,'°C')
#end if

# ---------------------------------------
# 4. Determination of the surface tension
# ---------------------------------------

NT=500
Nc=500
T=np.linspace(thermobinary.Tmono,thermobinary.Tc,NT)
x1=np.zeros(NT)
x2=np.zeros(NT)
xi=np.zeros(NT)
x1spi=np.zeros(NT)
x2spi=np.zeros(NT)
gamma=np.zeros(NT)
x1[0]=thermobinary.x1mono
x2[0]=thermobinary.x2mono
x1[NT-1]=thermobinary.xc
x2[NT-1]=thermobinary.xc
x1spi[NT-1]=thermobinary.xc
x2spi[NT-1]=thermobinary.xc

xi[0],gamma[0]=thermobinary.surfacetension(x1[0],x2[0],T[0])
x1spi[0],x2spi[0]=thermobinary.spinodal(x1[0],x2[0],T[0],Nc)
for i in range(1,NT-1):
    # Determination binodal curve
    x0=np.array([x1[i-1],x2[i-1]])
    x1[i],x2[i]=thermobinary.binodale(x0,T[i])
    
    # Determination of the spinodal curve
    x1spi[i],x2spi[i]=thermobinary.spinodal(x1[i],x2[i],T[i],Nc)
    
    # Determination of the surface tension
    xi[i],gamma[i]=thermobinary.surfacetension(x1[i],x2[i],T[i])
#end for
xi[NT-1],gamma[NT-1]=thermobinary.surfacetension(x1[NT-1],x2[NT-1],T[NT-1])

# ---------------------------
# 5. Graphical representation
# ---------------------------

# Fitting of the surface tension
# ------------------------------

popt,pcov=optimize.curve_fit(fgamma,T/thermobinary.Tc,gamma)
print('popt=',popt)
dgamma=pd.DataFrame(np.array([popt]),columns=['A','n'])
dgamma.to_csv('gamma'+thermobinary.system+'.csv')

plt.figure()
plt.plot(T,gamma,'k',label='Kaptay (2012)')
plt.plot(T,fgamma(T/thermobinary.Tc,*popt),'bo',markevery=20,markersize=5,label=r'$\gamma='+\
         str(np.round(popt[0],3))+'(1-T/T_c)^{'+str(np.round(popt[1],3))+'}$')
plt.xlabel(r'$T$ (K)',fontsize=14)
plt.ylabel(r'$\gamma$ (N/m)',fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc=0,fontsize=12)
if (savefig):
    plt.savefig('gamma'+thermobinary.system+'vsT.pdf',dpi=300,bbox_inches='tight')
#endif

# Saving of the results
# ---------------------

A=np.transpose(np.array([T,x1,x2,x1spi,x2spi,xi,gamma]))
data=pd.DataFrame(A,columns=np.array(['T','x1','x2','x1spi','x2spi','xi','gamma']))
data.to_csv('Results/'+thermobinary.system+'.csv')

plt.show()
