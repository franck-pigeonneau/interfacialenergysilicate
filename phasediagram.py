#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat March 14 2026

This program determines the miscibility gap in binary silicate glasses. The
chemical potentials of solute and solvant are defined according to Kim et al. [1].

Experimental data are used to compare with the numerical results. All methods are
gathered in the class BinarySystem.

References:

[1] S. S. Kim, J. Y. Park & T. H. Sanders Jr (2001) J. Alloys Compd., 321:84-90. 

@author: Franck Pigeonneau
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import zero_Celsius
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
Ref='Toropov et al. (1961)'
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
    thermobinary.Tmono=thermobinary.Tc-400
    #thermobinary.Tmono=1520.
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
#end for

# ---------------------------
# 5. Graphical representation
# ---------------------------

plt.figure()
plt.plot(x1,T,'k-',label='Binodal')
plt.plot(x2,T,'k-')
plt.plot(x1spi,T,'r-.',label='Spinodal')
plt.plot(x2spi,T,'r-.')
if (EXPE):
    fileexpe='BinodalData/xvsT-'+thermobinary.system+'.csv'
    dataexpe=pd.read_csv(fileexpe)
    plt.plot(dataexpe['x'].values,dataexpe['T'].values,'ko',label=thermobinary.ref)
#endif
plt.xlabel(r'$c$',fontsize=14)
plt.ylabel(r'$T$ (K)',fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc=0,fontsize=10)
plt.annotate(str(np.round(thermobinary.xc,3))+','+str(np.round(thermobinary.Tc,2)),\
             (thermobinary.xc,thermobinary.Tc))
if (savefig):
    plt.savefig('gapimmisible'+thermobinary.system+'vsT.png',dpi=300,bbox_inches='tight')
#endif

plt.show()
