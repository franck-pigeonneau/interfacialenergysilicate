#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 18:13:01 2023

This program determines the excess potential of Redlich-Kister upto the fourth
order. The coefficients A, function as a temperature, are taken from a function
given in [1]. The function to minimize is determined in a method of the class 
BinarySystem. The Jacobian vertor is also computed.

Reference:

[1] F. Pigeonneau & W. Blanc (2026). Interfacial energy in phase separated borate 
and silicate systems J. Am. Ceram. Soc., 10.1111/JACE.71109.


@author: Franck Pigeonneau
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import pandas as pd
from binarysystem import BinarySystem

# ----------------------------------------
# Definition of the binary system to study
# ----------------------------------------

solvent='B2O3'
solute='MgO'
OrderRK=4
NTRK=3
Optim=False
thermobinary=BinarySystem(solvent,solute,OrderRK,NTRK)

# --------------------------------------
# Experimental data of the binodal curve
# --------------------------------------

data=pd.read_csv('BinodalData/'+thermobinary.system+'.csv',index_col=0)
Texpe=data['T'].values
x1expe=data['x1'].values
x2expe=data['x2'].values

# ---------------------------------------------------
# Loading of the experimental critical condition data
# ---------------------------------------------------

dbcritical=pd.read_csv('critcondbinarysystem.csv',index_col=0)
isys=np.argwhere(dbcritical.index==thermobinary.system)[0][0]
Tcexpe=dbcritical['Tc'].values[isys]
xcexpe=dbcritical['xc'].values[isys]

# -----------------------------------------------------
# Loading of the experimental monotectic condition data
# -----------------------------------------------------

if (thermobinary.nature=='monotectic'):
    dbmono=pd.read_csv('monotecticbinarysystem.csv',index_col=0)
    isys=np.argwhere(dbmono.index==thermobinary.system)[0][0]
    x1mono=dbmono['x1'].values[isys]
    x2mono=dbmono['x2'].values[isys]
    Tmono=dbmono['Tmono'].values[isys]
else:
    x1mono=data['x1'].values[0]
    x2mono=data['x2'].values[0]
    Tmono=data['T'].values[0]
#endif

# -----------------------------------------------------------------------------------------
# Determination of the coefficient A of the Redlich-Kister potential by lesat-square method
# -----------------------------------------------------------------------------------------

if (Optim):
    # Initial values of A
    A=thermobinary.A.copy()
    Ninc=(OrderRK+1)*NTRK

    # Optimization procedure
    x0=np.zeros(Ninc)
    for i in range(OrderRK+1):
        for j in range(NTRK):
            x0[j+NTRK*i]=thermobinary.A[i,j]
        #end for
    #end for
    
    weightmonotectic=1.e0
    
    xRK=optimize.minimize(thermobinary.fRK,x0,args=(Texpe,x1expe,x2expe,xcexpe,Tcexpe,x1mono,x2mono,Tmono,weightmonotectic),\
                          method='COBYQA',jac=thermobinary.jacfRK,tol=1.e-14,options={'maxiter': 100000})
#endif

# Copy of the new coefficients of the Redlich-Kister
print('A=',thermobinary.A)

# ------------------------------------
# Determination of critical conditions
# ------------------------------------

thermobinary.criticalpoint(xcexpe,Tcexpe)
print('Critical conditions: xc=',thermobinary.xc,' at Tc=',thermobinary.Tc,'K')

# -----------------------------------------
# Determination of the monotectic condition
# -----------------------------------------

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
    thermobinary.Tmono=data['T'].values[0]
    x0=np.array([data['x1'].values[0],data['x2'].values[0]])
    thermobinary.x1mono,thermobinary.x2mono=thermobinary.binodale(x0,thermobinary.Tmono)
#endif

# ------------------------------------------------
# Determination of the binodal and spinodal curves
# ------------------------------------------------

NT=500
Nc=500
T=np.linspace(thermobinary.Tmono,thermobinary.Tc,NT)
x1=np.zeros(NT)
x2=np.zeros(NT)
x1spi=np.zeros(NT)
x2spi=np.zeros(NT)
x1[0]=thermobinary.x1mono
x2[0]=thermobinary.x2mono
x1[NT-1]=thermobinary.xc
x2[NT-1]=thermobinary.xc
x1spi[NT-1]=thermobinary.xc
x2spi[NT-1]=thermobinary.xc
x1spi[0],x2spi[0]=thermobinary.spinodal(x1[0],x2[0],T[0],Nc)
for i in range(1,NT-1):
    # Determination binodal curve
    x0=np.array([x1[i-1],x2[i-1]])
    x1[i],x2[i]=thermobinary.binodale(x0,T[i])
    
    # Determination of the spinodal curve
    x1spi[i],x2spi[i]=thermobinary.spinodal(x1[i],x2[i],T[i],Nc)
#end for

# --------
# Plotting
# --------

plt.figure()
plt.plot(x1expe,Texpe,'bo',label=thermobinary.ref)
plt.plot(x2expe,Texpe,'bo')
plt.plot(x1,T,'k-')
plt.plot(x2,T,'k-')
plt.plot(x1spi,T,'r-')
plt.plot(x2spi,T,'r-')
plt.annotate(str(np.round(thermobinary.xc,3))+','+str(np.round(thermobinary.Tc,2)),\
             (0.7*thermobinary.xc,1.001*thermobinary.Tc))
if (thermobinary.nature=='monotectic'):
    plt.plot([0.,thermobinary.x2mono],[thermobinary.Tmono,thermobinary.Tmono],'k-')
    #plt.plot([0,thermobinary.x1mono],[thermobinary.Tfsolvent,thermobinary.Tmono],'k-')
#endif
plt.xlabel(r'$c$',fontsize=12)
plt.ylabel(r'$T$ (K)',fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=10,loc='upper right')
plt.xlim((0,np.max([x2mono,thermobinary.x2mono])))
plt.savefig(thermobinary.system+'.pdf',dpi=300,bbox_inches='tight')
