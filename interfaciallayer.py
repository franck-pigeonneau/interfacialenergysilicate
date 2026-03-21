#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 11:23:51 2022

@author: fpigeonneau
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import zero_Celsius
from scipy import optimize
from scipy import integrate
import pandas as pd
from binarysystem import BinarySystem

#import sys
#eps=sys.float_info.epsilon

def surfacetension(theta,A,n):
    return A*(1.-theta)**n
#end surfacetension


# ---------------------------------
# Loading of the database of oxides
# ---------------------------------

solvant='SiO2'
solute='BaO'
thermobinary=BinarySystem(solvant,solute)

# Données sur la tension de surface
dbgamma=pd.read_csv('gammabynarysystem.csv',index_col=0)
isys=np.argwhere(dbgamma.index==thermobinary.system)[0][0]
gamma0=dbgamma['A'].values[isys]
n=dbgamma['n'].values[isys]

# Determination des conditions critiques
# --------------------------------------

dbcrit=pd.read_csv('critcondbinarysystem.csv',index_col=0)
isys=np.argwhere(dbcrit.index==thermobinary.system)[0][0]
xc=dbcrit['xc'].values[isys]
Tc=dbcrit['Tc'].values[isys]
thermobinary.criticalpoint(xc,Tc)
print('Critical conditions: xc=',thermobinary.xc,' at Tc=',thermobinary.Tc-zero_Celsius,'°C')

# Determination de la ligne monotectique
# --------------------------------------

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

# Température de Tc à Tm
NT=701
Nc=200
T=np.linspace(thermobinary.Tmono,0.99*thermobinary.Tc,NT)
zeta=np.zeros(NT)
x1=np.zeros(NT)
x2=np.zeros(NT)
x1[0]=thermobinary.x1mono
x2[0]=thermobinary.x2mono
x0=np.array([x1[0],x2[0]])

for i in range(NT):
    # Calcul des concentrations d'équilibre
    x1[i],x2[i]=thermobinary.binodale(x0,T[i])
    
    # Calcul de l'integrale
    I=integrate.quad(thermobinary.Integrand,x1[i],x2[i],args=(T[i],x1[i]))[0]
    
    # Recherche du maximum locale
    c=np.linspace(x1[i],x2[i],Nc)
    ddgdc=thermobinary.ddeltagdc(c,T[i],x1[i])
    for j in range(1,Nc-2):
        if (np.sign(ddgdc[j+1])!=np.sign(ddgdc[j])):
            argddgdc=j
        #end if
    #end for
    xdgmax=optimize.fsolve(thermobinary.ddeltagdc,c[argddgdc],args=(T[i],x1[i]))
    dgmax=thermobinary.deltag(xdgmax[0],T[i],x1[i])
    
    # Calcul de l'echelle de l'interface diffuse
    zeta[i]=thermobinary.Vmolsolvant*surfacetension(T[i]/thermobinary.Tc,gamma0,n)/(np.sqrt(2.*dgmax)*I)
    
    # New guest of x0
    x0=np.array([x1[i],x2[i]])
# end for

# ---------------------
# Non-linear regression
# ---------------------

popt,pcov=optimize.curve_fit(surfacetension,T[np.int64(NT//1.3):]/thermobinary.Tc,zeta[np.int64(NT//1.3):])

plt.figure()
plt.plot(T,zeta*1e9,'k',label='Eq. (26)',linewidth=2)
plt.plot(T,surfacetension(T/thermobinary.Tc,*popt)*1e9,'r-.',label=r'$\zeta='+\
         str(np.round(popt[0]*1e9,2))+'/(1-T/T_c)^{'+str(np.round(-popt[1],2))+'}$',linewidth=2)
plt.xlabel(r'$T$ (K)',fontsize=14)
plt.ylabel(r'$\zeta$ (nm)',fontsize=14)
plt.legend(loc=0,fontsize=12)
plt.savefig('zetavsT-'+thermobinary.system+'.pdf',dpi=300,bbox_inches='tight')

# Saving of the results
# ---------------------
A=np.transpose(np.array([T,zeta]))
data=pd.DataFrame(A,columns=np.array(['T','zeta']))
data.to_csv('zetavsT-'+thermobinary.system+'.csv')

plt.show()

