#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 10:27:48 2025

This class is developed to determine thermodynamic functions and their derivatives
of binary system. Methods to find equilibrium state, critical condition, 
monotectic points and interfacial energy are provided.

The variables are:

    solvant: string corresponding to the solvant of the system, in general SiO2
    solute: string corresponding to the solute of the system.
    system: string corresponding to the binary system equal to 'solute-solvant'
    nature: Type of the miscibility gap, either 'monotectic' or 'subliquidus'
    ref: Reference of the article of experimental data used for comparison.
    Msolvant: Molar mass [kg/mol] of the solvant. 
    rhosolvant: density [kg/m^3] of the solvant. 
    Vmolsolvant: Molar volume [m^3/mol] of the solvant.
    DHfsolvant: Fusion molar enthalpy [J/mol] of the solvant
    Tfsolvant: Fusion temperature [K] of the pure solvant.
    Msolute: Molar mass [kg/mol] of the solute.
    rhosolute: density [kg/m^3] of the solute.
    Vmolsolute: Molar volume [m^3/mol] of the solute.
    DHfsolute: Fusion molar enthalpy [J/mol] of the solute.
    Tfsolute: Fusion temperature [K] of the pure solute.
    A: Float array of (2X2) corresponding to the coefficients of the Redlich-Kister potential.
    xc: Float corresponding to the solute molar fraction of the critical point
    Tc: Float corresponding to the critical temperature [K]
    x1mono: Float corresponding to the solute molar fraction of the left point of the monotectic line.
    x2mono: Float corresponding to the solute molar fraction of the right point of the monotectic line.
    Tmono: Float corresponding to the monotectic temperature.

@author: fpigeonneau
"""

import numpy as np
from scipy import optimize
from scipy.constants import N_A
from scipy.constants import R
from molarmass import MolarMass
import pandas as pd

class BinarySystem():
    
    def __init__(self,solvant='SiO2',solute=None):
        
        """
        Initialization of the parameters useful for the thermodynamic functions.
        
        The names of the solvant and the solute are given as arguments of the
        initialization of the class.
        
        The dboxides.csv file gathers the properties of the oxides. 
        The dbexcesspotential.csv file gather the data of the Redlich-Kister potential.
        
        
        """
        
        # ----------------------
        # Variables of the class
        # ----------------------
        
        self.solvant=None
        self.solute=None
        self.system=None
        self.nature=None
        self.ref=None
        self.Msolvant=None
        self.rhosolvant=None
        self.Vmolsolvant=None
        self.DHfsolvant=None
        self.Tfsolvant=None
        self.Msolute=None
        self.rhosolute=None
        self.Vmolsolute=None
        self.DHfsolute=None
        self.Tfsolute=None
        self.A=None
        self.xc=None
        self.Tc=None
        self.x1mono=None
        self.x2mono=None
        self.Tmono=None

        # Names of the solvant and solute
        # -------------------------------
        self.solvant=solvant
        self.solute=solute
        
        # Reading of the data-set of oxides
        # ---------------------------------
        dboxides=pd.read_csv('dboxides.csv',index_col=0)
        
        # Reading of the coefficients of the Redlich-Kister
        # -------------------------------------------------
        dbRKexcesspotential=pd.read_csv('dbexcesspotential.csv',index_col=0)
        
        # Loading of the properties of solvant
        # ------------------------------------
        isolvant=np.argwhere(dboxides.index==self.solvant)[0][0]
        self.Msolvant=MolarMass(self.solvant)
        self.rhosolvant=dboxides['rho'].values[isolvant]
        self.Vmolsolvant=self.Msolvant/self.rhosolvant
        self.DHfsolvant=dboxides['DHf'].values[isolvant]
        self.Tfsolvant=dboxides['Tf'].values[isolvant]
        
        # Loading of the properties of solute
        #------------------------------------
        isolute=np.argwhere(dboxides.index==self.solute)[0][0]
        self.Msolute=MolarMass(self.solute)
        self.rhosolute=dboxides['rho'].values[isolute]
        self.Vmolsolute=self.Msolute/self.rhosolute
        self.DHfsolute=dboxides['DHf'].values[isolute]
        self.Tfsolute=dboxides['Tf'].values[isolute]
        
        # Coefficient A for the Redlich-Kister potential
        # ----------------------------------------------
        
        self.system=solute+'-'+solvant
        isys=np.argwhere(dbRKexcesspotential.index==self.system)[0][0]
        self.nature=dbRKexcesspotential['Nature'].values[isys]
        self.A=np.array([[dbRKexcesspotential['A00'].values[isys],dbRKexcesspotential['A01'].values[isys]],
                    [dbRKexcesspotential['A10'].values[isys],dbRKexcesspotential['A11'].values[isys]]])
        self.ref=dbRKexcesspotential['Ref'].values[isys]
        
    #end __init__
    
    # ------------------------------
    # Gibbs energy of pure component
    # ------------------------------
    
    def DGf(self,T,DHf,Tf):
        """
        Determine the Gibbs energy as a function of T for a pure component.
        
        Parameters:
            T: Temperature (K)
            DHf: Fusion enthalpy (J/mol)
            Tf: Fusion temperature (K)
            
            Return DGf in (J/mol)
        """
    
        return DHf*(1.-T/Tf)
    # end DGf
    
    # ----------------------------------------
    # Gibbs energy of the mixing of the system
    # ----------------------------------------
    
    def freeenergy(self,c,T):
        """
        Determination of the Gibb energy of a binary system with an interaction
        potential using the Redlich-Kister potential at the first order.
        
        Parameters
        ----------
        c : Float
            Molar fraction of the solute.
        T : Float
            Temperature (K).
        Returns
        -------
        Float
            Gibbs energy of the mixing solution.
            
        """
        A0=self.A[0,0]+self.A[0,1]*T
        A1=self.A[1,0]+self.A[1,1]*T
        return (1.-c)*self.DGf(T,self.DHfsolvant,self.Tfsolvant)+\
                c*self.DGf(T,self.DHfsolute,self.Tfsolute)+\
                R*T*(c*np.log(c)+(1.-c)*np.log(1.-c))+c*(1.-c)*(A0+(1.-2.*c)*A1)
    # end freeenergy
    
    def d2freeenergydc2(self,c,T):
        """
        Determination of the seconde derivative of the Gibb energy of a binary system 
        with an excess potential using the Redlich-Kister potential at the 
        first order.
        
        Parameters
        ----------
        c : Float
            Molar fraction of the solute.
            T : Float
            Temperature (K).
            
        Returns
        -------
            Float
                Second derivative respect to c of the Gibbs energy of the mixing solution.
            
        """
        A0=self.A[0,0]+self.A[0,1]*T
        A1=self.A[1,0]+self.A[1,1]*T
        return R*T/(c*(1.-c))-2*A0-6.*A1*(1.-2.*c)
    # end d2freeenergydc2
    
    def musolvant(self,T,c):
        """
        Chemical potential of the solvant in the binary system with an interaction
        potential using the Redlich-Kister potential at the first order.
        
        Parameters
        ----------
        T : Float
            Temperature (K).
            c : Float
            Molar fraction of the solvant.
        
        Returns
        -------
        Float
            chemical potential of the solvant.
            
        """
        A0=self.A[0,0]+self.A[0,1]*T
        A1=self.A[1,0]+self.A[1,1]*T
        return self.DGf(T,self.DHfsolvant,self.Tfsolvant)+R*T*np.log(1.-c)+c**2*(A0+(3.-4.*c)*A1)
    # end self.musolvant
    
    def musolute(self,T,c):
        """
            Chemical potential of the solute in the binary system with an interaction
            potential using the Redlich-Kister potential at the first order.
        
        Parameters
        ----------
            T : Float
            Temperature (K).
            c : Float
            Molar fraction of the solute.
        
        Returns
        -------
            Float
            chemical potential of the solute.
        
        """
        A0=self.A[0,0]+self.A[0,1]*T
        A1=self.A[1,0]+self.A[1,1]*T
        return self.DGf(T,self.DHfsolute,self.Tfsolute)+R*T*np.log(c)+(1.-c)**2*(A0+(1.-4.*c)*A1)
    # end musolute
    
    def dmusolvantdc(self,T,c):
        """
        First derivative respect to c of chemical potential of the solvant in 
        the binary system with an interaction potential using the Redlich-Kister 
        potential at the first order.
        
        Parameters
        ----------
        T : Float
            Temperature (K).
            c : Float
            Molar fraction of the solvant.
        
        Returns
        -------
        Float
            First derivative of chemical potential of the solvant respect to C.
            
        """
        A0=self.A[0,0]+self.A[0,1]*T
        A1=self.A[1,0]+self.A[1,1]*T
        return -R*T/(1.-c)+2.*c*A0+6.*c*(1.-2.*c)*A1
    # end dmusolvantdc
    
    def dmusolutedc(self,T,c):
        """
            First derivative of chemical potential of the solute respect to c 
            in the binary system with an interaction potential using the Redlich-
            Kister potential at the first order.
        
        Parameters
        ----------
            T : Float
            Temperature (K).
            c : Float
            Molar fraction of the solute.
        
        Returns
        -------
            Float
            First derivative of chemical potential of the solute respect to c.
        
        """
        A0=self.A[0,0]+self.A[0,1]*T
        A1=self.A[1,0]+self.A[1,1]*T
        return R*T/c-2.*(1.-c)*A0-6.*(1.-c)*(1.-2.*c)*A1
    # end dmusolutedc
    
    def fzerocritpoint(self,x):
        """
        Function determining of the zeros of relations defining the critical 
        point.
        
        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
        
        Returns
        -------
        y : Array of float (1x2)
            Zeros of the critical point definition.
        
        """
        xc=x[0]
        Tc=x[1]
        A0=self.A[0,0]+self.A[0,1]*Tc
        A1=self.A[1,0]+self.A[1,1]*Tc
        
        y=np.zeros(2)
        y[0]=R*Tc*(1./xc+1./(1.-xc))-2.*A0-6.*A1*(1.-2.*xc)
        y[1]=R*Tc*(1./(1.-xc)**2-1./xc**2)+12.*A1
        
        return y
    # end fzerocritpoint
    
    def criticalpoint(self,xc0,Tc0):
        """
        
        Determination of the critical condition.
        
        Parameters
        ----------
        xc0 : Float
            Initial value the molar fraction of the solute in critical condition.
        Tc0 : Float
            Initial value of the critical temperature.
        
        Returns
        -------
        None.
        """
        
        x=optimize.fsolve(self.fzerocritpoint,np.array([xc0,Tc0]))
        self.xc=x[0]
        self.Tc=x[1]
    #end criticalpoint
    
    def fmonotectic(self,x):
        """
        Function determining of the zeros of relations definiting the monotetic state
        
        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
        
        Returns
        -------
        y : Array of float (1x2)
            Zeros of the critical point definition.
        """
        x1=x[0]
        x2=x[1]
        Tm=x[2]
        y=np.zeros(3)
        y[0]=self.musolvant(Tm,x1)
        y[1]=self.musolvant(Tm,x2)
        y[2]=self.musolute(Tm,x1)-self.musolute(Tm,x2)
        return y
    #end fmonotectic
    
    
    def monotectic(self,x1mono,x2mono,Tmono):
        """
        

        Parameters
        ----------
        x1mono : Float
            Guest value of the left point of the monotectic conditions.
        x2mono : Float
            Guest value of the right point of the monotectic conditions.
        Tmono : Float
            Guest value of the monotectic temperature.

        Returns
        -------
        None.

        """
        x0=np.array([x1mono,x2mono,Tmono])
        x=optimize.fsolve(self.fmonotectic,x0)
        self.x1mono=x[0]
        self.x2mono=x[1]
        self.Tmono=x[2]
    #end monotectic
    
    def fgap(self,x,T):
        """
        Function determining the balance between the chemical potentials of the 
        equilibrium between the two phases.
        
        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
        T : Float
            Temperature (K).
        
        Returns
        -------
        y : Array of floats (1x2)
            zeros of the relations to satisfy the thermodynamic equilibrium 
            between the two phases.
        
        """
        y=np.zeros(2)
        
        y[0]=self.musolvant(T,x[0])-self.musolvant(T,x[1])
        y[1]=self.musolute(T,x[0])-self.musolute(T,x[1])
        return y
    #end fgap
    
    def binodale(self,x0,T):
        """
        Determination of the molar fractions corresponding to the equilibrium
        between the two phases on the binodal conditions.

        Parameters
        ----------
        x0 : Array of float of size 2.
            guest values of the equilibrium molar fractions.
        T : Float
            Temperature [K].

        Returns
        -------
        Float
            Molar fraction of the left point.
        Float
            Molar fraction of the right point.

        """
        
        x=optimize.fsolve(self.fgap,x0,args=(T,))
        
        # Return of the molar fractions of the equilibrium conditions:
        return x[0],x[1]
    #end binodale
    
    def fsurftenion(self,x,xalpha,T):
        """
        Determination of the difference of the chemical potentials in the interface
        and at the interface corresponding to the interfacial energy.
        
        Parameters
        ----------
        x : Float
            molar fraction of the solute in the interface.
        xalpha : Float
            molar fraction of the solute in the thermodynamic equilibrium state
            between the two states in the continuous phase.
            T : Float
                Temperature (K).
        
        Returns
        -------
        y : Float
            Difference between the two definitions of the surface tension.
        
        """
        y=(self.musolvant(T,x)-self.musolvant(T,xalpha))/self.Vmolsolvant**(2./3.)-\
          (self.musolute(T,x)-self.musolute(T,xalpha))/self.Vmolsolute**(2./3.)
        
        return y
    #end fsurftenion
    
    def surfacetension(self,xalpha,xbeta,T):
        """
        Determination of the surface tension of the binary system according to 
        Kaptay's model.
        
        Parameters
        ----------
        xalpha : Float
            molar fraction of the solute in the thermodynamic equilibrium state
            between the two states in the continuous phase.
        xbeta : Float
            molar fraction of the solute in the thermodynamic equilibrium state
            between the two states in the dispersed phase.
        T : Float
            Temperature (K).
        
        Returns
        -------
        xi : Float
            Molar fraction of the solute in the interface.
        gamma : Float
            Interfacial tension between the two phases (N/m).
        
        """
        # Determination of the average
        x0=0.5*(xalpha+xbeta)
        xi=optimize.fsolve(self.fsurftenion,x0,args=(xalpha,T))
        
        # Computation of the surface tension
        c=(np.pi/6.)**(1./3.)
        gamma=(self.musolvant(T,xi[0])-self.musolvant(T,xalpha))/(c*self.Vmolsolvant**(2./3.)*N_A**(1./3.))
        
        # Return of results
        return xi[0],gamma
    #end surfacetension
    
    
    def spinodal(self,x1,x2,T,N):
        """
        Determination of the curve of the spinodal condition.
        
        Parameters
        ----------
        x1 : Float
            Binodal location less than the critical concentration.
        x2 : Float
            Binodal location larger than the critical concentration.
        T : Float
            Temperature (K).
        N : Integer
            Numbers of nodes in linear space in [x1,x2]
        
        Returns
        -------
        x1spi : Float
            Spinodal location less than the critical concentration.
        x2spi : Float
            Spinodal location larger than the critical concentration..
        
        """
        # Range of c between the binonal points
        c=np.linspace(x1,x2,N)
        
        # Computation of the second derivative of the Gibbs energy
        d2fdc2=self.d2freeenergydc2(c,T)
        
        # Research of the initial values of c where the second derivative changes of sign
        argd2fdc2=[]
        for j in range(N-1):
            if (np.sign(d2fdc2[j+1])!=np.sign(d2fdc2[j])):
                argd2fdc2=np.append(argd2fdc2,j)
            #end if
        #end for
        argd2fdc2=np.int64(argd2fdc2)
        x1spi=c[argd2fdc2[0]]
        x2spi=c[argd2fdc2[1]]
        
        # Determination of accurate values of cancellation of the second derivative
        x=optimize.fsolve(self.d2freeenergydc2,x1spi,args=(T,))
        x1spi=x[0]
        x=optimize.fsolve(self.d2freeenergydc2,x2spi,args=(T,))
        x2spi=x[0]
        
        # Reurn to the solutions
        return x1spi,x2spi
    # end spinodal

    def MolarVolume(self,x):
        return self.Vmolsolvant*(1.-x)+self.Vmolsolute*x
    # end MolarVolume
    
    def criticalradius(self,x0,x2,gamma,T):
        """
        Determination of the critical radius.
        
        Parameters
        ----------
        x0 : Float
            Molar fraction of the homogeneous liquid.
        x2 : Float
            Molar fraction of the dispersed phase beta.
        T : Float
            Temperature (K).
        
        Returns
        -------
        rcrit : Float
            Critical radius of the nucleii.
        DGcrit : Float
            Critical Gibbs energy
            
        """
        
        DeltamuA=self.musolvant(T,x0)-self.musolvant(T,x2)
        DeltamuB=self.musolute(T,x0)-self.musolute(T,x2)
        Deltatot=(DeltamuA*(1.-x2)+DeltamuB*x2)/self.MolarVolume(x2)
        if (Deltatot!=0.):
            rcrit=2.*gamma/Deltatot
            DGcrit=16.*np.pi*gamma**3/(3.*Deltatot**2)
        else:
            rcrit=0.
            DGcrit=0.
        #endif
        return rcrit,DGcrit,Deltatot
    # end criticalradius
    
    def deltag(self,x,T,xalpha):
        """
        Determination of the Gibbs energy minus the equilibrium state line

        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
        T : Float
            Temperature (K).
        xalpha : Float
            Molar fraction in the phase alpha at the thermodynamic equilibrium.
        
        Returns
        -------
        Float
            The Gibbs energy minus the equilibrium state line.
        
        """
        return (1.-x)*(self.musolvant(T,x)-self.musolvant(T,xalpha))+\
                x*(self.musolute(T,x)-self.musolute(T,xalpha))
    # end deltag
    
    def ddeltagdc(self,x,T,xalpha):
        """
        Determination of the first derivative of the Gibbs energy minus the 
        equilibrium state line respect to x.

        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
        T : Float
            Temperature (K).
        xalpha : Float
            Molar fraction in the phase alpha at the thermodynamic equilibrium.
        
        Returns
        -------
        Float
            The first derivative of the Gibbs energy minus the equilibrium 
            state line.
        
        """
        return self.musolute(T,x)-self.musolute(T,xalpha)-\
                (self.musolvant(T,x)-self.musolvant(T,xalpha))
    # end ddeltagdc
    
    def ddeltamudc(self,x,x0,T):
        """
        First derivative of the deltamu respect to x.

        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
        x0 : Float
            Molar fraction of the solute in the homogeneous liquid.
        T : Float
            Temperature of the system (K).

        Returns
        -------
        First derivative of deltamu.

        """
        
        Vmol=self.MolarVolume(x)
        dVmoldc=self.Vmolsolute-self.Vmolsolvant
        Deltamu=self.deltag(x,T,x0)/Vmol
        return self.ddeltagdc(x,T,x0)/Vmol-Deltamu*dVmoldc/Vmol
    # end ddeltamudc
    
    def d2deltamudc2(self,x,x0,T):
        """
        Second derivative of the deltamu respect to x.

        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
            x0 : Float
            Molar fraction of the solute in the homogeneous liquid.
        T : Float
        Temperature of the system (K).

        Returns
        -------
        d2dmudc2 : Float
        Second derivative of deltamu.

        """
        
        Vmol=self.MolarVolume(x)
        dVmoldc=self.Vmolsolute-self.Vmolsolvant
        Deltamu=self.deltag(x,T,x0)/Vmol
        d2dmudc2=(self.dmusolutedc(T,x)-self.dmusolvantdc(T,x))/Vmol-\
                 (self.ddeltagdc(x,T,x0)/Vmol+self.ddeltamudc(x,x0,T))*dVmoldc/Vmol+\
                 Deltamu*(dVmoldc/Vmol)**2
        return d2dmudc2
    # end ddeltamudc


    def Integrand(self,x,T,xalpha):
        """
        Determination of the integrand of the computation of the interfacial
        energy according to the Cahn-Hilliard model.

        Parameters
        ----------
        x : Float
            Molar fraction.
        T : Float
            Temperature [K].
        xalpha : Float
            Molar fraction of the equilibrium condition.

        Returns
        -------
        Float
            Integrand of the integral to determine the interfacial energy.
        """
        return np.sqrt(self.deltag(x,T,xalpha))
    #end Integrand
    
#end class