#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 10:27:48 2025

This class is developed to determine thermodynamic functions and their derivatives
of binary system. Methods to find equilibrium state, critical condition, 
monotectic points and interfacial energy are provided.

The variables are:

    solvent: string corresponding to the solvent of the system, in general SiO2 or B2O3.
    solute: string corresponding to the solute of the system.
    OrderRK: Order of the Redlich-Kister potential.
    NTRK: Number of coefficients of function of T An.
    system: string corresponding to the binary system equal to 'solute-solvent'
    nature: Type of the miscibility gap, either 'monotectic' or 'subliquidus'
    ref: Reference of the article of experimental data used for comparison.
    Msolvent: Molar mass [kg/mol] of the solvent. 
    rhosolvent: density [kg/m^3] of the solvent. 
    Vmolsolvent: Molar volume [m^3/mol] of the solvent.
    DHfsolvent: Fusion molar enthalpy [J/mol] of the solvent
    Tfsolvent: Fusion temperature [K] of the pure solvent.
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

For more details see [1].

References
----------

F. Pigeonneau & W. Blanc (2026). Interfacial energy in phase separated borate 
and silicate systems J. Am. Ceram. Soc., 10.1111/JACE.71109.

@author: fpigeonneau
"""

import numpy as np
from scipy import optimize
from scipy.constants import N_A
from scipy.constants import R
from molarmass import MolarMass
import pandas as pd

class BinarySystem():
    
    def __init__(self,solvent='SiO2',solute=None,OrderRK=1,NTRK=2):
        
        """
        Initialization of the parameters useful for the thermodynamic functions.
        
        The names of the solvent and the solute are given as arguments of the
        initialization of the class. The order of the Redlich-Kister potential is
        also given as the thrird argument as well as the number of coefficient of
        the temperature-dependence function.
        
        The dboxides.csv file gathers the properties of the oxides. 
        The dbexcesspotential.csv file gather the data of the Redlich-Kister potential.
        
        """
        
        # ----------------------
        # Variables of the class
        # ----------------------
        
        self.solvent=None
        self.solute=None
        self.OrderRK=None
        self.NTRK=None
        self.alphaRK=None
        self.system=None
        self.nature=None
        self.ref=None
        self.Msolvent=None
        self.rhosolvent=None
        self.Vmolsolvent=None
        self.DHfsolvent=None
        self.Tfsolvent=None
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
        
        # Names of the solvent and solute
        # -------------------------------
        self.solvent=solvent
        self.solute=solute
        
        # Type of the Redlich-Kister potential
        # ------------------------------------
        
        self.OrderRK=OrderRK
        self.NTRK=NTRK
        
        # Definition of the exponents of the thermal coefficients in the 
        # Redlich-Kister potential
        self.alphaRK=np.zeros(self.NTRK)
        for i in range(self.NTRK-1):
            self.alphaRK[i]=i
        #end for
        self.alphaRK[self.NTRK-1]=-0.5
        
        # Reading of the data-set of oxides
        # ---------------------------------
        dboxides=pd.read_csv('dboxides.csv',index_col=0)
        
        # Reading of the coefficients of the Redlich-Kister
        # -------------------------------------------------
        dbRKexcesspotential=pd.read_csv('dbexcesspotential.csv',index_col=0)
        
        # Loading of the properties of solvent
        # ------------------------------------
        isolvent=np.argwhere(dboxides.index==self.solvent)[0][0]
        self.Msolvent=MolarMass(self.solvent)
        self.rhosolvent=dboxides['rho'].values[isolvent]
        self.Vmolsolvent=self.Msolvent/self.rhosolvent
        self.DHfsolvent=dboxides['DHf'].values[isolvent]
        self.Tfsolvent=dboxides['Tf'].values[isolvent]
        
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
        
        self.system=solute+'-'+solvent
        isys=np.argwhere(dbRKexcesspotential.index==self.system)[0][0]
        self.nature=dbRKexcesspotential['Nature'].values[isys]
        self.A=np.zeros((self.OrderRK+1,self.NTRK))
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                self.A[i,j]=dbRKexcesspotential['A'+str(i)+str(j)].values[isys]
            #end for
        #end for
        
        # Reference used for validation on the phase diagram
        # --------------------------------------------------
        
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
            
        Return:
            DGf in (J/mol)
        """
    
        return DHf*(1.-T/Tf)
    # end DGf
    
    def RedlichKister(self,T):
        """
        Determination of the coefficients of the Redlich-Kister potential as a
        function of the temperature.

        Parameters
        ----------
        T : Float
            Temperature in K.

        Returns
        -------
        Array od float
            Coefficient of the Redlich-Kister potential.
        """
        
        ARK=np.zeros(self.OrderRK+1)
        for i in range(self.OrderRK+1):
            ARK[i]=0.
            for j in range(self.NTRK):
                ARK[i]+=self.A[i,j]*T**self.alphaRK[j]
            #end for
        #end for
        
        # Return of the coefficient
        return ARK
    #end RedlichKister
    
    def RedlichP(self,n,x):
        """
        Function of the Redlich-Kister at the order n.

        Parameters
        ----------
        n : Float
            Exponent of the function (1-2x).
        x : Float
            Molar concentrzation of the solute.

        Returns
        -------
        Float
            Function of the Redlich-Kister at the order n.
        """
        return (1.-2.*x)**n
    #end RedlichP
    
    def dmRedlichPdxm(self,n,m,x):
        """
        Derivative as a respect of x at the order m of the Redlich-Kister at the
        order n.
        
        Parameters
        ----------
        n: Float
            Exponent of the function (1-2x).
        m: Float
            order of the derivative.
        x : Float
            Molar concentrzation of the solute.
        
        Returns
        -------
        Derivative as a respect of x at the order m.
        """
        if (n<m):
            return 0.
        else:
            c=1.
            for i in range(m):
                c*=-2.*(n-i)
            #end for
            return c*(1.-2.*x)**(n-m)
        #endif
    #end dmRedlichPdxm
    
    def Omega(self,c,T):
        """
        Determination of the excess potential.
        
        Parameters
        ----------
        c: Float
            Molar concentration of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Excess potential from the Redlich-Kister function.
        """
        # Determination of the coefficients of the excess potential
        # ---------------------------------------------------------
        
        ARK=self.RedlichKister(T)
        
        # Compute of the excess potential
        # -------------------------------
        
        w=0.
        for i in range(self.OrderRK+1):
            w+=ARK[i]*self.RedlichP(i,c)
        #end for
        
        return w
    #end Omega

    def dOmegadc(self,c,T):
        """
        First derivative of the Redlich-Kister potential.
        
        Parameters
        ----------
        c: Float
            Molar concentration of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        First derivative of the Redlich-Kister excess potential.
        
        """
        # Determination of the coefficients of the Redlich-Kister potential
        # -----------------------------------------------------------------
        
        ARK=self.RedlichKister(T)
        
        # Compute of the first derivative of the Redlich-Kister potential
        # ---------------------------------------------------------------
        dwdc=0.
        for i in range(self.OrderRK+1):
            dwdc+=ARK[i]*self.dmRedlichPdxm(i,1,c)
        #end for
        
        return dwdc
    #end dOmegadc
    
    def d2Omegadc2(self,c,T):
        """
        Second derivative of the Redlich-Kister potential.
       
        Parameters
        ----------
        c: Float
            Molar concentration of the solute.
        T: Float
           Temperature in K.
       
        Returns
        -------
        Second derivative of the Redlich-Kister excess potential.
       
        """
      
        # Determination of the coefficients of the Redlich-Kister potential
        # -----------------------------------------------------------------
        
        ARK=self.RedlichKister(T)
        
        # Compute of the first derivative of the Redlich-Kister potential
        # ---------------------------------------------------------------
        d2wdc2=0.
        for i in range(self.OrderRK+1):
            d2wdc2+=ARK[i]*self.dmRedlichPdxm(i,2,c)
        #end for
        
        return d2wdc2
    #end d2Omegadc2
    
    def d3Omegadc3(self,c,T):
        """
        Second derivative of the Redlich-Kister potential.
      
        Parameters
        ----------
        c: Float
           Molar concentration of the solute.
        T: Float
           Temperature in K.
      
        Returns
        -------
        Second derivative of the Redlich-Kister excess potential.
        
        """
       
        # Determination of the coefficients of the Redlich-Kister potential
        # -----------------------------------------------------------------
        
        ARK=self.RedlichKister(T)
        
        # Compute of the first derivative of the Redlich-Kister potential
        # ---------------------------------------------------------------
        d3wdc3=0.
        for i in range(self.OrderRK+1):
            d3wdc3+=ARK[i]*self.dmRedlichPdxm(i,3,c)
        #end for
        
        return d3wdc3
    #end d3Omegadc3
    
    
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
        
        # Return of the Gibbs free energy
        # -------------------------------
        
        return (1.-c)*self.DGf(T,self.DHfsolvent,self.Tfsolvent)+\
                c*self.DGf(T,self.DHfsolute,self.Tfsolute)+\
                R*T*(c*np.log(c)+(1.-c)*np.log(1.-c))+\
                c*(1.-c)*self.Omega(c,T)
    # end freeenergy
    
    def d2freeenergydc2(self,c,T):
        """
        Determination of the seconde derivative of the Gibb energy of a binary system 
        with an excess potential using the Redlich-Kister potential.
        
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
        
        # Compute and return of the second derivative respect to c of the Gibbs free energy
        return R*T/(c*(1.-c))-2*self.Omega(c,T)+2.*(1.-2.*c)*self.dOmegadc(c,T)+\
               c*(1.-c)*self.d2Omegadc2(c,T)
    # end d2freeenergydc2
    
    def d3freeenergydc3(self,c,T):
        """
        Determination of the third derivative of the Gibb energy of a binary system 
        with an excess potential using the Redlich-Kister potential.
        
        Parameters
        ----------
        c : Float
            Molar fraction of the solute.
            T : Float
            Temperature (K).
            
        Returns
        -------
            Float
                Third derivative respect to c of the Gibbs energy of the mixing solution.
            
        """
        
        # Compute and return of the second derivative respect to c of the Gibbs free energy
        return R*T*(1./(1.-c)**2-1./c**2)-6*self.dOmegadc(c,T)+3.*(1.-2.*c)*self.d2Omegadc2(c,T)+\
               c*(1.-c)*self.d3Omegadc3(c,T)
    # end d3freeenergydc3
    
    def musolvent(self,T,c):
        """
        Chemical potential of the solvent in the binary system with an interaction
        potential using the Redlich-Kister potential at the first order.
        
        Parameters
        ----------
        T : Float
            Temperature (K).
            c : Float
            Molar fraction of the solvent.
        
        Returns
        -------
        Float
            chemical potential of the solvent.
            
        """
        
        # Compute and return of the chemical potential of the solvent
        # -----------------------------------------------------------
        
        return self.DGf(T,self.DHfsolvent,self.Tfsolvent)+R*T*np.log(1.-c)+\
               c**2*(self.Omega(c,T)-(1.-c)*self.dOmegadc(c,T))
    # end self.musolvent
    
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
                
        return self.DGf(T,self.DHfsolute,self.Tfsolute)+R*T*np.log(c)+\
               (1.-c)**2*(self.Omega(c,T)+c*self.dOmegadc(c,T))
    # end musolute
    
    def dmusolventdc(self,T,c):
        """
        First derivative respect to c of chemical potential of the solvent in 
        the binary system with an interaction potential using the Redlich-Kister 
        potential.
        
        Parameters
        ----------
        T : Float
            Temperature (K).
            c : Float
            Molar fraction of the solvent.
        
        Returns
        -------
        Float
            First derivative of chemical potential of the solvent respect to C.
            
        """
        
        # Compute and return of the first derivate respect to c of the chemical
        # potential of the solvent
        
        return -R*T/(1.-c)+2.*c*(self.Omega(c,T)-(1.-2.*c)*self.dOmegadc(c,T)-\
                                 0.5*c*(1.-c)*self.d2Omegadc2(c,T))
    # end dmusolventdc
    
    def dmusolutedc(self,T,c):
        """
            First derivative of chemical potential of the solute respect to c 
            in the binary system with an interaction potential using the Redlich-
            Kister potential.
        
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
        
        # Compute and return of the first derivate respect to c of the chemical
        # potential of the solute
        
        return R*T/c-2.*(1.-c)*(self.Omega(c,T)-(1.-2.*c)*self.dOmegadc(c,T)-\
                                 0.5*c*(1.-c)*self.d2Omegadc2(c,T))
    # end dmusolutedc
    
    def fzerocritpoint(self,x):
        """
        Function determining of the zeros of relations defining the critical 
        point.
        
        Parameters
        ----------
        x : Array of 2 floats
            x[0]: molar fraction in the critical condition
            x[1]: temperature in K in the critical condition 
        Returns
        -------
        y : Array of float (1x2)
            Zeros of the critical point definition.
            y[0]: Second derivatives of the Gibbs free energy
            y[1]: Third derivative of the Gibbs free energy 
        """
        xc=x[0]
        Tc=x[1]
                
        # Determination of the second and the third derivatives of the Gibbs
        # free energy to cancel in the critical condition
        
        y=np.zeros(2)
        y[0]=self.d2freeenergydc2(xc,Tc)
        y[1]=self.d3freeenergydc3(xc,Tc)
        
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
        x : Array of 3 floats
            x[0]: molar fraction in the left part of the binodal limit
            x[1]: molar fraction in the right part of the binodal limit
            x[1]: temperature in K in the monotectic condition 
        
        Returns
        -------
        y : Array of float (1x3)
            Zeros of the monotectic points definition.
        """

        x1=x[0]
        x2=x[1]
        Tm=x[2]
        y=np.zeros(3)
        y[0]=self.musolvent(Tm,x1)
        y[1]=self.musolvent(Tm,x2)
        y[2]=self.musolute(Tm,x1)-self.musolute(Tm,x2)
        return y
    #end fmonotectic
    
    
    def monotectic(self,x1mono,x2mono,Tmono):
        """
        Determination of the monotectic conditions

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
        
        y[0]=self.musolvent(T,x[0])-self.musolvent(T,x[1])
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
    
    def fsurftension(self,x,xalpha,T):
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
        y=(self.musolvent(T,x)-self.musolvent(T,xalpha))/self.Vmolsolvent**(2./3.)-\
          (self.musolute(T,x)-self.musolute(T,xalpha))/self.Vmolsolute**(2./3.)
        
        return y
    #end fsurftension
    
    def gammapureoxide(self,T,gamma0,Agamma,Tgamma):
        """
        Determination of the surface tension of a pure oxide function of T.
        
        Parameters
        ----------
        T : Float
            Temperature in [K]
    
        gamma0 : Float
            Prefactor of the surface tension en [N/m]
            
        Agamma : Float
            Thermal coefficient of the variation of the surface tension [N/(mK)]
        Tgamma : Float
            Temperature of the reference in which gamma0 is determined [K].
            
        Returns
        
        Float corresponding to the surface tension of a pure oxide [N/m]
        """
        
        return gamma0+Agamma*(T-Tgamma)
    #end gammapureoxide
    
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
        xi=optimize.fsolve(self.fsurftension,x0,args=(xalpha,T))
        
        # Computation of the surface tension
        c=(np.pi/6.)**(1./3.)
        gamma=2.*(self.musolvent(T,xi[0])-self.musolvent(T,xalpha))/(c*self.Vmolsolvent**(2./3.)*N_A**(1./3.))
        
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
        return self.Vmolsolvent*(1.-x)+self.Vmolsolute*x
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
        
        DeltamuA=self.musolvent(T,x0)-self.musolvent(T,x2)
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
        return (1.-x)*(self.musolvent(T,x)-self.musolvent(T,xalpha))+\
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
                (self.musolvent(T,x)-self.musolvent(T,xalpha))
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
        dVmoldc=self.Vmolsolute-self.Vmolsolvent
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
        dVmoldc=self.Vmolsolute-self.Vmolsolvent
        Deltamu=self.deltag(x,T,x0)/Vmol
        d2dmudc2=(self.dmusolutedc(T,x)-self.dmusolventdc(T,x))/Vmol-\
                 (self.ddeltagdc(x,T,x0)/Vmol+self.ddeltamudc(x,x0,T))*dVmoldc/Vmol+\
                 Deltamu*(dVmoldc/Vmol)**2
        return d2dmudc2
    # end ddeltamudc
    
    def activity(self,x,T):
        """
        Determination of the activites of the solvent and the solute.

        Parameters
        ----------
        x : Float
            Molar fraction of the solute.
        T : Float
            Temperature in Kelvin.

        Returns
        -------
        asolvent : Float
            Activity of the solvent.
        asolute : Float
            Activity of the solute.
        """
        asolvent=(1.-x)*np.exp(x**2*(self.Omega(x,T)-\
                                     (1.-x)*self.dOmegadc(x,T))/(R*T))
        asolute=x*np.exp((1.-x)**2*(self.Omega(x,T)+x*self.dOmegadc(x,T))/(R*T))
        return asolvent,asolute
    #end activity
    
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
    
    def JacOmega(self,c,T):
        """
        Determination of the derivative respect of the coefficients of the
        Redlich-Kister potential.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of the Redlich-Kister potential
        """
        jOmega=np.zeros((self.OrderRK+1,self.NTRK))
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                jOmega[i,j]=T**self.alphaRK[j]*(1.-2.*c)**i
            #end for
        #end for
        return jOmega
    #end JacOmega
    
    def JacdOmegadc(self,c,T):
        """
        Determination of the derivative respect of the coefficients of the
        first derivative as a respect of c of the Redlich-Kister potential.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of the first derivative respect of c of the Redlich-Kister 
            potential.
        """
        jdOmegadc=np.zeros((self.OrderRK+1,self.NTRK))
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                jdOmegadc[i,j]=-2.*i*T**self.alphaRK[j]*(1.-2.*c)**(i-1)
            #end for
        #end for
        return jdOmegadc
    #end JacdOmegadc
    
    def Jacd2Omegadc2(self,c,T):
        """
        Determination of the derivative respect of the coefficients of the
        second derivative as a respect of c of the Redlich-Kister potential.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of the second derivative respect of c of the Redlich-Kister 
            potential.
        """
        jd2Omegadc2=np.zeros((self.OrderRK+1,self.NTRK))
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                jd2Omegadc2[i,j]=4.*i*(i-1.)*T**self.alphaRK[j]*(1.-2.*c)**(i-2)
            #end for
        #end for
        return jd2Omegadc2
    #end Jacd2Omegadc2
    
    def Jacd3Omegadc3(self,c,T):
        """
        Determination of the derivative respect of the coefficients of the
        third derivative as a respect of c of the Redlich-Kister potential.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of the third derivative respect of c of the Redlich-Kister 
            potential.
        """
        jd3Omegadc3=np.zeros((self.OrderRK+1,self.NTRK))
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                jd3Omegadc3[i,j]=-8.*i*(i-1.)*(i-2.)*T**self.alphaRK[j]*(1.-2.*c)**(i-3)
            #end for
        #end for
        return jd3Omegadc3
    #end Jacd3Omegadc3
    
    def Jacmusolvent(self,c,T):
        """
        Determination of the Jacobian of the chemical potential of the solvent.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of the chemical potential of the solvent.
        """
        return c**2*(self.JacOmega(c,T)-(1.-c)*self.JacdOmegadc(c,T))
    #end Jacmusolvent
    
    def Jacmusolute(self,c,T):
        """
        Determination of the Jacobian of the chemical potential of the solute.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of the chemical potential of the solute.
        """
        return (1.-c)**2*(self.JacOmega(c,T)+c*self.JacdOmegadc(c,T))
    #end Jacmusolute
    
    def Jacd2Gdc2(self,c,T):
        """
        Determination of the Jacobian of the second derivative respact of c of
        the Gibbs free energy.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of of the second derivative respact of c of
            the Gibbs free energy.
        """
        return -2.*self.JacOmega(c,T)+2.*(1.-2.*c)*self.JacdOmegadc(c,T)+\
                c*(1.-c)*self.Jacd2Omegadc2(c,T)
    #end Jacd2Gdc2
    
    def Jacd3Gdc3(self,c,T):
        """
        Determination of the Jacobian of the third derivative respact of c of
        the Gibbs free energy.
        
        Parameters
        ----------
        c: Float
            Molar fraction of the solute.
        T: Float
            Temperature in K.
        
        Returns
        -------
        Array of floats of dimension (OrderRK+1,NTRK)
            Jacobian of of the third derivative respact of c of
            the Gibbs free energy.
        """
        return -6.*self.JacdOmegadc(c,T)+3.*(1.-2.*c)*self.Jacd2Omegadc2(c,T)+\
                c*(1.-c)*self.Jacd3Omegadc3(c,T)
    #end Jacd3Gdc3

    
    def fRK(self,x,T,x1,x2,xc,Tc,x1mono,x2mono,Tmono,weightmono):
        """
        This function determines the function to minimize to compute a new set
        of Radlich-Kister coefficients.
        
        Parameters
        ----------
        x : Array of float
            Unknowns corresponding to the coefficient of A in the excess potential.
            
        T : array of Float
            Temperature of the thermodynamic system in equilibrium in K.
        x1 : Array of Float
            Molar concentration of the solute on the left of the binonal curve.
        x2 : Array of Float
            Molar concentration of the solute on the left of the binonal curve.

        Returns
        -------
        f : Float
            Root Mean square error of binonal curve, critical point and 
            monotectic line.
        """
        
        # Reshape of unknows to be a matrix
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                self.A[i,j]=x[j+self.NTRK*i]
            #end for
        #end for
        
        # Determination of the function to minimize
        N=np.size(T)
        f=0.
        for i in range(N):
            muA1=self.musolvent(T[i],x1[i])
            muA2=self.musolvent(T[i],x2[i])
            muB1=self.musolute(T[i],x1[i])
            muB2=self.musolute(T[i],x2[i])
            f+=(muA1-muA2)**2+(muB1-muB2)**2
        #end for
        
        # Addition of the critical conditions
        # -----------------------------------
                
        f+=(self.d2freeenergydc2(xc,Tc))**2
        f+=(self.d3freeenergydc3(xc,Tc))**2
        
        # Addition of the monotectic conditions
        # -------------------------------------
        
        if (self.nature=='monotectic'):
            f+=weightmono*self.musolvent(Tmono,x1mono)**2
            f+=weightmono*self.musolvent(Tmono,x2mono)**2
            f+=weightmono*(self.musolute(Tmono,x1mono)-self.musolute(Tmono,x2mono))**2
        #endif
        
        # Return of the quadratic difference
        return f
    #end fRK

    def jacfRK(self,x,T,x1,x2,xc,Tc,x1mono,x2mono,Tmono,weightmono):
        """
        This function determines the jacobienne of the function to minimize to compute a new set
        of Redlich-Kister coefficients.
        
        Parameters
        ----------
        x : Array of float
            Unknowns corresponding to the coefficient of A in the excess potential.
            
        T : array of Float
            Temperature of the thermodynamic system in equilibrium in K.
        x1 : Array of Float
            Molar concentration of the solute on the left of the binonal curve.
        x2 : Array of Float
            Molar concentration of the solute on the left of the binonal curve.

        Returns
        -------
        jacf : Array of Float
            Quadratic difference between each chemical potential of A and B on the
            binonal curve.
        """
        
        # Reshape of unknows to be a matrix
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                self.A[i,j]=x[j+self.NTRK*i]
            #end for
        #end for
        
        # Determination of the jacobian of the function to minimize
        N=np.size(T)
        jacobf=np.zeros(self.NTRK*(self.OrderRK+1))
        
        for l in range(N):
            muA1=self.musolvent(T[l],x1[l])
            muA2=self.musolvent(T[l],x2[l])
            muB1=self.musolute(T[l],x1[l])
            muB2=self.musolute(T[l],x2[l])
            dmuAdx1=self.Jacmusolvent(x1[l],T[l])
            dmuAdx2=self.Jacmusolvent(x2[l],T[l])
            dmuBdx1=self.Jacmusolute(x1[l],T[l])
            dmuBdx2=self.Jacmusolute(x2[l],T[l])
            for i in range(self.OrderRK+1):
                for j in range(self.NTRK):
                    jacobf[j+self.NTRK*i]+=(muA1-muA2)*(dmuAdx1[i,j]-dmuAdx2[i,j])
                    jacobf[j+self.NTRK*i]+=(muB1-muB2)*(dmuBdx1[i,j]-dmuBdx2[i,j])
                #end for
            #end for
        #end for
        
        # Addition of the critical conditions
        # -----------------------------------
        d2Gdc2=self.d2freeenergydc2(xc,Tc)
        jd2Gdc2=self.Jacd2Gdc2(xc,Tc)
        d3Gdc3=self.d3freeenergydc3(xc,Tc)
        jd3Gdc3=self.Jacd3Gdc3(xc,Tc)
        for i in range(self.OrderRK+1):
            for j in range(self.NTRK):
                jacobf[j+self.NTRK*i]+=d2Gdc2*jd2Gdc2[i,j]
                jacobf[j+self.NTRK*i]+=d3Gdc3*jd3Gdc3[i,j]
            #end for
        #end for
        
        
        # Addition of the monotectic conditions
        # -------------------------------------
        
        if (self.nature=='monotectic'):
            muA1=self.musolvent(Tmono,x1mono)
            muA2=self.musolvent(Tmono,x2mono)
            muB1=self.musolute(Tmono,x1mono)
            muB2=self.musolute(Tmono,x2mono)
            dmuAdx1=self.Jacmusolvent(x1mono,Tmono)
            dmuAdx2=self.Jacmusolvent(x2mono,Tmono)
            dmuBdx1=self.Jacmusolute(x1mono,Tmono)
            dmuBdx2=self.Jacmusolute(x2mono,Tmono)
            for i in range(self.OrderRK+1):
                for j in range(self.NTRK):
                    jacobf[j+self.NTRK*i]+=weightmono*muA1*dmuAdx1[i,j]
                    jacobf[j+self.NTRK*i]+=weightmono*muA2*dmuAdx2[i,j]
                    jacobf[j+self.NTRK*i]+=weightmono*(muB1-muB2)*(dmuBdx1[i,j]-dmuBdx2[i,j])
                #end for
            #end for
        #end if
        
        # Return of the jacobian
        return 2.*jacobf
    #end jacfRK
    
    def OxygenNumber(self,oxide):
        """
        Determination of the number of oxygen in the list of oxides constituting
        the glass composition.
        
        Parameters
        ----------
        oxide: String
            name of the oxide.
        
        Returns
        -------
        nbO: Integer
            Number of oxygen in the oxide.
        cation: String
            Name of the cation in the oxide.
        nbcation: Interger
            Number of cation.
        """
        
        if (oxide.find('O')>0):
            if len(oxide)==(oxide.find('O')+1):
                nbO=1
            else:
                nbO=int(oxide[-1])
            #end if
        #end if
        
        if (nbO==1):
            cation=oxide[:-1]
        else:
            cation=oxide[:-2]
        #end if
        
        nbcation=1.
        for j in cation:
            if (j.isdigit()):
                nbcation=np.float64(j)
            #end if
        #end for
        
        if (nbcation>1.):
            cation=cation[:-1]
        #end if
        
        return nbO,cation,nbcation
    #end OxygenNumber
    
    def ComputeNBOsT(self,x):
        """
        Determination of the NBOsT of the binary system.
        
        """
        
        # Determination of the number of oxygen, the name of cation  and number
        # of cation in the solute
        
        nbOsolute,cationsolute,nbcationsolute=self.OxygenNumber(self.solute)
        
        # Computation of the atomic fraction of Si
        xatomSi=(1.-x)/(3.*(1.-x)+x*(nbOsolute+nbcationsolute))
        
        # Computation of the atomic fraction of O
        xatomO=(2.*(1.-x)+nbOsolute*x)/(3.*(1.-x)+x*(nbOsolute+nbcationsolute))
        
        # Computation and return of the NBO/T
        return (2.*xatomO-4.*xatomSi)/xatomSi
    #end ComputeNBOsT

#end class