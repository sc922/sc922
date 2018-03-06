# -*- coding: utf-8 -*-
"""
Created on Tue Aug 08 15:54:58 2017

@author: sc922
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.interpolate import interp1d

def PolyRI(wls,PolyPhase):
    wlsPolyData = np.array([400, 460, 520, 580, 640, 700, 760, 820, 880, 940, 1000])
    if(PolyPhase == 'cold'):
        nPolyData = np.array([1.35, 1.345,1.342, 1.339, 1.3365, 1.335, 1.334, 1.333, 1.3315, 1.331, 1.3307])
    elif(PolyPhase == 'hot'):
        nPolyData = np.array([1.446, 1.435,1.426, 1.420, 1.415, 1.413,  1.410, 1.407, 1.405,  1.404, 1.403])
    elif(PolyPhase == 'dry'):
        nPolyData = np.array([1.446, 1.435,1.426, 1.420, 1.415, 1.413,  1.410, 1.407, 1.405,  1.404, 1.403])+0.12
    else:
        print("PolyPhase must be 'cold', 'hot', or 'dry'")
        raise Exception('exit')
    nPoly = np.interp(wls,wlsPolyData,nPolyData)
    cPoly = nPoly+1e-6j
    return cPoly
    
def AuRI(wls):
    AuData = np.genfromtxt('AuRI.csv', delimiter=',')
    nAu = np.interp(wls,AuData[::-1,0],AuData[::-1,1])
    kAu = np.interp(wls,AuData[::-1,0],AuData[::-1,2])
    cAu = nAu+kAu*1j
    return cAu
    
def SubstrateRI(wls,SubstrateMat):
    if(SubstrateMat == 'Si'):
        SubstrateData = np.genfromtxt('SiRI.csv', delimiter=',')
    elif(SubstrateMat == 'SiO2'):
        SubstrateData = np.genfromtxt('SiO2RI.csv', delimiter=',')
    #elif(SubstrateMat == 'Al2O3'):
    #    SubstrateData = np.genfromtxt('Al2O3RI.csv', delimiter=',')
    elif(SubstrateMat == 'H2O'):
        SubstrateData = np.genfromtxt('H2O.csv', delimiter=',')
    else:
        print("PolyPhase must be 'Si', 'SiO2', 'Al2O3', 'H2O'")
        raise Exception('exit')
    nSubstrate= np.interp(wls,SubstrateData[::-1,0],SubstrateData[::-1,1])
    kSubstrate = np.interp(wls,SubstrateData[::-1,0],SubstrateData[::-1,2])
    cSubstrate = nSubstrate + kSubstrate*1j
    return cSubstrate
    
def EffectiveMedium(EMA_Model,ff,cPoly,cAu):
    epsPoly = cPoly**2
    epsAu = cAu**2
    if(EMA_Model == 'MG'):
        epsEMA = epsPoly*(2*ff*(epsAu-epsPoly)+epsAu+2*epsPoly)/(2*epsPoly+epsAu+ff*(epsPoly-epsAu))
        cEMA = np.sqrt(1./2.*((np.sqrt(epsEMA.real**2+epsEMA.imag**2))+epsEMA.real))+np.sqrt(1./2.*(np.sqrt(epsEMA.real**2+epsEMA.imag**2)-epsEMA.real))*1j
    else:
        print("EMA_Model must be 'MG'")#, 'Brugg', 'LL', 'Par',or 'Ser'")
        raise Exception('exit')    
    cEMA = cEMA.real + cEMA.imag*1j    
    return cEMA

def FresnelCoef(n0,n1,theta0):
    theta1 = np.arcsin(n0/n1*np.sin(theta0))
    rp = (n1*np.cos(theta0)-n0*np.cos(theta1))/(n0*np.cos(theta1)+n1*np.cos(theta0))
    rs = (n0*np.cos(theta0)-n1*np.cos(theta1))/(n0*np.cos(theta0)+n1*np.cos(theta1))
    tp = 2*n0*np.cos(theta0)/(n0*np.cos(theta1)+n1*np.cos(theta0))
    ts = 2*n0*np.cos(theta0)/(n0*np.cos(theta0)+n1*np.cos(theta1))
    return(rp,rs,tp,ts)
    
def RouardsMethod(wls,n0,cEMA,cSubstrate,inc_ang,thickness):
    trans_ang = np.arcsin(n0/cEMA*np.sin(inc_ang))
    trans_ang2 = np.arcsin(cEMA/cSubstrate*np.sin(trans_ang))
    r01p,r01s,t01p,t01s=FresnelCoef(n0,cEMA,inc_ang)
    r12p,r12s,t12p,t12s=FresnelCoef(cEMA,cSubstrate,trans_ang)
    phi = 2*(2*np.pi/wls)*(cEMA)*thickness*np.cos(trans_ang)

    rcoef_p = (r01p+r12p*np.exp(1j*phi))/(1+r01p*r12p*np.exp(1j*phi))
    rcoef_s = (r01s+r12s*np.exp(1j*phi))/(1+r01s*r12s*np.exp(1j*phi))
    
    tcoef_p = (t01p*t12p*np.exp(1j*phi))/(1+t01p*t12p*np.exp(1j*phi))
    tcoef_s = (t01s*t12s*np.exp(1j*phi))/(1+t01s*t12s*np.exp(1j*phi))
    
    Rp = rcoef_p*np.conj(rcoef_p)
    Rs = rcoef_s*np.conj(rcoef_s)
    
    Tp = cSubstrate/n0*np.cos(trans_ang2)/np.cos(inc_ang)*(tcoef_p*np.conj(tcoef_p))
    Ts = cSubstrate/n0*np.cos(trans_ang2)/np.cos(inc_ang)*(tcoef_s*np.conj(tcoef_s))
    
    return Rp.real,Rs.real,Tp.real,Ts.real
    #return (rcoefp.real,rcoef_p.real,rcoef_p.imag,rcoef_p.imag)
    #return rcoef_p.imag,-rcoef_s.imag,tcoef_p.imag,tcoef_s.imag
    
#Main Script    
plt.close("all")    
#Parameters
n0 = 1
inc_ang = 0*np.pi/180
thickness = 1000 #nm
EMA_Model = 'MG'

num_ff = 301
num_wls = 301
ff  = np.linspace(0,0.7405,num_ff)
wls = np.linspace(420,800,num_wls)

cPoly = PolyRI(wls,'dry') #Choose 'cold','hot', or 'dry'
cSubstrate = SubstrateRI(wls,'Si') #Choose 'Si', 'SiO2', 'Al2O3', 'H2O'
cAu = AuRI(wls)

M_cEMA = np.zeros((num_wls,num_ff),dtype=complex)
for i in range (0,num_ff-1):
    cEMA = EffectiveMedium(EMA_Model,ff[i],cPoly,cAu)
    M_cEMA[:,i] = cEMA
        
M_nEMA = M_cEMA.real
M_kEMA = M_cEMA.imag

#Initialize matrices
M_RpEMA = np.zeros((num_wls,num_ff))
M_RsEMA = np.array(M_RpEMA) 
M_TpEMA = np.array(M_RpEMA)
M_TsEMA = np.array(M_RpEMA)

#Solving reflectance and transmittance
for i in range (0,num_ff-1):
    cEMA = M_cEMA[:,i]
    M_RpEMA[:,i],M_RsEMA[:,i],M_TpEMA[:,i],M_TsEMA[:,i] = RouardsMethod(wls,n0,cEMA,cSubstrate,inc_ang,thickness)

Reflectance = (M_RpEMA+M_RsEMA)/2
Transmittance = (M_TpEMA+M_TsEMA)/2

#Plotting Refrcative indices, Reflectance and transmittance
X, Y = np.meshgrid(ff, wls)
fig,((ax1, ax2), (ax3, ax4))  = plt.subplots(2,2,figsize=(20,12))
z1 = ax1.pcolor(X,Y, M_nEMA, cmap='viridis' )
plt.colorbar(z1,ax=ax1)
z2 = ax2.pcolor(X,Y, M_kEMA, cmap='viridis' )
plt.colorbar(z2,ax=ax2)
ax1.set_xlim([min(ff),max(ff)])
ax2.set_xlim([min(ff),max(ff)])
ax1.set_ylim([min(wls),max(wls)])
ax2.set_ylim([min(wls),max(wls)])
ax1.set_title('Real Refractive Index, '+EMA_Model)
ax2.set_title('Imaginary Refractive Index, '+EMA_Model)

z3 = ax3.pcolor(X,Y, Reflectance, cmap='viridis' )
plt.colorbar(z3,ax=ax3)
z4 = ax4.pcolor(X,Y, Transmittance, cmap='viridis' )
plt.colorbar(z4,ax=ax4)
ax3.set_xlim([min(ff),max(ff)])
ax4.set_xlim([min(ff),max(ff)])
ax3.set_ylim([min(wls),max(wls)])
ax4.set_ylim([min(wls),max(wls)])
ax3.set_title('Reflectance, '+EMA_Model)
ax4.set_title('Transmittance, '+EMA_Model)