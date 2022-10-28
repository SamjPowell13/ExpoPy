# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:08:29 2022

@author: samjp
"""
import numpy as np
z = np.pi*2
mag = np.array([9.0,13.2,17.5])
C14 = np.array([35,0.8,0.9,0.55,1000])
TL60 = np.array([60,0.9,0.9,0.55,1000])

def atmo_extinction(k, z):
    """Calculates the atmospheric extinction given the wav. dependent 
    extinction coefficient and the zenith distance of the object"""
    X = 1/np.cos(z) # airmass
    return 10 ** (-0.4 * k * 1.5) 

def sourceflux(mag):
    """ Calculates the photon flux at the top of the TOA"""
    Fvega = 1000
    mvega = 0
    return Fvega * 10 ** (-0.4 * (mag - mvega))

def source_countrate(mag, k, z, aperture, T, R, QE, bandwidth):
    """Calculates number of counts per unit time. Input apperture in cm!!
    bandwidth in Armstrongs"""
    F = sourceflux(mag)
    a = atmo_extinction(k, z)
    A = np.pi * (aperture/2) ** 2
    return F * a * A * T * R * QE * bandwidth

def sky_countrate(msky,Ssky, aperture, T, R, QE, bandwidth):
    F = sourceflux(msky)
    A = np.pi * (aperture/2) ** 2
    return F * A * T * R * QE * bandwidth *Ssky

def background_countrate(sky_countrate, dark_countrate, npix):
    return sky_countrate + dark_countrate * npix
    

def expotime(SNR, source_countrate, sky_countrate):
    return (SNR ** 2)*(source_countrate + 2 * sky_countrate )/(source_countrate**2)


NsC14 = source_countrate(mag, 0.2, z, C14[0], C14[1], C14[2], C14[3], C14[4])
sky_countrate_valueC14 = sky_countrate(17,12 ,C14[0], C14[1], C14[2], C14[3], C14[4])
NbC14 = background_countrate(sky_countrate_valueC14, 0.2, 16)#



print("Source photo-electron flux for C14:",NsC14)

NsTL60 = source_countrate(mag, 0.2, z, TL60[0], TL60[1], TL60[2], TL60[3], TL60[4])
sky_countrate_valueTL60 = sky_countrate(19,8, TL60[0], TL60[1], TL60[2], TL60[3], TL60[4])
NbTL60= background_countrate(sky_countrate_valueTL60, 0.2, 20)
tTL60 = expotime(333, NsTL60, NbTL60)

print("Source photo-electron flux for TL60:",NsTL60)
print("Nb for C14 is:",NbC14*60)
print("Ns for C14 is:",NsC14*60)
print(sky_countrate_valueC14*60/(NbC14*60))
print(0.2)

SNR = NsC14*60/np.sqrt(NsC14*60+2*NbC14*60)
print("SNR for given magnitudes on C14 are:",SNR)

SNRTL60 = NsTL60*60/np.sqrt(NsTL60*60 + 2*NbTL60*60)
print("SNR for TL60 is:", SNRTL60)

tC14 = expotime(10, NsC14, NbC14)
tTL60 = expotime(10,NsTL60,NbTL60)
print("Exposure times on C14:",tC14)
print("Exposure times on TL60:",tTL60)

