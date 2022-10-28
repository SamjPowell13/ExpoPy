import numpy as np
import scipy.optimize as opt

mag,  zdegrees, aperture, T,R,QE,bandwidth = np.loadtxt("\\Users\\samjp\\Documents\\Year 3\\Observatory\\Target_data.csv", unpack = True, delimiter = ",")
z = zdegrees/360*np.pi*2

def atmo_extinction(k, z):
    """Calculates the atmospheric extinction given the wav. dependent 
    extinction coefficient and the zenith distance of the object"""
    X = 1/np.cos(z) # airmass
    return 10 ** (-0.4 * k * X) 

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


Ns = source_countrate(mag, 0.2, z, aperture, T, R, QE, bandwidth)
sky_countrate_value = sky_countrate(17,12 ,aperture, T, R, QE, bandwidth)
Nb = background_countrate(sky_countrate_value, 0.2, 16)

noise = np.sqrt(Ns+2*Nb)
SNR = Ns*60/np.sqrt(Ns*60+2*Nb*60)
print("SNR for given magnitudes on C14 are:",SNR)
SskyNeb = 400*3600
npixNeb = SskyNeb*16/12

sky_countrate_valueNeb = sky_countrate(17, SskyNeb, aperture, T, R, QE, bandwidth)
NbNeb = background_countrate(sky_countrate_valueNeb, 0.2, npixNeb)


tKepler = expotime(3333, Ns[0], Nb[0])
tOri = expotime(100, Ns[1], Nb[1])
tHAT = expotime(40, Ns[2], Nb[2])
tGaia = expotime(100, Ns[3], Nb[3])
tSN = expotime(10, Ns[4], Nb[4])
tNeb = expotime(10, Ns[5], NbNeb[5])
print("Supernova exposure time is:",tSN)
print("Nebula exposure time is:",tNeb)
print("Gaia16bnz exposure time is:",tGaia)
print("V2841 Ori exposure time is:",tOri)
print("HAT-P-54b exposure time is:",tHAT)
print("Kepler 1530 exposure time is:",tKepler)

print(1/0.0003)

 


