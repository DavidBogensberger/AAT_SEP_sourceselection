#In this script, I will try to identify stars and AGN from the catalogue that Mara suggested. Actually, there are two, and I want to see the comparison between them.

#In version1, I will select the WDs for the 9 fields selected for the feb22 run. 

import numpy as np
import math
from astropy.io import fits
from astropy.table import Table
import os
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title('ABC123 vs $\mathrm{ABC123}^{123}$')

#CobsRA = [89.55933908896876, 85.58109939350186, 85.70617035492553, 90.24172763147342, 92.00231465201824, 87.75361963217966]
#CobsDec = [-64.56901953351532, -65.66592132133434, -67.66235231957467, -67.55626141512766, -65.9686064529718, -66.12743183064512]
#For Chandreyee's and Frank's LMC field, the center is located at:
#CobsRA = [83.6757208333]
#CobsDec = [-67.20414722222]
CobsRA = [92.3301558180629, 89.75234057721435, 87.91299627794228, 93.12441226666928, 94.95790067508487, 91.94321419933092, 87.26050805857497, 84.9671513072304, 87.74470210202175]
CobsDec = [-66.95389859711848, -65.56576742094907, -67.13235116259148, -68.15930317079865, -66.29420073504171, -64.72803323753627, -64.9072191697707, -66.67071859250139, -68.36459073645581]


def angbet(A1, A2):
    #Returns the separation of two ra, dec pointings, in terms of arcsec. 
    A = [math.pi/180 * A1[i] for i in range(2)]
    B = [math.pi/180 * A2[i] for i in range(2)]
    Ab0 = math.cos(A[1])*math.cos(B[1])*math.cos(B[0]-A[0])+math.sin(A[1])*math.sin(B[1])
    if Ab0 > 1.0:
        print('problem calculating angle: ', A, B, Ab0)
        Ab0 = 1.0
    Ab = math.acos(Ab0)
    Ab *= 180 * 60 * 60 / math.pi
    return Ab

def radegtohms(ra):
    #Assuming RA is positive
    h = math.floor(ra/15)
    m = math.floor(((ra/15) - h)*60)
    s = round((ra/15 - h)*3600 - m*60, 3)
    if h < 10:
        h1 = '+0'+str(h)
    else:
        h1 = '+'+str(h)
    if m < 10:
        m1 = '0'+str(m)
    else:
        m1 = str(m)
    if s < 10:
        s1 = '0'+str(s)
    else:
        s1 = str(s)
    rhms = str(h1+' '+m1+' '+s1)
    return rhms

def decdegtodms(dec):
    #Assuming Dec is negative
    h = math.floor(abs(dec))
    m = math.floor((abs(dec) - h)*60)
    s = round((abs(dec) - h)*3600 - m*60, 3)
    if h < 10:
        h1 = '-0'+str(h)
    else:
        h1 = '-'+str(h)
    if m < 10:
        m1 = '0'+str(m)
    else:
        m1 = str(m)
    if s < 10:
        s1 = '0'+str(s)
    else:
        s1 = str(s)
    rhms = str(h1+' '+m1+' '+s1)
    return rhms

#Open the cataogue:
F = fits.open('WD_gaia.fits')
Nam = F[1].data['WD']
RA = F[1].data['_RAJ2000']
Dec = F[1].data['_DEJ2000']
rmag = F[1].data['RPmag'] 
pmRA = F[1].data['pmRA']
pmDec = F[1].data['pmDE']
F.close()

N = len(Nam)

#for i in range(4, 10):
for i in range(len(CobsRA)):
    vars()['Nwd'+str(i)], vars()['RAwd'+str(i)], vars()['Decwd'+str(i)], vars()['rmagwd'+str(i)], vars()['pmRwd'+str(i)], vars()['pmDwd'+str(i)], vars()['Tp'+str(i)] = [], [], [], [], [], [], []
#Find which WDs lie in the 6 fields:

#for i in range(4, 10):
#When dealing with the LMC field: 
for i in range(len(CobsRA)):
    for j in range(N):
        #if angbet([CobsRA[i-4], CobsDec[i-4]], [RA[j], Dec[j]]) < 60*60:
        if angbet([CobsRA[i], CobsDec[i]], [RA[j], Dec[j]]) < 60*60:
            vars()['Nwd'+str(i)].append(Nam[j])
            vars()['RAwd'+str(i)].append(RA[j])
            vars()['Decwd'+str(i)].append(Dec[j])
            vars()['rmagwd'+str(i)].append(rmag[j])
            vars()['pmRwd'+str(i)].append(round(pmRA[j] / 1000, 6))
            vars()['pmDwd'+str(i)].append(round(pmDec[j] / 1000, 6))
            #vars()['Tp'+str(i)].append(Tpwd[i])

#for i in range(4, 10):
for i in range(len(CobsRA)):    
    print('\n\nWhite Dwarfs that lie in SEP field'+str(i+7)+':')
    for j in range(len(vars()['Nwd'+str(i)])):
        print(vars()['Nwd'+str(i)][j] + '\t' + radegtohms(vars()['RAwd'+str(i)][j]) + '\t' + decdegtodms(vars()['Decwd'+str(i)][j]) + '\tP\t9\t' + str(vars()['rmagwd'+str(i)][j]) + '\t0\t' + str(vars()['pmRwd'+str(i)][j]) + '\t' + str(vars()['pmDwd'+str(i)][j]))


