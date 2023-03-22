#In this file, I will determine which SEP field a certain set of sources are located in, and the exact configuration of their observation (observed once, or twice?)

#In the clusters version of this file, I will see what sky tiles to put the cluster targets in.

#In version1, I will repeat this, but for the new regions for the feb22 run. 

import math
from astropy.io import fits
from astropy.table import Table
import os

CobsRA = [92.3301558180629, 89.75234057721435, 87.91299627794228, 93.12441226666928, 94.95790067508487, 91.94321419933092, 87.26050805857497, 84.9671513072304, 87.74470210202175]
CobsDec = [-66.95389859711848, -65.56576742094907, -67.13235116259148, -68.15930317079865, -66.29420073504171, -64.72803323753627, -64.9072191697707, -66.67071859250139, -68.36459073645581]


F = fits.open('SEP_members_targets.210826.fits')
RAo = F[1].data['RA']
Deco = F[1].data['Dec']
rmag = F[1].data['rmag']
Prty = [3]*len(RAo)
Nam = [0]*len(RAo)
for i in range(len(RAo)):
    Nam[i] = 'erocluster'+str(i+1)

#Determine the decimal degree RA and dec from this:
#RAo, Deco = [0]*len(RAoh), [0]*len(RAoh)
#for i in range(len(RAoh)):
#    RAo[i] = float(RAoh[i].split()[0])*(360/24) + float(RAoh[i].split()[1])*(360/(24*60)) + float(RAoh[i].split()[2])*(360/(24*60*60))
#    Deco[i] = float(Decod[i].split()[0]) - float(Decod[i].split()[1])/60 - float(Decod[i].split()[2])/3600   #Assumes that the Dec is negative.

#print('Decimal RA, Dec:')
#print('RA = ', RAo)
#print('Dec = ', Deco)    

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

def erosnamcon(ra, dec):
    hr = math.floor(ra/15)
    mr = math.floor(((ra/15) - hr)*60)
    sr = round((ra/15 - hr)*3600 - mr*60, 1)
    if hr < 10:
        hr1 = '0'+str(hr)
    else:
        hr1 = str(hr)
    if mr < 10:
        mr1 = '0'+str(mr)
    else:
        mr1 = str(mr)
    if sr < 10:
        sr1 = '0'+str(sr)
    else:
        sr1 = str(sr)
    hd = math.floor(abs(dec))
    md = math.floor((abs(dec) - hd)*60)
    sd = round((abs(dec) - hd)*3600 - md*60)
    if hd < 10:
        hd1 = '-0'+str(hd)
    else:
        hd1 = '-'+str(hd)
    if md < 10:
        md1 = '0'+str(md)
    else:
        md1 = str(md)
    if sd < 10:
        sd1 = '0'+str(sd)
    else:
        sd1 = str(sd)

    #For AGN, I will instead use U, rather than t:
    nam = 'eRASSUJ'+hr1+mr1+sr1+hd1+md1+sd1
    return nam

#Check whether I can associate the hh mm ss.s and dd mm ss.s to degrees, and vice versa:
#print(RAoh[0], RAo[0], radegtohms(RAo[0]))
#print(Decod[0], Deco[0], decdegtodms(Deco[0]))

#Print out the hh mm ss, dd mm ss of all the centers of the SEP regions of interest:
print('For the 9 centers of fields we will be observing now:')
RA, Dec = [89.55933908896876, 85.58109939350186, 85.70617035492553, 90.24172763147342, 92.00231465201824, 87.75361963217966], [-64.56901953351532, -65.66592132133434, -67.66235231957467, -67.55626141512766, -65.9686064529718, -66.12743183064512]
RA, Dec = [92.3301558180629, 89.75234057721435, 87.91299627794228, 93.12441226666928, 94.95790067508487, 91.94321419933092, 87.26050805857497, 84.9671513072304, 87.74470210202175], [-66.95389859711848, -65.56576742094907, -67.13235116259148, -68.15930317079865, -66.29420073504171, -64.72803323753627, -64.9072191697707, -66.67071859250139, -68.36459073645581]

for i in range(9):
    print(radegtohms(RA[i]))
    print(decdegtodms(Dec[i]))

#Now associate these sources with the 6 sky fields we will be observing:
for i in range(7, 16): #because I want it to start with SEP4, and then continue to SEP9.
    vars()['Nam'+str(i)] = [] 
    vars()['RAo'+str(i)] = []
    vars()['RAohms'+str(i)] = []
    vars()['Deco'+str(i)] = []
    vars()['Decodms'+str(i)] = []
    vars()['rmag'+str(i)] = []
    vars()['Prty'+str(i)] = []

for i in range(len(RAo)):
    iSF = [] #contains list of SEPfields that the source is in.
    for j in range(7, 16):
        if angbet([RAo[i], Deco[i]], [CobsRA[j-7], CobsDec[j-7]]) < 60 * 60:
            iSF.append(j)
    print('Source', i, 'lies in the following SEP fields:', iSF)
    print('It has RA in decimals and hh:mm:ss of:', RAo[i], radegtohms(RAo[i]))
    print('It has Dec in decimals and dd:mm:ss of:', Deco[i], decdegtodms(Deco[i]))
    #If a source is sufficiently bright, then only observe it once, otherwise observe it twice.But not more than twice.
    if len(iSF) == 1:
        vars()['Nam'+str(iSF[0])].append(Nam[i])
        vars()['RAo'+str(iSF[0])].append(RAo[i])
        vars()['RAohms'+str(iSF[0])].append(radegtohms(RAo[i]))
        vars()['Deco'+str(iSF[0])].append(Deco[i])
        vars()['Decodms'+str(iSF[0])].append(decdegtodms(Deco[i]))
        vars()['rmag'+str(iSF[0])].append(rmag[i])
        vars()['Prty'+str(iSF[0])].append(Prty[i])
        print('Associated with SEP field:', iSF[0])
    elif len(iSF) > 1:
        if rmag[i] < 20: # If it is sufficiently bright, only observe once.
#        if r[i] < 0: # For WDs we need all the regions to contain them.           
            iSF = [iSF[i%len(iSF)]]
            print('Source lying in more than one field, selected to be placed in', iSF, 'because of its rmag of ', rmag[i])
        else:
            iSF = [iSF[i%len(iSF)], iSF[(i+1)%len(iSF)]]
            print('Source lying in more than one field, selected to be placed in', iSF, 'because of its rmag of ', rmag[i])
        for k in range(len(iSF)):
            vars()['Nam'+str(iSF[k])].append(Nam[i])
            vars()['RAo'+str(iSF[k])].append(RAo[i])
            vars()['RAohms'+str(iSF[k])].append(radegtohms(RAo[i]))
            vars()['Deco'+str(iSF[k])].append(Deco[i])
            vars()['Decodms'+str(iSF[k])].append(decdegtodms(Deco[i]))
            vars()['rmag'+str(iSF[k])].append(rmag[i])
            vars()['Prty'+str(iSF[k])].append(Prty[i])

#To make sure, delete duplicate sources:
for i in range(7, 16):
    k = 1
    while k < len(vars()['RAo'+str(i)]):
        hd = 0
        for j in range(k):
            if vars()['RAo'+str(i)][k] == vars()['RAo'+str(i)][j]:
                if vars()['Deco'+str(i)][k] == vars()['Deco'+str(i)][j]:
                    print('Duplicate source, deleting higher integer one')
                    del vars()['Nam'+str(i)][k], vars()['RAo'+str(i)][k], vars()['RAohms'+str(i)][k], vars()['Deco'+str(i)][k], vars()['Decodms'+str(i)][k], vars()['rmag'+str(i)][k], vars()['Prty'+str(i)][k]
                    hd = 1
                    break
        if hd == 0:
            k += 1

#Check once more, whether there are any duplicate sources left:
print('\n\n')

for i in range(7, 16):
    for j in range(1, len(vars()['RAo'+str(i)])):
        for k in range(j):
            if vars()['RAo'+str(i)][j] == vars()['RAo'+str(i)][k]:
                if vars()['Deco'+str(i)][j] == vars()['Deco'+str(i)][k]:
                    print('Duplicate source in SEP'+str(i)+':')
                    print(vars()['Nam'+str(i)][j], vars()['Nam'+str(i)][k], vars()['RAo'+str(i)][j], vars()['RAo'+str(i)][k], vars()['Deco'+str(i)][j], vars()['Deco'+str(i)][k])
            
for i in range(7, 16):
    print('\n\nSEP'+str(i)+'\n')
#    print('# \tR. Ascention \tDeclination \t\t\t\tProg \tProper Motion \tComments')
#    print('# Name hh mm ss.sss \tdd mm ss.sss \t\t\tmag \tID \tra \tdec')
    for j in range(len(vars()['RAo'+str(i)])):
        print(str(vars()['Nam'+str(i)][j]) + '\t' + str(vars()['RAohms'+str(i)][j]) + '\t' + str(vars()['Decodms'+str(i)][j]) + '\tP\t' + str(vars()['Prty'+str(i)][j]) + '\t' + str(vars()['rmag'+str(i)][j]) + '\t0\t0\t0')


