#In this code, I will generate the optical positions and all the other information I need to observe my list of variable sources in the 6 SEP fields we will be focusing on to begin with.

#In version 1 of the code, I will keep all variable sources, even if they lie in an overlapping region featuring two sources. I will manually remove duplicate sources if they lie in two regions, both of which will be observed in the same night.

#In version 2 of this code, I will update the list of variable sources to those that lie above 3 sigma only. But these are new sources based on an improved variable source detection pipeline, and improved light curves. A lot of the code will be rewritten, as it is based on results of my analysis and matching of the variable sources. But most of this code is still based on Varsrcprep_AAT1.py of course.

#In version 3, I will update the AGN/star selection. I will only exclude the sources that have a parallax of larger than 5 sigma. I will keep all other sources, just to be sure that I don't exclude sources that could be AGNs. I decided against changing the rmag limits to 16-21.5, as the distribution of sources really is centered around 17-22.5

#In version 4, I will update the code to generate the data for the feb22 run of obsevations. This means changing the fields observed. 

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

    nam = 'eRASStJ'+hr1+mr1+sr1+hd1+md1+sd1
    return nam

#Open up all the relevant properties of the variable sources.

#F = fits.open('/home/david/Documents/SharedDocuments/eROSITA/eRASS/946eRASS123/FVS1a_prop_varsrc3_allvarsrc_matchnonvar_e123_3sigma_optAAT.fits')
F = fits.open('/home/david/Documents/SharedDocuments/eROSITA/eRASS/946eRASS123/FVS2_prop_varsrc3_allvarsrc_matchnonvar_e123_3sigma_optAAT_3.fits')
RAx = F[1].data['RAx7'].tolist()
Decx = F[1].data['Decx7'].tolist()
RAo = F[1].data['RAo'].tolist()
Deco = F[1].data['Deco'].tolist()
#prlxoe = F[1].data['propmot_div_error'].tolist()
prlxoe = F[1].data['parallax_div_error'].tolist()
AoS = F[1].data['AGNorStar_photandspec'].tolist() # This is the criterion I am going to base this selection on. It is very generous to the AGNs.
QOP = F[1].data['AAT_QOP'].tolist()
rmag = F[1].data['rmag'].tolist()
SR = F[1].data['AAT_SEPreg'].tolist()
ID = F[1].data['AAT_ID'].tolist()
F.close()

#Remove sources that were identified as stars, with a QOP of 6:

N = len(RAo)
print('\nNumber of variable sources before removing spectroscopically identified stars from the list:', N)

k = 0
while k < len(SR):
    if QOP[k] == 6:
        print('Variable source with a stellar spectrum removed from list of variable sources to reobserve:', SR[k], ID[k])
        del RAx[k], Decx[k], RAo[k], Deco[k], prlxoe[k], AoS[k], QOP[k], rmag[k], SR[k], ID[k]
    else:
        k += 1
N = len(RAo)
print('Number of variable sources after spectroscopically identified stars were removed:', N)

#Now remove sources that have a parralax over error of more than 5.
print('\nNumber of sources before deleting sources with parallax over error large than 5:', len(RAo))
k = 0
while k < len(RAo):
    if prlxoe[k] > 5:
        del RAx[k], Decx[k], RAo[k], Deco[k], prlxoe[k], AoS[k], QOP[k], rmag[k], SR[k], ID[k]
    else:
        k += 1

print('Number of sources after deleting sources with parallax over error large than 5:', len(RAo))

#Determine what range of r magnitudes would contain the greatest number of variable sources.
#start off with a grid. That denotes the lower of the two magnitudes (so the brighter limit) defining the range. 
gr = np.arange(12, 18.51, 0.01)
#See how many sources there are within the starting point specified in the grid, and the end point, which is start point + 5.5.

nsir = [0]*len(gr)
bgr, nbgr = 18.51, 0 # saving the start/stop index where most sources are kept.
#start from the bottom up.
for i in range(len(gr)):
    for j in range(len(RAo)):
        if gr[len(gr) - 1 - i] < rmag[j] < gr[len(gr) - 1 - i] + 5.5:
            nsir[len(gr) - 1 - i] += 1
    if nsir[len(gr) - 1 - i] > nbgr:
        nbgr = nsir[len(gr) - 1 - i]
        bgr = gr[len(gr) - 1 - i]

print('The best choice of r magnitude range to keep the greatest number of variable sources is:', bgr, bgr+5.5, '. When choosing that range, we would keep ', nbgr, 'variable sources.')

#Find the number of variable sources between 16 and 17, and between 21.5 and 22.5:
n1617, n215225 = 0, 0
for i in range(len(RAo)):
    if 16 < rmag[i] < 17:
        n1617 += 1
    elif 21.5 < rmag[i] < 22.5:
        n215225 += 1

print('Number of variable sources between rmag 16 and 17, and between 21.5 and 22.5:', n1617, n215225)

#Next, remove sources that are too bright or too faint.
print('\nLength of list of variable sources before deleting too bright and too faint sources, or sources that have an r mag of nan or inf:', len(RAo))
k = 0
while k < len(RAo):
    if rmag[k] < 17:
        del RAx[k], Decx[k], RAo[k], Deco[k], prlxoe[k], AoS[k], QOP[k], rmag[k], SR[k], ID[k]
    elif rmag[k] > 22.5:
        del RAx[k], Decx[k], RAo[k], Deco[k], prlxoe[k], AoS[k], QOP[k], rmag[k], SR[k], ID[k]
    elif np.isnan(rmag[k]) + np.isinf(rmag[k]) > 0:
        del RAx[k], Decx[k], RAo[k], Deco[k], prlxoe[k], AoS[k], QOP[k], rmag[k], SR[k], ID[k]
    else:
        k += 1
print('Length of list of variable sources after deleting too bright and too faint sources, or sources that have an r mag of nan or inf:', len(RAo))

#Now assign the priorities. I start off with a priority of 9. If they were already previously observed, assign a lower priority to them. (-1). If they have an AoS of unknown, also subtract 1. If they have an AoS of star, subtract 2.

Prty = [9]*len(RAo)
for i in range(len(RAo)):
    if QOP[i] > 0:
        Prty[i] -= 1
    if AoS[i] == 1:
        Prty[i] -= 1 # unknown sources.
    elif AoS[i] == 0:
        Prty[i] -= 2 # stars 

#Generate names for all of these files as well:

Nam = [0]*len(RAo)
for i in range(len(RAo)):
    Nam[i] = erosnamcon(RAx[i], Decx[i]) # I will use the X-ray position, as this is an X-ray telescope position.

#test to see whether the deletion is going correctly, such that the relevant lists are equally long.
print('\n', len(RAo), len(Deco), len(rmag), len(Prty), len(Nam)) 

#Now associate these sources with the 6 sky fields we will be observing:
for i in range(7, 16): #because I want it to start with SEP4, and then continue to SEP9.
    vars()['Nam'+str(i)] = [] 
    vars()['RAo'+str(i)] = []
    vars()['RAohms'+str(i)] = []
    vars()['Deco'+str(i)] = []
    vars()['Decodms'+str(i)] = []
    vars()['rmag'+str(i)] = []
    vars()['Prty'+str(i)] = []
    vars()['nv215_225_'+str(i)], vars()['nv16_17_'+str(i)] = 0, 0
for i in range(len(RAo)):
    iSF = [] #contains list of SEPfields that the source is in.
    #Also make sure that the magnitude of the source is within the ranges that we specified. 
    #for j in range(4, 7):
    for j in range(7, 16):
        if angbet([RAo[i], Deco[i]], [CobsRA[j-7], CobsDec[j-7]]) < 60 * 60:
            #if 17.0 < rmag[i] < 22.5:
            iSF.append(j)
            #if 21.5 < rmag[i] < 22.5:
            #    vars()['nv215_225_'+str(j)] += 1
            #if 16 < rmag[i] < 17:
            #    vars()['nv16_17_'+str(j)] += 1
    #for j in range(7, 10):
    #    if angbet([RAo[i], Deco[i]], [CobsRA[j-4], CobsDec[j-4]]) < 60 * 60:
    #        if 16.0 < rmag[i] < 21.5:
    #            iSF.append(j)
    #        if 21.5 < rmag[i] < 22.5:
    #            vars()['nv215_225_'+str(j)] += 1
    #        if 16 < rmag[i] < 17:
    #            vars()['nv16_17_'+str(j)] += 1
    #print('Variable source', i, 'lies in the following SEP fields:', iSF)
    #print('It has RA in decimals and hh:mm:ss of:', RAo[i], radegtohms(RAo[i]))
    #print('It has Dec in decimals and dd:mm:ss of:', Deco[i], decdegtodms(Deco[i]))
    #Add the source to all the skytiles it is containd within. Even if it appears twice, that is good.
    #for k in range(len(iSF)):
    #    vars()['Nam'+str(iSF[k])].append(Nam[i])
    #    vars()['RAo'+str(iSF[k])].append(RAo[i])
    #    vars()['RAohms'+str(iSF[k])].append(radegtohms(RAo[i]))
    #    vars()['Deco'+str(iSF[k])].append(Deco[i])
    #    vars()['Decodms'+str(iSF[k])].append(decdegtodms(Deco[i]))
    #    vars()['rmag'+str(iSF[k])].append(round(rmag[i], 2))
    #    vars()['Prty'+str(iSF[k])].append(Prty[i])

    #No, only include sources once, unless they are too faint, then include them twice.
    if len(iSF) == 1:
        vars()['Nam'+str(iSF[0])].append(Nam[i])
        vars()['RAo'+str(iSF[0])].append(RAo[i])
        vars()['RAohms'+str(iSF[0])].append(radegtohms(RAo[i]))
        vars()['Deco'+str(iSF[0])].append(Deco[i])
        vars()['Decodms'+str(iSF[0])].append(decdegtodms(Deco[i]))
        vars()['rmag'+str(iSF[0])].append(round(rmag[i], 2))
        vars()['Prty'+str(iSF[0])].append(Prty[i])
        #print('Associated with SEP field:', iSF[0])
    elif len(iSF) > 1:
        if rmag[i] < 20: # If it is sufficiently bright, only observe once.
            iSF = [iSF[i%len(iSF)]]
            #print('Source lying in more than one field, selected to be placed in', iSF, 'because of its rmag of ', r[i])
        #else:
        #    iSF = [iSF[i%len(iSF)], iSF[(i+1)%len(iSF)]]
        #    print('Source lying in more than one field, selected to be placed in', iSF, 'because of its rmag of ', r[i])
        for k in range(len(iSF)):
            vars()['Nam'+str(iSF[k])].append(Nam[i])
            vars()['RAo'+str(iSF[k])].append(RAo[i])
            vars()['RAohms'+str(iSF[k])].append(radegtohms(RAo[i]))
            vars()['Deco'+str(iSF[k])].append(Deco[i])
            vars()['Decodms'+str(iSF[k])].append(decdegtodms(Deco[i]))
            vars()['rmag'+str(iSF[k])].append(round(rmag[i], 2))
            vars()['Prty'+str(iSF[k])].append(Prty[i])
    
    
#Check whether the regions are defined correctly by plotting them:

#Also check whether there are any sources within 30":
for i in range(7, 16):
    for j in range(len(vars()['RAo'+str(i)])):
        for k in range(j):
            if angbet([vars()['RAo'+str(i)][j], vars()['Deco'+str(i)][j]], [vars()['RAo'+str(i)][k], vars()['Deco'+str(i)][k]]) < 30:
                print('Two sources in mas of SEP', i, ' are too close together:', j, k, vars()['RAo'+str(i)][j], vars()['RAo'+str(i)][k], vars()['Deco'+str(i)][j], vars()['Deco'+str(i)][k], angbet([vars()['RAo'+str(i)][j], vars()['Deco'+str(i)][j]], [vars()['RAo'+str(i)][k], vars()['Deco'+str(i)][k]]), 30)
        
        
#Add in the SEP AAT observing regions.
RAao = [90, 90.49395819453254, 94.7218336866107, 93.91137202703467, 89.55933908896876, 85.58109939350186, 85.70617035492553] #already observed
Decao = [-66.560708333333, -68.55181449691935, -67.32389564748709, -65.35130972297905, -64.56901953351532, -65.66592132133434, -67.66235231957467]
RAwo = [92.3301558180629, 89.75234057721435, 87.91299627794228, 93.12441226666928, 94.95790067508487, 91.94321419933092, 87.26050805857497, 84.9671513072304, 87.74470210202175] #Will observe these in future. 
Decwo = [-66.95389859711848, -65.56576742094907, -67.13235116259148, -68.15930317079865, -66.29420073504171, -64.72803323753627, -64.9072191697707, -66.67071859250139, -68.36459073645581]

Crdao = [[0]*1000 for i in range(7)]
Crrpao = [[0]*1000 for i in range(7)]
Crrnao = [[0]*1000 for i in range(7)]
Crdwo = [[0]*1000 for i in range(9)]
Crrpwo = [[0]*1000 for i in range(9)]
Crrnwo = [[0]*1000 for i in range(9)]

for i in range(7):
    Crdao[i] = [Decao[i]+1.0*math.sin(math.pi*j/999-0.5*math.pi) for j in range(1000)]
    for j in range(1000):
        a = (math.cos(1.0*math.pi/180)-math.sin(Decao[i]*math.pi/180)*math.sin(Crdao[i][j]*math.pi/180))/(math.cos(Decao[i]*math.pi/180)*math.cos(Crdao[i][j]*math.pi/180))
        if a > 1:
            a = 1
        Crrpao[i][j] = RAao[i] + (180/math.pi) * math.acos(a)
        Crrnao[i][j] = RAao[i] - (180/math.pi) * math.acos(a)

for i in range(9):
    Crdwo[i] = [Decwo[i]+1.0*math.sin(math.pi*j/999-0.5*math.pi) for j in range(1000)]
    for j in range(1000):
        a = (math.cos(1.0*math.pi/180)-math.sin(Decwo[i]*math.pi/180)*math.sin(Crdwo[i][j]*math.pi/180))/(math.cos(Decwo[i]*math.pi/180)*math.cos(Crdwo[i][j]*math.pi/180))
        if a > 1:
            a = 1
        Crrpwo[i][j] = RAwo[i] + (180/math.pi) * math.acos(a)
        Crrnwo[i][j] = RAwo[i] - (180/math.pi) * math.acos(a)

#Also draw a circle around the central 3 deg around the SEP:

C3d = [Decao[0]+3.0*math.sin(math.pi*j/999-0.5*math.pi) for j in range(1000)]
C3rp = [0]*1000
C3rn = [0]*1000
for j in range(1000):
    a = (math.cos(3.0*math.pi/180)-math.sin(Decao[0]*math.pi/180)*math.sin(C3d[j]*math.pi/180))/(math.cos(Decao[0]*math.pi/180)*math.cos(C3d[j]*math.pi/180))
    if a > 1:
        a = 1

    C3rp[j] = RAao[0] + (180/math.pi) * math.acos(a)
    C3rn[j] = RAao[0] - (180/math.pi) * math.acos(a)


col = ['#f58231', '#911eb4', '#f032e6']
fig, ax = plt.subplots(1,1, figsize=(10,6))
ax.plot(C3rp, C3d, linewidth=4, color='0.5', linestyle='dashed', label='3 deg SEP')
ax.plot(C3rn, C3d, linewidth=4, color='0.5', linestyle='dashed')
ax.plot(Crrpao[0], Crdao[0], linewidth=4, color='green', label='AAT already observed')
ax.plot(Crrpao[0], Crdao[0], linewidth=4, color='red', label='AAT will observe')
for i in range(7):
    ax.plot(Crrpao[i], Crdao[i], linewidth=4, color='black')
    ax.plot(Crrnao[i], Crdao[i], linewidth=4, color='black')
for i in range(9):
    ax.plot(Crrpwo[i], Crdwo[i], linewidth=4, color='black')
    ax.plot(Crrnwo[i], Crdwo[i], linewidth=4, color='black')
for i in range(7):
    ax.plot(Crrpao[i], Crdao[i], linewidth=2.5, color='green')
    ax.plot(Crrnao[i], Crdao[i], linewidth=2.5, color='green')
for i in range(9):
    ax.plot(Crrpwo[i], Crdwo[i], linewidth=2.5, color='red')
    ax.plot(Crrnwo[i], Crdwo[i], linewidth=2.5, color='red')

col = ['red', 'blue', 'orange', 'green', 'magenta', 'cyan', '0.5', 'yellow', 'purple']
for i in range(7, 16):
    ax.plot(vars()['RAo'+str(i)], vars()['Deco'+str(i)], '.', markersize=5, color=col[i-7], alpha=0.4, label='SEP'+str(i))
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_title('Variable sources identified and distributed between the different SEPfields')
plt.legend()
name = 'RaDec_Varsrc_AAT4.png'
plt.savefig(name, format='png')

#Now print out the results in the format needed for the fld files, and also save the result in a .txt file:

print('\n\n\n')

for i in range(7, 16):
    print('\n\nSEP'+str(i)+'\n')
    print('# \tR. Ascention \tDeclination \t\t\t\tProg \tProper Motion \tComments')
    print('# Name hh mm ss.sss \tdd mm ss.sss \t\t\tmag \tID \tra \tdec')
    for j in range(len(vars()['RAo'+str(i)])):
        print(str(vars()['Nam'+str(i)][j]) + '\t' + str(vars()['RAohms'+str(i)][j]) + '\t' + str(vars()['Decodms'+str(i)][j]) + '\tP\t' + str(vars()['Prty'+str(i)][j]) + '\t' + str(vars()['rmag'+str(i)][j]) + '\t0 \t0 \t0')


    F1 = open('VS_AATfldformat_SEP'+str(i)+'_4.txt', 'w+')
    F1.write('# \tR. Ascention \tDeclination \t\t\t\tProg \tProper Motion \tComments\n')
    F1.write('# Name hh mm ss.sss \tdd mm ss.sss \t\t\tmag \tID \tra \tdec\n')
    for j in range(len(vars()['RAo'+str(i)])):
        F1.write(str(vars()['Nam'+str(i)][j]) + '\t' + str(vars()['RAohms'+str(i)][j]) + '\t' + str(vars()['Decodms'+str(i)][j]) + '\tP\t' + str(vars()['Prty'+str(i)][j]) + '\t' + str(vars()['rmag'+str(i)][j]) + '\t0 \t0 \t0 \n')
    F1.close()



plt.show()
