#In this script, I will test whether the fld files are configured correctly. It will test whether the names are correctly defined, and whether the sources are correctly positioned.

#In version 1, I include the LMC region from Chandreyee.

#In version 2, I will also look at the distribution of the different categories of sources; P, F, and S.

#In version 3, I will do this for the new selection of sources for feb22. 

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from astropy.io import fits
from astropy.table import Table
import os.path
import matplotlib
import math
import numpy as np

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

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

#Read in the files in question:
rts = ['7a', '8a', '9a', '10', '11', '12', '13', '14', '15']

for i in range(len(rts)):
    vars()['Name'+rts[i]], vars()['RA'+rts[i]], vars()['Dec'+rts[i]], vars()['Prty'+rts[i]], vars()['rmag'+rts[i]], vars()['tp'+rts[i]] = [], [], [], [], [], []
    #Open the .txt file:
    F = open('SEP'+rts[i]+'fld_simp.txt', 'r')
    Fl = F.readlines()
    print(i)
    for x in Fl:
        vars()['Name'+rts[i]].append(x.split()[0])
        vars()['RA'+rts[i]].append(float(x.split()[1])*(360/24) + float(x.split()[2])*(360/(24*60)) + float(x.split()[3])*(360/(24*60*60)))
        vars()['Dec'+rts[i]].append(float(x.split()[4]) - float(x.split()[5])/60 - float(x.split()[6])/3600)
        vars()['tp'+rts[i]].append(x.split()[7])
        vars()['rmag'+rts[i]].append(float(x.split()[9]))
        vars()['Prty'+rts[i]].append(float(x.split()[8]))

    F.close()

    if os.path.isfile('SEP'+rts[i]+'fldasfits_feb22.fits') == 1:
        os.remove('SEP'+rts[i]+'fldasfits_feb22.fits')

    #Save to fits file:
    t = Table([vars()['RA'+rts[i]], vars()['Dec'+rts[i]], vars()['tp'+rts[i]], vars()['rmag'+rts[i]]], names=('RA', 'Dec', 'type', 'rmag'))
    t.write('SEP'+rts[i]+'fldasfits_feb22.fits', format='fits')

#Make a brief check that this all works by plotting the sources:

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
fig, ax = plt.subplots(1,1, figsize=(10,8))
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

ax.plot(C3rp, C3d, linewidth=4, color='0.5', linestyle='--', label='3 deg from SEP')
ax.plot(C3rn, C3d, linewidth=4, color='0.5', linestyle='--')
    
col = ['red', 'blue', 'orange', 'green', 'magenta', 'cyan', '0.5', 'yellow', 'purple']
for i in range(len(rts)):
    ax.plot(vars()['RA'+rts[i]][0], vars()['Dec'+rts[i]][0], '.', markersize=5, color=col[i], alpha=0.1*vars()['Prty'+rts[i]][0], label='SEP'+rts[i])
    for j in range(1, len(vars()['RA'+rts[i]])):
        ax.plot(vars()['RA'+rts[i]][j], vars()['Dec'+rts[i]][j], '.', markersize=5, color=col[i], alpha=0.1*vars()['Prty'+rts[i]][j])
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_title('AAT sources to be observed')
plt.legend()
name = 'RaDec_allsrc_AATfldfilesfinal2_SEP3.png'
plt.savefig(name, format='png')

#Now investigate the sources that lie outside of the field of these regions:
for i in range(len(rts)):
    for j in range(len(vars()['RA'+rts[i]])):
        if angbet([CobsRA[i], CobsDec[i]], [vars()['RA'+rts[i]][j], vars()['Dec'+rts[i]][j]]) > 60*60:
            print('In SEP'+str(i)+', source: ', vars()['Name'+rts[i]][j], ' lies outside the region that the field covers:', (angbet([CobsRA[i], CobsDec[i]], [vars()['RA'+rts[i]][j], vars()['Dec'+rts[i]][j]])/3600), vars()['RA'+rts[i]][j], vars()['Dec'+rts[i]][j])

#Now determine whether there are any duplicate sources in the lists.
for i in range(len(rts)):
    for j in range(1, len(vars()['RA'+rts[i]])):
        for k in range(j):
            if vars()['RA'+rts[i]][j] == vars()['RA'+rts[i]][k]:
                if vars()['Dec'+rts[i]][j] == vars()['Dec'+rts[i]][k]:
                    print('Duplicate source in SEP'+rts[i]+':')
                    print(vars()['Name'+rts[i]][j], vars()['Name'+rts[i]][k], vars()['RA'+rts[i]][j], vars()['RA'+rts[i]][k])

#Distinguish between the three different types: P, F, S
tc = ['P', 'F', 'S']
for r in range(len(rts)):
    for i in range(3):
        vars()['RA'+rts[r]+tc[i]] = []
        vars()['Dec'+rts[r]+tc[i]] = []
    for j in range(len(vars()['RA'+rts[r]])):
        t = -1
        for i in range(3):
            if vars()['tp'+rts[r]][j] == tc[i]:
                t = i
        if t == -1:
            #In this case, it is P_w2, and counts for P as well. 
            #print('Something went wrong with the type association:', vars()['tp'+str(r)][j])
            t = 0
        else:
            vars()['RA'+rts[r]+tc[t]].append(vars()['RA'+rts[r]][j])
            vars()['Dec'+rts[r]+tc[t]].append(vars()['Dec'+rts[r]][j])

#Now plot these separately, and then together, with the following colors:
colt = ['red', 'blue', 'green']
tn = ['Science targets', 'Guide stars', 'Sky fibres']

for u in range(3):
    fig, ax = plt.subplots(1,1, figsize=(10,8))
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

    ax.plot(C3rp, C3d, linewidth=4, color='0.5', linestyle='--', label='3 deg from SEP')
    ax.plot(C3rn, C3d, linewidth=4, color='0.5', linestyle='--')
    
    for i in range(len(rts)):
        ax.plot(vars()['RA'+rts[i]+tc[u]], vars()['Dec'+rts[i]+tc[u]], '.', markersize=5, color=col[i], label='SEP'+rts[i])
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.set_title('AAT sources to be observed, '+tn[u])
    plt.legend()
    name = 'RaDec_allsrc_AATfldfilesfinal2_SEP_tp'+str(tc[u])+'_3.png'
    plt.savefig(name, format='png')

#Determine the number of sources selected:
print('\nNumber of science targets, guide stars, and sky fibres selected:')
#print(np.sum(np.array([len(vars()['RA'+rts[i]+'P']) for i in range(len(rts))])), np.sum(np.array([len(vars()['RA'+rts[i]+'F']) for i in range(len(rts))])), np.sum(np.array([len(vars()['RA'+rts[i]+'S']) for i in range(len(rts))])))
print(len(RA7aP)+len(RA8aP)+len(RA9aP)+len(RA10P)+len(RA11P)+len(RA12P)+len(RA13P)+len(RA14P)+len(RA15P), len(RA7aF)+len(RA8aF)+len(RA9aF)+len(RA10F)+len(RA11F)+len(RA12F)+len(RA13F)+len(RA14F)+len(RA15F), len(RA7aS)+len(RA8aS)+len(RA9aS)+len(RA10S)+len(RA11S)+len(RA12S)+len(RA13S)+len(RA14S)+len(RA15S))


fig, ax = plt.subplots(1,1, figsize=(10,8))
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

ax.plot(C3rp, C3d, linewidth=4, color='0.5', linestyle='--', label='3 deg from SEP')
ax.plot(C3rn, C3d, linewidth=4, color='0.5', linestyle='--')
for u in range(3):
    ax.plot(vars()['RA7a'+tc[u]], vars()['Dec7a'+tc[u]], '.', markersize=5, color=colt[u], label=tn[u])
    for i in range(len(rts)):
        ax.plot(vars()['RA'+rts[i]+tc[u]], vars()['Dec'+rts[i]+tc[u]], '.', markersize=5, color=colt[u])
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_title('AAT sources to be observed, all source types')
plt.legend()
name = 'RaDec_allsrc_AATfldfilesfinal2_SEP_alltypes_3.png'
plt.savefig(name, format='png')


plt.show()
         
