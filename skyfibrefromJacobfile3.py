#In this script, I will print out the fld format for the sky fibres selected by Jacob in his file.

#In version 1 of this script, I will also plot the LMC region

#In 1a, I will just select some other positions randomly to fill out the entire LMC field. I have given up trying to solve the issues with the code. I just want the result now.

#In version 3, I will generate the sky fibre positions for the feb22 run, and I will see which fields need some assistance. In various versions of this, I will generate the missing regions as well. 


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

F = fits.open('sky.sep_feb22.fit')
RA0 = F[1].data['RA']
Dec0 = F[1].data['Dec']
is100_0 = F[1].data['in_sample_100']
is500_0 = F[1].data['in_sample_500']
reg0 = F[1].data['region']
igr0 = F[1].data['in_good_region']
gc0 = F[1].data['gcoverage']
rc0 = F[1].data['rcoverage']
ic0 = F[1].data['icoverage']
zc0 = F[1].data['zcoverage']
yc0 = F[1].data['ycoverage']
F.close()

#Only take sources that have is100_0 == 1

RA, Dec, is100, reg = [], [], [], []
for i in range(len(RA0)):
    if is100_0[i] == 1:
        RA.append(RA0[i])
        Dec.append(Dec0[i])
        is100.append(is100_0[i])
        reg.append(reg0[i])

print(len(reg), reg[0], reg[-1])
        

#Since the above doesn't work, I need to try something else:
#for i in range(len(RA0)):
#    if igr0[i] == 1:
#        if gc0[i] * rc0[i] * ic0[i] * zc0[i] * yc0[i] == 1:
#            RA.append(RA0[i])
#            Dec.append(Dec0[i])
#            reg.append(reg0[i])
        
#print(reg)

for i in range(9):
    vars()['RA'+str(i+7)] = []
    vars()['RAhms'+str(i+7)] = []
    vars()['Dec'+str(i+7)] = []
    vars()['Decdms'+str(i+7)] = []

for i in range(len(RA)):
    for j in range(9):
        #if reg[i] == 'region_'+str(j):
        if reg[i] == j:
            vars()['RA'+str(j+7)].append(RA[i])
            vars()['RAhms'+str(j+7)].append(radegtohms(RA[i]))
            vars()['Dec'+str(j+7)].append(Dec[i])
            vars()['Decdms'+str(j+7)].append(decdegtodms(Dec[i]))


#Now write the output in the format of the fld file:
for i in range(7, 16):
    print('\n\nSEP'+str(i)+'\n')
    for j in range(len(vars()['RA'+str(i)])):
        print('skyfibre'+str(j+1) + '\t' + str(vars()['RAhms'+str(i)][j]) + '\t' + str(vars()['Decdms'+str(i)][j]) + '\tS\t9\t30.0\t0\t0\t0')

##Because the code doesn't work as I want it to, select some random sources, and select them as sky fibre positions.
#First, do SEP9:
#I need 20, so generate 60 here, to be sure.

print('\n\n')
S9rata, S9decta = [0]*60, [0]*60
S9ratah, S9dectad = [0]*60, [0]*60
for i in range(60):
    a, t = 0, 0
    while a == 0: #no selection has been made.
        print('trial', t, 'for additional sky field source ', i)
        t += 1
        j = np.random.randint(0, len(RA0))
        if reg0[j] == 2: #SEP reg 9 is 2 here.
            if igr0[j] == 1:
                if RA0[j] > 88.5:
                    if Dec0[j] > -67.6:
                        #Check that there is no other selected source within 5 arcminutes.
                        mina = 5*60
                        sma = 0
                        for k in range(len(S9rata)):
                            if angbet([S9rata[k], S9decta[k]], [RA0[j], Dec0[j]]) < mina:
                                sma = 1
                        if sma == 0:
                            S9rata[i] = RA0[j]
                            S9decta[i] = Dec0[j]
                            S9ratah[i] = radegtohms(RA0[j])
                            S9dectad[i] = decdegtodms(Dec0[j])
                            a = 1

#Now print out the extra regions to take:

print('\nExtra sources to select from for SEP9:')
for i in range(60):
    print(S9ratah[i], S9dectad[i])

#And here are the ones that I will add:
RAhms9add = ['+06 00 08.442', '+05 55 34.225', '+05 54 13.826', '+06 01 07.516', '+06 00 49.783', '+05 54 13.725', '+06 00 13.688', '+05 58 12.976', '+05 55 42.849', '+05 56 22.195', '+05 57 47.219', '+05 55 39.792', '+05 58 57.065', '+05 55 01.809', '+05 58 27.758', '+06 00 28.65', '+06 00 37.278', '+05 57 29.061', '+05 56 13.551', '+05 57 29.252']
RA9add = [float(RAhms9add[i].split()[0])*360/24 + float(RAhms9add[i].split()[1])*360/(24*60) + float(RAhms9add[i].split()[2])*360/(24*60*60) for i in range(len(RAhms9add))]
Decdms9add = ['-67 20 58.346', '-67 11 01.731', '-67 19 02.43', '-67 11 51.372', '-66 57 38.07', '-67 29 34.531', '-67 06 02.942', '-66 54 23.291', '-66 15 18.979', '-67 15 39.902', '-67 29 51.512', '-67 02 23.155', '-66 38 05.131', '-66 43 57.72', '-67 16 45.764', '-66 44 32.833', '-67 32 58.905', '-67 20 54.835', '-66 56 26.735', '-66 24 37.323']
Dec9add = [float(Decdms9add[i].split()[0]) - float(Decdms9add[i].split()[1])/60 + float(Decdms9add[i].split()[2])/3600 for i in range(len(Decdms9add))]

for i in range(len(RA9add)):
    RA9.append(RA9add[i])
    Dec9.append(Dec9add[i])
    RAhms9.append(RAhms9add[i])
    Decdms9.append(Decdms9add[i])

#Now repeat this for SEP12: 
S12rata, S12decta = [0]*20, [0]*20
S12ratah, S12dectad = [0]*20, [0]*20
for i in range(20):
    a, t = 0, 0
    while a == 0: #no selection has been made.
        print('trial', t, 'for additional sky field source ', i)
        t += 1
        j = np.random.randint(0, len(RA0))
        if reg0[j] == 5: #SEP reg 12 is 5 here.
            if igr0[j] == 1:
                if RA0[j] > 93:
                    if Dec0[j] > -64.8:
                        #Check that there is no other selected source within 3 arcminutes.
                        mina = 3*60
                        sma = 0
                        for k in range(len(S12rata)):
                            if angbet([S12rata[k], S12decta[k]], [RA0[j], Dec0[j]]) < mina:
                                sma = 1
                        if sma == 0:
                            S12rata[i] = RA0[j]
                            S12decta[i] = Dec0[j]
                            S12ratah[i] = radegtohms(RA0[j])
                            S12dectad[i] = decdegtodms(Dec0[j])
                            a = 1

#Now print out the extra regions to take:

print('\nExtra sources to select from for SEP12:')
for i in range(20):
    print(S12ratah[i], S12dectad[i])

#And here are the ones that I will add:
RAhms12add = ['+06 12 27.734', '+06 15 13.714', '+06 14 08.708', '+06 13 53.761', '+06 13 23.346']
RA12add = [float(RAhms12add[i].split()[0])*360/24 + float(RAhms12add[i].split()[1])*360/(24*60) + float(RAhms12add[i].split()[2])*360/(24*60*60) for i in range(len(RAhms12add))]
Decdms12add = ['-64 23 48.644', '-64 47 32.978', '-64 27 45.87', '-64 38 48.868', '-64 30 57.904']
Dec12add = [float(Decdms12add[i].split()[0]) - float(Decdms12add[i].split()[1])/60 + float(Decdms12add[i].split()[2])/3600 for i in range(len(Decdms12add))]

for i in range(len(RA12add)):
    RA12.append(RA12add[i])
    Dec12.append(Dec12add[i])
    RAhms12.append(RAhms12add[i])
    Decdms12.append(Decdms12add[i])

#Now repeat for SEP10:
#Now repeat this for SEP12: 
S10rata, S10decta = [0]*200, [0]*200
S10ratah, S10dectad = [0]*200, [0]*200
for i in range(150):
    a, t = 0, 0
    while a == 0: #no selection has been made.
        print('trial', t, 'for additional sky field source ', i)
        t += 1
        j = np.random.randint(0, len(RA0))
        if reg0[j] == 3: #SEP reg 10 is 3 here.
            if igr0[j] == 1:
                if Dec0[j] > -0.61111*RA0[j] -12.57:
                    if Dec0[j] < -0.588235*RA0[j] -12.01:
                        #Check that there is no other selected source within 3 arcminutes.
                        mina = 3*60
                        sma = 0
                        for k in range(len(S10rata)):
                            if angbet([S10rata[k], S10decta[k]], [RA0[j], Dec0[j]]) < mina:
                                sma = 1
                        if sma == 0:
                            S10rata[i] = RA0[j]
                            S10decta[i] = Dec0[j]
                            S10ratah[i] = radegtohms(RA0[j])
                            S10dectad[i] = decdegtodms(Dec0[j])
                            a = 1

#Now print out the extra regions to take:

print('\nExtra sources to select from for SEP10:')
for i in range(150):
    print(S10ratah[i], S10dectad[i])

#And here are the ones that I will add:
RAhms10add = ['+06 14 55.921', '+06 13 09.646', '+06 13 21.459', '+06 13 05.711', '+06 15 11.573', '+06 07 33.06', '+06 10 10.249', '+06 12 39.987', '+06 14 42.939', '+06 11 11.81', '+06 13 28.314', '+06 15 39.583', '+06 09 59.298', '+06 09 23.283', '+06 12 42.18', '+06 06 42.675', '+06 20 07.842', '+06 08 22.955', '+06 07 38.041', '+06 18 12.748', '+06 12 18.471', '+06 10 31.304', '+06 20 33.737', '+06 10 50.006', '+06 14 31.408', '+06 16 11.514', '+06 07 58.808', '+06 14 12.212', '+06 09 09.545', '+06 15 48.923', '+06 11 30.064', '+06 04 43.509', '+06 09 41.532', '+06 11 09.033', '+06 17 56.431', '+06 14 10.072', '+06 14 04.159', '+06 20 05.274', '+06 06 37.199', '+06 17 07.431', '+06 12 34.738', '+06 08 38.636', '+06 10 28.292', '+06 15 45.303', '+06 12 57.744', '+06 16 09.582', '+06 18 38.6', '+06 19 39.232', '+06 17 33.44', '+06 08 50.061']
RA10add = [float(RAhms10add[i].split()[0])*360/24 + float(RAhms10add[i].split()[1])*360/(24*60) + float(RAhms10add[i].split()[2])*360/(24*60*60) for i in range(len(RAhms10add))]
Decdms10add = ['-67 55 18.804', '-68 11 22.826', '-68 33 55.598', '-67 41 36.84', '-67 46 15.58', '-68 13 33.872', '-68 02 53.451', '-67 30 26.592', '-68 31 26.514', '-67 20 10.654', '-67 26 16.532', '-68 25 41.645', '-68 22 02.745', '-67 41 49.813', '-68 56 45.719', '-68 18 30.304', '-68 35 26.953', '-68 21 47.164', '-67 36 04.905', '-68 32 56.944', '-68 59 07.11', '-68 44 17.635', '-68 49 03.828', '-68 58 35.603', '-67 34 05.974', '-67 51 41.617', '-67 47 19.443', '-67 47 25.635', '-67 21 41.036', '-68 35 13.679', '-67 45 00.642', '-68 05 48.987', '-68 40 15.401', '-68 09 26.62', '-68 16 20.411', '-68 30 33.016', '-68 02 15.592', '-68 21 11.265', '-68 13 54.933', '-67 46 37.585', '-68 43 00.418', '-67 20 40.65', '-67 29 16.482', '-68 52 11.353', '-67 54 06.838', '-67 32 22.289', '-68 15 37.124', '-68 07 02.056', '-68 43 28.779', '-67 29 37.901']
Dec10add = [float(Decdms10add[i].split()[0]) - float(Decdms10add[i].split()[1])/60 + float(Decdms10add[i].split()[2])/3600 for i in range(len(Decdms10add))]

for i in range(len(RA10add)):
    RA10.append(RA10add[i])
    Dec10.append(Dec10add[i])
    RAhms10.append(RAhms10add[i])
    Decdms10.append(Decdms10add[i])

#Now print out all the data again, so that I can add it to the files:
for i in range(7, 16):
    print('\n\nSEP'+str(i)+'\n')
    for j in range(len(vars()['RA'+str(i)])):
        print('skyfibre'+str(j+1) + '\t' + str(vars()['RAhms'+str(i)][j]) + '\t' + str(vars()['Decdms'+str(i)][j]) + '\tS\t9\t30.0\t0\t0\t0')


#for i in range(50): #Select 50 more.
#    a, t = 0, 0
#    while a == 0: #no selection has been made
#        print('trial', t, 'for additional sky field source ', i)
#        t += 1
#        j = np.random.randint(0, len(RA0))
#        if angbet([RA0[j], Dec0[j]], [90, -66.560708333333]) > 3 * 60 * 60: #Need to make sure that this source is outside the 3deg region, for which we already have sufficient sources.
#            print('Source is inside of 3degSEP')
#            if igr0[j] == 1:
#                #if ils0[j] == 0:
#                if ils0[j] > -2: # Just a useless command that is always true. 
#                    print('Source is in good region and not in LMC/SMC')
#                    #Now check whether there are any nearby sources:
#                    ns = 0
#                    for k in range(len(RA10)):
#                        if angbet([RA0[j], Dec0[j]], [RA10[k], Dec10[k]]) < 2 * 60: # choose the minimum separation to be 5 arcminutes
#                            ns += 1
#                            break
#                    if ns == 0:
#                        RA10.append(RA0[j])
#                        RAhms10.append(radegtohms(RA0[j]))
#                        Dec10.append(Dec0[j])
#                        Decdms10.append(decdegtodms(Dec0[j]))
#                        a = 1
#        elif Dec0[j] < -67.9: # this is another region that lacks any sky fibres. 
#            if igr0[j] == 1:
#                #if ils0[j] == 0:
#                if ils0[j] > -2: # Just a useless command that is always true. 
#                    #Now check whether there are any nearby sources:
#                    ns = 0
#                    for k in range(len(RA10)):
#                        if angbet([RA0[j], Dec0[j]], [RA10[k], Dec10[k]]) < 5 * 60: # choose the minimum separation to be 5 arcminutes
#                            ns += 1
#                            break
#                    if ns == 0:
#                        RA10.append(RA0[j])
#                        RAhms10.append(radegtohms(RA0[j]))
#                        Dec10.append(Dec0[j])
#                        Decdms10.append(decdegtodms(Dec0[j]))
#                        a = 1

#Check that my selection based on the above actually works:
#RAhms10add = ['+05 27 08.712', '+05 25 18.534', '+05 26 39.697', '+05 38 51.921', '+05 25 52.158', '+05 29 29.428', '+05 28 57.181', '+05 30 37.17', '+05 27 26.44', '+05 25 37.215', '+05 31 53.388', '+05 31 29.07', '+05 28 24.624', '+05 26 31.168', '+05 25 17.272', '+05 30 37.53', '+05 27 30.536', '+05 29 28.898', '+05 25 46.234', '+05 34 04.096', '+05 31 14.681', '+05 29 25.845', '+05 33 36.276', '+05 28 42.971', '+05 29 25.228', '+05 34 43.414', '+05 28 08.718', '+05 26 31.096', '+05 33 26.301', '+05 29 27.676' ]
#RA10add = [float(RAhms10add[i].split()[0])*360/24 + float(RAhms10add[i].split()[1])*360/(24*60) + float(RAhms10add[i].split()[2])*360/(24*60*60) for i in range(len(RAhms10add))]
#Decdms10add = ['-67 03 55.721', '-67 01 37.14', '-67 35 03.525', '-68 06 21.859', '-66 48 12.659', '-67 55 55.119', '-67 50 04.977', '-67 32 40.37', '-67 49 33.761', '-67 29 06.474', '-68 01 43.574', '-68 01 59.536', '-67 09 34.648', '-66 44 46.024', '-67 13 30.459', '-67 53 09.573', '-66 48 39.153', '-67 22 26.644', '-67 43 01.871', '-68 12 05.09', '-67 53 29.368', '-67 06 19.512', '-68 11 09.305', '-66 35 08.217', '-67 43 50.222', '-68 02 34.876', '-67 27 21.641', '-67 00 16.397', '-67 54 40.264', '-66 49 24.808']
#Dec10add = [float(Decdms10add[i].split()[0]) - float(Decdms10add[i].split()[1])/60 + float(Decdms10add[i].split()[2])/3600 for i in range(len(Decdms10add))]

#for i in range(len(RA10add)):
#    RA10.append(RA10add[i])
#    Dec10.append(Dec10add[i])
#    RAhms10.append(RAhms10add[i])
#    Decdms10.append(Decdms10add[i])

#Plot the positions:
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
    ax.plot(vars()['RA'+str(i)], vars()['Dec'+str(i)], '.', markersize=5, color=col[i-7], alpha=0.4, label='SEP'+str(i))
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_title('AAT sources to be observed')
plt.legend()
name = 'RaDec_skyfibresSEPLMC_in100_3.png'
plt.savefig(name, format='png')




plt.show()
