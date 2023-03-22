#In this file, I will sort out the brightest stars that are not high proper motion.
#I have downloaded the relevant data from simbad into text files.

#In version 1 I will get this to work for the feb22 dataset. 

import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title('ABC123 vs $\mathrm{ABC123}^{123}$')

#Open the file in question, and select the columns of importantce:
ID, Nam, Tp, RA, Dec, pmRA, pmDec, rmag = [], [], [], [], [], [], [], []
tst = []

#F = open('SIMBAD_LMC.txt', 'r')
F = open('SIMBAD_SEP7a.txt', 'r')
Fl = F.readlines()
for x in Fl:
    ID.append(x.split()[0])
    Nam.append(x.split()[2])
    Tp.append(x.split()[3])
    if float(x.split()[6]) < 10:
        RA.append(x.split()[4] + ' ' + x.split()[5] + ' 0' + str(round(float(x.split()[6]), 2)))
    else:
        RA.append(x.split()[4] + ' ' + x.split()[5] + ' ' + str(round(float(x.split()[6]), 2)))
    if float(x.split()[9]) < 10:
        Dec.append(x.split()[7] + ' ' + x.split()[8] + ' 0' + str(round(float(x.split()[9]), 2)))
    else:
        Dec.append(x.split()[7] + ' ' + x.split()[8] + ' ' + str(round(float(x.split()[9]), 2)))
    pmRA.append(x.split()[10])
    pmDec.append(x.split()[11])
    rmag.append(x.split()[12])
    tst.append(x.split()[4])

N = len(Tp)

#Because so many sources have complicated names, I will instead get rid of all sources that have '~' in their 10, 11, or 12 column.

print('Length of array at the start:', len(Tp), len(ID), len(Nam), len(RA), len(Dec), len(pmRA), len(pmDec), len(rmag), len(tst))

k = 0
while k < len(Tp):
    if pmRA[k] == '~':
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    elif pmDec[k] == '~':
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    elif rmag[k] == '~':
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    else:
        k += 1

N = len(Tp)
print('Length of array after deleting initially problematic/ not useful entries:', len(Tp))

for i in range(N):
    if tst[i] != '05':
        if tst[i] != '06':
            #there must be an error:
            print('Line', ID[i], 'has an error, probably because the name contains a space')

#Now remove the sources that have an rmag of greater than 13.5. And a proper motion of more than 20, taking the squashing of lines of constant declination into account. So pm of 1 in ra is equal to a proper motion of 0.3977 in dec. 


print('\nLength before deleting sources that are too dim, or have too large proper motions:', len(ID))
k = 0
while k < len(Tp):
    pm = math.sqrt((0.3977*float(pmRA[k]))**2 + float(pmDec[k])**2)
    if pmRA[k] == '~':
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    elif pmDec[k] == '~':
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    elif rmag[k] == '~':
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    #elif pm > 20:
    elif pm > 100: # relax condition from what I had before. 
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    #elif float(rmag[k]) > 13.5:
    #For the LMC region, take the lower limit as 12.0
    #elif float(rmag[k]) > 12.0:
    #since we want more sources, take the limit of 13.0, and 14.0
    elif float(rmag[k]) > 14.0:
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    elif float(rmag[k]) < 12.0: # remove this if there are too few sources. 
        del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    #elif float(rmag[k]) < 10.0: # remove this if there are too few sources. 
    #    del ID[k], Nam[k], Tp[k], RA[k], Dec[k], pmRA[k], pmDec[k], rmag[k], tst[k]
    else:
        k += 1
print('Length after deleting sources that are too dim, too bright, or have too large proper motions:', len(ID))
        

rmag1 = [float(rmag[i]) for i in range(len(rmag))]
#Now sort the resulting sources by rmag, and then print them out:
#This is where there was a problem, due to multiple sources having the same rmag
iin0 = [i for i in range(len(rmag1))]
iin = [x for _,x in sorted(zip(rmag1, iin0))]
#ID1 = [x for _,x in sorted(zip(rmag1, ID))]
#Nam1 = [x for _,x in sorted(zip(rmag1, Nam))]
#Tp1 = [x for _,x in sorted(zip(rmag1, Tp))]
#RA1 = [x for _,x in sorted(zip(rmag1, RA))]
#Dec1 = [x for _,x in sorted(zip(rmag1, Dec))]
#pmRA1 = [x for _,x in sorted(zip(rmag1, pmRA))]
#pmDec1 = [x for _,x in sorted(zip(rmag1, pmDec))]
#rmag = [x for _,x in sorted(zip(rmag1, rmag1))]
#print(rmag1)
ID1 = [ID[iin[j]] for j in range(len(rmag1))]
Nam1 = [Nam[iin[j]] for j in range(len(rmag1))]
Tp1 = [Tp[iin[j]] for j in range(len(rmag1))]
RA1 = [RA[iin[j]] for j in range(len(rmag1))]
Dec1 = [Dec[iin[j]] for j in range(len(rmag1))]
pmRA1 = [pmRA[iin[j]] for j in range(len(rmag1))]
pmDec1 = [pmDec[iin[j]] for j in range(len(rmag1))]
rmag = [rmag1[iin[j]] for j in range(len(rmag1))]



print('As a preliminary check:')
#for i in range(len(ID)):
#    print(str(Nam1[i]) + '\t' + str(Tp1[i]) + '\t+' + str(RA1[i]) + '\t' + str(Dec1[i]) + '\tF\t9\t' + str(rmag[i]) + '\t0\t' + str(round(float(pmRA1[i])/1000 ,6)) + '\t' + str(round(float(pmDec1[i])/1000 ,6)))

#But I think I will only take sources that are classififed as stars, otherwise that could be difficult. So delete sources that are not labeled as "*"

print('\n\n')

print('\nLength before deleting sources that are not stars:', len(ID1))

k = 0
while k < len(ID1):
    #if Tp1[k] != '*':
    #If there are not many sources there,just delete the high PM stars:
    #if Tp1[k] == 'PM*':
    #No, alsways run this, there should always just be stars there.
    a = 0
    if Tp1[k] == '*':
        a = 1
    elif Tp1[k] == 'HB*':
        a = 1
    elif Tp1[k] == 'Ae*':
        a = 1
    elif Tp1[k] == 'Em*':
        a = 1
    elif Tp1[k] == 'Be*':
        a = 1
    elif Tp1[k] == 'BS*':
        a = 1
    elif Tp1[k] == 'RG*':
        a = 1
    elif Tp1[k] == 'AB*':
        a = 1
    elif Tp1[k] == 'C*':
        a = 1
    elif Tp1[k] == 'S*':
        a = 1
    elif Tp1[k] == 'sg*':
        a = 1
    elif Tp1[k] == 's*r':
        a = 1
    elif Tp1[k] == 's*y':
        a = 1
    elif Tp1[k] == 'HS*':
        a = 1
#    #For fields that have very few stars, also include these:
#    elif Tp1[k] == 'LP*':
#        a = 1
#    elif Tp1[k] == 'RR*':
#        a = 1
#    elif Tp1[k] == 'EB*':
#        a = 1
#    elif Tp1[k] == 'WD*':
#        a = 1
    if a == 0:
        del ID1[k], Nam1[k], Tp1[k], RA1[k], Dec1[k], pmRA1[k], pmDec1[k], rmag[k]
    else:
        k += 1
        
print('\nLength after deleting sources that are not stars:', len(ID1))

print('\n\nAnd now, the final output, to be put into the fld file:')
for i in range(len(ID1)):
    print(str(Nam1[i]) + '\t+' + str(RA1[i]) + '\t' + str(Dec1[i]) + '\tF\t9\t' + str(rmag[i]) + '\t0\t' + str(round(float(pmRA1[i])/1000 ,6)) + '\t' + str(round(float(pmDec1[i])/1000 ,6)))


