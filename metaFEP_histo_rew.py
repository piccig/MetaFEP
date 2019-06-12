#!/usr/bin/env python

import numpy as np
from numpy import *
import itertools
from itertools import dropwhile

def is_comment(s):
    # return true if a line starts with #
    return s.startswith('#')

smin=-3.0
smax=+5.0
nbin=600
hist=np.zeros(nbin)
deltahist=np.zeros(nbin)
fep=np.zeros(nbin)

beta=1.0/2.519282

input_filename='COLVAR'
cv_col=1
rbias_col=2
u1_col=3
u2_col=4

fcol=open(input_filename,'r')
ffep=open('MetaFEP.dat','w')

data= np.loadtxt(dropwhile(is_comment,fcol),usecols=(u1_col,u2_col))

u1_min = np.mean(data[:,0])
u2_min = np.mean(data[:,1])


fcol.seek(0, 0)

for line in dropwhile(is_comment,fcol):
#for line in fcol:
    line=line.strip()
    columns=line.split()
    scv=float(columns[cv_col])
    binning=int(((scv-smin)/(smax-smin))*nbin)
    rbias=float(columns[rbias_col])
    u2=float(columns[u2_col])-u2_min
    u1=float(columns[u1_col])-u1_min
    hist[binning]+=np.exp(+beta*rbias)
    deltahist[binning]+=np.exp(+beta*(rbias-u2+u1))


indexes=np.zeros(nbin)
data_delta=np.zeros(nbin)

for i in range(nbin):
    
    if(deltahist[i]==0.0):
      delta_norm=0.0
    else:
      delta_norm=deltahist[i]/hist[i]

    fep[i]=float(-np.log(hist[i])/beta)+float(-np.log(delta_norm)/beta)


where_are_NaNs = isnan(fep)
fep[where_are_NaNs] = 350.0

Nmov=20
txt='# CV   MetaFEP(moving average)  95% confindence(lower bound)  95% confindence(upper bound) MetaFEP(raw)'+'\n'
ffep.write(txt)
for i in range(0,nbin-Nmov):

    fep_mov=np.mean(fep[i:i+Nmov])
    fep_std=np.std(fep[i:i+Nmov])


    txt=str((i+int(Nmov/2.0))*(smax-smin)/nbin+smin)+'\t'+str(fep_mov)+'\t'+str(fep_mov-1.9600*fep_std/np.sqrt(Nmov))+'\t'+str(fep_mov+1.9600*fep_std/np.sqrt(Nmov))+'\t'+str(fep[i+int(Nmov/2)])+'\n'

    ffep.write(txt)

quit()

