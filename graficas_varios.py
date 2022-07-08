import matplotlib.pyplot as plt
import numpy as np
import math as mt
from scipy.special import eval_legendre
from scipy import integrate
import pandas as pd

SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 14
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
x = np.linspace(-1, 1, num=20001)
ang=np.arccos(x)
grad=ang*180/(mt.pi)

g1120=np.loadtxt('100120gmcgb.txt',usecols=2,skiprows=1,delimiter=', ')
g19120=np.loadtxt('190120gmcgb.txt',usecols=2,skiprows=1,delimiter=', ')
g22120=np.loadtxt('220120gmcgb.txt',usecols=2,skiprows=1,delimiter=', ')

b1120=np.loadtxt('100120gmcgb.txt',usecols=3,skiprows=1,delimiter=', ')
b130120=np.loadtxt('130120gmcgb.txt',usecols=3,skiprows=1,delimiter=', ')
b160120=np.loadtxt('160120gmcgb.txt',usecols=3,skiprows=1,delimiter=', ')
b19120=np.loadtxt('190120gmcgb.txt',usecols=3,skiprows=1,delimiter=', ')
b22120=np.loadtxt('220120gmcgb.txt',usecols=3,skiprows=1,delimiter=', ')

plt.plot(grad,b1120,'-b',label='s=1.0')
plt.plot(grad,b130120,'-c',label='s=1.30')
plt.plot(grad,b160120,'-m',label='s=1.60')
plt.plot(grad,b19120,'-g',label='s=1.90')
plt.plot(grad,b22120,'-r',label='s=2.20')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$b(\theta)$')
plt.title('n=120')
plt.legend()
#plt.title('n='+t0+ ', s='+str(s))
plt.show()
#plt.savefig(t0+'g.png')
plt.clf()
