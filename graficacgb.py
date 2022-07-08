from matplotlib import pyplot as plt
import numpy as np
import math as mt
from scipy.special import eval_legendre
from scipy import integrate
import pandas as pd
import csv
import os

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


dir_name='C:\\Users\52333\Documents\liquidos\programas_python\facil_mezclacte\gfig'
#plt.rcParams["savefig.directory"] = os.chdir(os.path.dirname(dir_name))

x = np.linspace(-1, 1, num=20001)
ang=np.arccos(x)
grad=ang*180/(mt.pi)
c_dict={}
g_dict={}
b_dict={}

#es el n anterior con el que generamos los datos para gamma en oz.py
n0=20
au=20
#este es el n al que queremos llegar mas una iteración (queremos 120)
nf=120+au


fin=int((nf-n0)/au)

aus=0.05
si=1.0
sf=2.6+aus
finsc=int((sf-si)/aus)
s=si
n0=20


for nc in range(0,finsc):
  for i in range(0,fin):
      name=round(s*100)
      name2=round(s,2)
      sstr=str(name)+str(n0)
      sstrc='cfig/'+str(name)+'_'+str(n0)
      sstrg='gfig/'+str(name)+'_'+str(n0)
      sstrb='bfig/'+str(name)+'_'+str(n0)
      if(i<=finsc-1):
          c_dict[i]=np.loadtxt(sstr+'gmcgb.txt',skiprows=1,usecols=1,delimiter=',  ')
          g_dict[i]=np.loadtxt(sstr+'gmcgb.txt',skiprows=1,usecols=2,delimiter=',  ')
          b_dict[i]=np.loadtxt(sstr+'gmcgb.txt',skiprows=1,usecols=3,delimiter=',  ')
      if(i==finsc):
          c_dict[finsc]=c_dict[finsc-1]
          g_dict[finsc]=g_dict[finsc-1]
          b_dict[finsc]=b_dict[finsc-1]
      c=c_dict[i]
      g=g_dict[i]
      b=b_dict[i]
      plt.plot(grad,g,'-',markersize=5,label=s)
      plt.xlabel(r'$\theta$')
      plt.ylabel(r'$g(\theta)$')
      plt.title('n='+str(n0)+ ', s='+str(name2))
      plt.savefig(sstrg+'g.png')
      plt.clf()

      plt.plot(grad,c,'-',markersize=5,label=s)
      plt.title('n='+str(n0)+ ', s='+str(name2))
      plt.xlabel(r'$\theta$')
      plt.ylabel(r'$c(\theta)$')
      plt.savefig(sstrc+'c.png')
      plt.clf()

      plt.plot(grad,b,'-',markersize=5,label=s)
      plt.title('n='+str(n0)+ ', s='+str(name2))
      plt.xlabel(r'$\theta$')
      plt.ylabel(r'$b(\theta)$')
      plt.savefig(sstrb+'b.png')
      plt.clf()
      n0=n0+au

  s=s+aus
  n0=20
  g_dict.clear()
  c_dict.clear()
  b_dict.clear()


g=np.loadtxt("14040gmcgb.txt",skiprows=1,usecols=2,delimiter=',  ')
gant=np.loadtxt("14040gant.txt",skiprows=1,usecols=1,delimiter=',  ')
gsig=np.loadtxt("14040gsig.txt",skiprows=1,usecols=1,delimiter=',  ')

g2=np.loadtxt("13540gmcgb.txt",skiprows=1,usecols=2,delimiter=',  ')
gant2=np.loadtxt("13540gant.txt",skiprows=1,usecols=1,delimiter=',  ')
gsig2=np.loadtxt("13540gsig.txt",skiprows=1,usecols=1,delimiter=',  ')

plt.plot(grad,g,label='g')
plt.plot(grad,gant,label='gant')
plt.plot(grad,gsig,label='gsig')

plt.plot(grad,g2,label='g2')
plt.plot(grad,gant2,label='gant2')
plt.plot(grad,gsig2,label='gsig2')
plt.legend()
plt.show()

'''
#g120=np.loadtxt('120gant1.txt',skiprows=1,usecols=0)
g=np.loadtxt('20gant70.txt',skiprows=1,usecols=0)
gsig=np.loadtxt('7020gsig.txt',skiprows=1,usecols=1)
gant=np.loadtxt('7020gant.txt',skiprows=1,usecols=1)
g0=np.loadtxt('7020gmcgb.txt',skiprows=1,usecols=2,delimiter=',  ')
plt.plot(grad,g)
plt.plot(grad,gsig)
plt.plot(grad,gant)
plt.plot(grad,g0)
plt.show()
'''
