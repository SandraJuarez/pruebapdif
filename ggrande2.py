import matplotlib.pyplot as plt
import numpy as np
import math as mt
from scipy.special import eval_legendre
from scipy import integrate
import pandas as pd
import csv



#vamos a utilizar los resultados de gaf con n=30 que se calculó en oz.py
#cargamos la gamma obtenida para cada valor de s.
#el nombre del achivo es s*100gaf.txt
gam_dict={}
num=100
nf=235
au1=5
fin1=int((nf-num)/au1)
for i in range(0,fin1):
    file_name=str(num)+'gam.txt'
    gam_dict[i]=np.loadtxt(file_name,skiprows=1)
    num=num+5

gam_dict[fin1]=gam_dict[fin1-1]#np.loadtxt('200gam.txt',skiprows=1)
print(gam_dict[fin1])


#cargamos los parámetros que utilizamos en  oz.py
#bt=pd.read_csv("parametros.csv",usecols=['bt'],delimiter=None,header='infer')
bt=np.loadtxt("parametros.csv",skiprows=1,usecols=0,delimiter=',')
#radio de la partícula
a=np.loadtxt("parametros.csv",skiprows=1,usecols=1,delimiter=',')
#para el potencial
k=np.loadtxt("parametros.csv",skiprows=1,usecols=2,delimiter=',')
b=np.loadtxt("parametros.csv",skiprows=1,usecols=3,delimiter=',')
al=np.loadtxt("parametros.csv",skiprows=1,usecols=4,delimiter=',')

#parametros=pd.read_csv('parametros.csv')


x = np.linspace(-1, 1, num=20001)
ang=np.arccos(x)
grad=ang*180/(mt.pi) #cambiamos de radianes a grados
r=2*a*np.sin(ang/2)

grad=ang*180/(mt.pi) #cambiamos de radianes a grados
np.savetxt('grad.txt', np.transpose([grad]))

#u=(500*a*np.exp(-0.15*(r-a)))/r

#utilizaremos los valores del potencial que calculamos en oz.py
u=np.loadtxt("potencial.txt",skiprows=1,usecols=0)

#u=2*np.exp(-2*r)/r
#plt.plot(r,u)
#print(u)
sum=0.0



#es el n anterior con el que generamos los datos para gamma en oz.py
n0=20
#este es el n al que queremos llegar mas una iteración (queremos 120)
nf=125
#queremos aumentar n  en pasos de 10
au=5
fin=int((nf-n0)/au)
aus=0.05
si=1.0
sf=2.4+aus
finsc=int((sf-si)/aus)
s=si

'''
for sc in range(0,finsc):
    gaf=gam_dict[sc] #gamma funcion
    n0=20
    for nc in range(0,fin):
        for w in range(0,750):
        #calculamos c(x) con la primera aproximación para gamma(x)
            cfh=np.exp(gaf-bt*u)-1-gaf
            cfp=np.exp(-bt*u)*(gaf+1)-1-gaf
            cf=s*cfh+(1-s)*cfp
        #print(cf,w)
        #guardamos el valor anterior de gama(x)
            anterior=gaf
        #ciclo para calcular las transformadas de Legendre
            for m in range(0,150):
                #vamos a calcular los coeficientes cm
                pol=(eval_legendre(m, x))
                y = pol*cf
                cm=(2*m+1)*(integrate.simpson(y, x))/2
                #print(cm,w)
                #y con esos coeficientes calculamos gamma_m
                gam=cm*(n0/(2*m+1))*cm*(1.0/(1-n0*cm/(2*m+1)))
                #ahora vamos  a calcular gamma(x) con la formula que tiene la sumatoria
                pgm=pol*gam
                #sumatoria
                sum=sum+pgm
            gaf=sum
            sum=0.0
    #comparamos la gamma anterior con la nueva
            resta=np.allclose(gaf, anterior, rtol=1e-05, atol=1e-07, equal_nan=True)
            if(resta==True):
                print(s,n0,sep=',     ')
                break

        if(n0 % 10 == 0):
            name=round(s*100)
            t0=str(name)+str(n0)
            g=gaf+cf+1
            pu=np.log(s*np.exp(gaf)+(1-s)*(1+gaf))-gaf
            np.savetxt(t0+'gmcgb.txt', np.transpose([gaf,cf,g,pu]),delimiter=',   ',header=t0)
            #plt.plot(grad,g)
        n0=n0+5
    s=s+aus
'''
####################################################################################################
n0=20
dn=0.01
n1=n0+dn+dn
df = pd.read_csv("parametros.csv")
df['n'] = n0
df['dn'] = dn
df.to_csv('parametros.csv',index=False)



gaf_dict={}
sum=0.0
s=si
au2=20
nf2=nf+au2
fin2=int((nf2-n0)/au2)

for sc in range(0,finsc):
    n0=20
    n1=n0+dn+dn
    for i in range(0,fin2):
        name=round(s*100)
        sstr=str(name)+str(n0)+'gmcgb.txt'
        if(i<=fin2-1):
            gaf_dict[i]=np.loadtxt(sstr,skiprows=1,usecols=0,delimiter=',  ')
        if(i==fin2):
            gaf_dict[fin2]=gaf_dict[fin2-1]
        gaf=gaf_dict[i]
        for w in range(0,750):
                #calculamos c(x) con la primera aproximación para gamma(x)
            cfh=np.exp(gaf-bt*u)-1-gaf
            cfp=np.exp(-bt*u)*(gaf+1)-1-gaf
            cf=s*cfh+(1-s)*cfp

            anterior=gaf
                #ciclo para calcular las transformadas de Legendre
            for m in range(0,150):
                    #vamos a calcular los coeficientes cm
                pol=(eval_legendre(m, x))
                y = pol*cf
                cm=(2*m+1)*(integrate.simpson(y, x))/2
                #y con esos coeficientes calculamos gamma_m
                gam=cm*(n1/(2*m+1))*cm*(1.0/(1-n1*cm/(2*m+1)))
                    #ahora vamos  a calcular gamma(x) con la formula que tiene la sumatoria
                pgm=pol*gam
                    #sumatoria
                sum=sum+pgm
            gaf=sum
            g1=gaf+cf+1
            sstr2=str(name)+str(n0)+'gant.txt'
            np.savetxt(sstr2, np.transpose([grad,g1]),delimiter=',   ',header=str(n1))
            sum=0.0
            #comparamos la gamma anterior con la nueva
            resta=np.allclose(gaf, anterior, rtol=1e-05, atol=1e-07, equal_nan=True)
            if(resta==True):
                print(s,n1,sep=',     ')
                break
        n1=n1+au2
        n0=n0+au2
    gaf_dict.clear()
    s=s+aus


gaf_dict.clear()

'''
################################################################################
#para la diferencia adelantada
sum=0.0
n0=20
n2=n0+dn
bt=1
s=si
gaf_dict={}


#queremos aumentar n  en pasos de 10


for sc in range(0,finsc):
    n0=20
    n2=n0+dn
    for j in range(0,fin2):
        name=round(s*100)
        sstr=str(name)+str(n0)+'gmcgb.txt'
        if(j<=fin2-1):
            gaf_dict[j]=np.loadtxt(sstr,skiprows=1,usecols=0,delimiter=',  ')
        if(j==fin2):
            gaf_dict[fin2]=gaf_dict[fin2-1]
        gaf=gaf_dict[j]
        for w in range(0,750):
                #calculamos c(x) con la primera aproximación para gamma(x)
            cfh=np.exp(gaf-bt*u)-1-gaf
            cfp=np.exp(-bt*u)*(gaf+1)-1-gaf
            cf=s*cfh+(1-s)*cfp
            #np.savetxt(t0+'cf.txt',np.transpose([cf]))
            anterior=gaf
                #ciclo para calcular las transformadas de Legendre
            for m in range(0,150):
                    #vamos a calcular los coeficientes cm
                pol=(eval_legendre(m, x))
                y = pol*cf
                cm=(2*m+1)*(integrate.simpson(y, x))/2
                #y con esos coeficientes calculamos gamma_m
                gam=cm*(n2/(2*m+1))*cm*(1.0/(1-n2*cm/(2*m+1)))
                    #ahora vamos  a calcular gamma(x) con la formula que tiene la sumatoria
                pgm=pol*gam
                    #sumatoria
                sum=sum+pgm
            gaf=sum
            g2=gaf+cf+1
            sstr2=str(name)+str(n0)+'gsig.txt'
            np.savetxt(sstr2, np.transpose([grad,g2]),delimiter=',   ',header=str(n2))
            sum=0.0
            #comparamos la gamma anterior con la nueva
            resta=np.allclose(gaf, anterior, rtol=1e-05, atol=1e-07, equal_nan=True)
            if(resta==True):
                print(s,n2,sep=',     ')
                break
        n2=n2+au2
        n0=n0+au2
    gaf_dict.clear()
    s=s+aus

#############################################################################
plt.plot(grad,g,label='gant')
plt.plot(grad,g1,label='g')
plt.plot(grad,g2,label='gsig')
plt.legend()

plt.show()
'''
