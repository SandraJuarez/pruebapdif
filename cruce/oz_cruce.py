import matplotlib.pyplot as plt
import numpy as np
import math as mt
from scipy.special import eval_legendre
from scipy import integrate
import pandas as pd



gaf=0.0
bt=1#2.4987e20 #kb*T lo tomamos como unidad t=290
a=100

x = np.linspace(-1, 1, num=20001)
ang=np.arccos(x)
grad=ang*180/(mt.pi) #cambiamos de radianes a grados
r=2*a*np.sin(ang/2)
np.place(r, r==0, [1e-9])
k=1.0/9.6
b=1510.01 #amplitud del potencial
#b=0
n=0 #número de partículas
al=1.0/0.1 #alfa
u=b*np.exp(-k*r)*(1-np.exp(-al*r))/(r)
np.savetxt('potencial.txt', np.transpose([u,r]),header='Potencial')

#guardamos nuestros parámetros
sum=0.0

##############################################################################
s=2.54
nf=40
nf2=nf+1
################################################################################
    #ciclo que va aumentando N en cada iteración
for nc in range(0,nf2):
        #ciclo para iterar cada vez con una nueva gamma
        #la letra f significa función y la m cuando sea un coeficiente
    for w in range(0,750):
            #calculamos c(x) con la primera aproximación para gamma(x)
        cfh=np.exp(gaf-bt*u)-1-gaf
        cfp=np.exp(-bt*u)*(gaf+1)-1-gaf
        cf=s*cfh+(1-s)*cfp
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
            gam=cm*(n/(2*m+1))*cm*(1.0/(1-n*cm/(2*m+1)))
                #ahora vamos  a calcular gamma(x) con la formula que tiene la sumatoria
            pgm=pol*gam
                #sumatoria
            sum=sum+pgm
        gaf=sum
        sum=0.0
        resta=np.allclose(gaf, anterior, rtol=1e-05, atol=1e-07, equal_nan=True)
        if(resta==True):
            print(w,n,sep=',     ')
            break
    n=n+1

#para poner en el encabezado el valor de n


#ahora calculamos g(x)=gamma(x)+c(x)+1
#g=gaf+cf+1
#grad=ang*180/(mt.pi) #cambiamos de radianes a grados

#guardamos los resultados en un archivo de texto
t0=str(nf)
g=gaf+cf+1
pu=np.log(s*np.exp(gaf)+(1-s)*(1+gaf))-gaf
np.savetxt(t0+'gmcgb.txt', np.transpose([gaf,cf,g,pu]),delimiter=',   ',header=t0)
plt.plot(grad,g,'-r')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$g(\theta)$')
plt.title('n='+t0+ ', s='+str(s))
plt.savefig(t0+'g.png')
plt.clf()
plt.plot(grad,cf,'-r')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$c(\theta)$')
plt.title('n='+t0+', s='+str(s))
plt.savefig(t0+'c.png')
plt.clf()
plt.plot(grad,pu,'-r')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$b(\theta)$')
plt.title('n='+t0+', s='+str(s))
plt.savefig(t0+'b.png')
