import numpy as np
import matplotlib.pyplot as plt
from scipy.special import eval_legendre
from scipy import integrate

dim=1001
sep=dim
m=dim

D=np.zeros((dim,dim))
x1=-1
dx=2/sep

for n in range(0,sep):
    w=sep-1
    for l in range(0,m):
        pol=eval_legendre(l,x1)
        fun=pol
        m1=dx*(2*l+1)/6
        #m1=dx/3
        mod=n%2
        if(mod==0)and(n!=0)and(n!=w):
            D[l,n]=m1*2*fun
        elif(mod!=0)and(n!=0)and(n!=w):
            D[l,n]=m1*4*fun
        elif(n==0):
            D[l,n]=m1*fun
        elif(n==w):
            D[l,n]=m1*fun
    x1=x1+dx
bt=1
u=1
x = np.linspace(-1, 1, num=sep)
gaf=np.loadtxt('20gmcgb.txt',usecols=0,skiprows=1,delimiter=', ')
print(np.shape(gaf))
cf=np.exp(gaf-bt*u)-1-gaf
ex=np.exp(gaf-bt*u)
    #print(cf)
    #calculamos las sumatorias:

c=np.zeros((m))
gam=np.zeros((m))
cg=np.zeros((m))

for l in range(0,m):
    sum=0.0
    sum2=0.0
    sum3=0.0
    for n in range(0,sep):
        sum+=D[l,n]
        sum2+=D[l,n]*(ex[n]-1)
        sum3+=D[l,n]*gaf[n]
    c[l]=sum
    cg[l]=sum2
    gam[l]=sum3

print(c[123])

c=np.zeros((m))
gam=np.zeros((m))
for l in range(0,m):
    #vamos a calcular los coeficientes cm
    pol=(eval_legendre(l, x))
    y = pol*cf
    c[l]=(2*l+1)*(integrate.simpson(y, dx=dx))/2
    #y con esos coeficientes calculamos gamma_m
    gam[l]=c[l]*(N/(2*l+1))*c[l]*(1.0/(1-N*c[l]/(2*l+1)))
print('Usando simpson obtenemos:',c[1])
