import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8})
import numpy as np
import math as mt
from scipy import integrate
import pandas as pd


SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


datos=pd.read_csv('resultados.csv')
df1=pd.DataFrame(datos)
#n=df1[['n']]
#
#co=df1[['camino compresibilidad']]
#vi=df1[['camino virial']]

ind20 = []
ind40 = []
ind60 = []
ind80 = []
ind100 = []
ind120 = []

for i in range(len(df1.n)):
    if 20 == df1.n[i]:
        ind20.append(i)
    if 40 == df1.n[i]:
        ind40.append(i)
    if 60 == df1.n[i]:
        ind60.append(i)
    if 80 == df1.n[i]:
        ind80.append(i)
    if 100 == df1.n[i]:
        ind100.append(i)
    if 120 == df1.n[i]:
        ind120.append(i)

df20 =pd.DataFrame()
df40 = pd.DataFrame()
df60 = pd.DataFrame()
df80 = pd.DataFrame()
df100 = pd.DataFrame()
df120 = pd.DataFrame()


for indexes in ind20:
    df20 = df20.append(df1.iloc[indexes])
for indexes in ind40:
    df40 = df40.append(df1.iloc[indexes])
for indexes in ind60:
    df60 = df60.append(df1.iloc[indexes])
for indexes in ind80:
    df80 = df80.append(df1.iloc[indexes])
for indexes in ind100:
    df100 = df100.append(df1.iloc[indexes])
for indexes in ind120:
    df120 = df120.append(df1.iloc[indexes])
#df2 = df2.where(df1.n == 40.1)

df20.dropna()
df40.dropna()
df60.dropna()
df80.dropna()
df100.dropna()
df120.dropna()


svec=np.linspace(1,2.65,num=1000)

s=df60[['s']].values
s=s.ravel()
#s=np.linspace(0,1,num=)
print(s)

c20=df20[['camino compresibilidad']].values
c20=c20.ravel()
pc20=np.polyfit(s,c20,4)
poly_coeff = np.polynomial.polynomial.polyfit(s,c20, 4)
valc20=np.polyval(pc20,svec)
v20=df20[['camino virial']].values
v20=v20.ravel()
pv20=np.polyfit(s,v20,4)
x020=2.322
y020=np.polyval(pv20,x020)
valv20=np.polyval(pv20,svec)
poly_coeff2 = np.polynomial.polynomial.polyfit(s,v20, 4)
root20 = np.polynomial.polynomial.polyroots(poly_coeff - poly_coeff2)
print('la interseccion en 20', root20)


c40=df40[['camino compresibilidad']].values
c40=c40.ravel()
pc40=np.polyfit(s,c40,3)
poly_coeff = np.polynomial.polynomial.polyfit(s,c40, 3)
valc40=np.polyval(pc40,svec)
v40=df40[['camino virial']].values
v40=v40.ravel()
pv40=np.polyfit(s,v40,3)
valv40=np.polyval(pv40,svec)
x040=2.5446
y040=np.polyval(pv40,x040)
poly_coeff2 = np.polynomial.polynomial.polyfit(s,v40, 3)
root40 = np.polynomial.polynomial.polyroots(poly_coeff - poly_coeff2)
print('la interseccion en 40',root40)


c60=df60[['camino compresibilidad']].values
c60=c60.ravel()
pc60=np.polyfit(s,c60,3)
poly_coeff = np.polynomial.polynomial.polyfit(s,c60, 3)
valc60=np.polyval(pc60,svec)
v60=df60[['camino virial']].values
v60=v60.ravel()
pv60=np.polyfit(s,v60,3)
x060=2.85
y060=np.polyval(pv60,x060)
poly_coeff2 = np.polynomial.polynomial.polyfit(s,v60, 3)
root60 = np.polynomial.polynomial.polyroots(poly_coeff - poly_coeff2)
valv60=np.polyval(pv60,svec)
print('la interseccion en 60',root60)

nans=df80['camino compresibilidad'].isnull().sum(axis=0)
print('el numero de nans es:',nans)
df80.drop(df80.tail(nans).index,inplace=True)
s=df80[['s']].values
s=s.ravel()
c80=df80[['camino compresibilidad']].values
c80=c80.ravel()
pc80=np.polyfit(s,c80,4)
poly_coeff = np.polynomial.polynomial.polyfit(s,c80, 4)
valc80=np.polyval(pc80,svec)
v80=df80[['camino virial']].values
v80=v80.ravel()
pv80=np.polyfit(s,v80,4)
x080=2.6021
y080=np.polyval(pv80,x080)
poly_coeff2 = np.polynomial.polynomial.polyfit(s,v80, 4)
root80 = np.polynomial.polynomial.polyroots(poly_coeff - poly_coeff2)
valv80=np.polyval(pv80,svec)
print('la interseccion en 80',root80)



nans=df100['camino compresibilidad'].isna().sum()
df100.drop(df100.tail(nans).index,inplace=True)
s=df100[['s']].values
s=s.ravel()
c100=df100[['camino compresibilidad']].values
c100=c100.ravel()
print(c100)
pc100=np.polyfit(s,c100,3)
poly_coeff = np.polynomial.polynomial.polyfit(s,c100, 3)
valc100=np.polyval(pc100,svec)
v100=df100[['camino virial']].values
v100=v100.ravel()
pv100=np.polyfit(s,v100,3)
x0100=2.7109
y0100=np.polyval(pv100,x080)
poly_coeff2 = np.polynomial.polynomial.polyfit(s,v100, 3)
root100 = np.polynomial.polynomial.polyroots(poly_coeff - poly_coeff2)
valv100=np.polyval(pv100,svec)
print('la interseccion en 100',root100)

nans=df120['camino compresibilidad'].isna().sum()
df120.drop(df120.tail(nans).index,inplace=True)
s=df120[['s']].values
s=s.ravel()
c120=df120[['camino compresibilidad']].values
c120=c120.ravel()
pc120=np.polyfit(s,c120,3)
poly_coeff = np.polynomial.polynomial.polyfit(s,c120, 3)
valc120=np.polyval(pc120,svec)
v120=df120[['camino virial']].values
v120=v120.ravel()
pv120=np.polyfit(s,v120,3)
poly_coeff2 = np.polynomial.polynomial.polyfit(s,v120, 3)
x0120=2.67
y0120=np.polyval(pv120,x0120)
root120 = np.polynomial.polynomial.polyroots(poly_coeff - poly_coeff2)
valv120=np.polyval(pv120,svec)
print('la interseccion en 120',root120)


s=df60[['s']].values
fig, ax = plt.subplots(2,3,sharex=True)

ax[0][0].plot(svec,valv20,'--c')
ax[0][0].plot(svec,valc20,'--y')
ax[0][0].plot(s,v20,'.b', markersize=5,label='20v')
ax[0][0].plot(s,c20,'.g', markersize=5,label='20v')
ax[0][0].plot(x020,y020,'.r',markersize=8)
#ax[0][0].legend(loc='center right')
ax[0][0].set_title('N=20')
#ax[0][0].set_xlabel('s')
ax[0][0].set_ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
#ax[0][0].savefig('comparacion/20.png')
#ax[0][0].clf()

ax[0][1].plot(svec,valc40,'--y')
ax[0][1].plot(svec,valv40,'--c')
ax[0][1].set_title('N=40')
ax[0][1].plot(s,c40,'.g', markersize=5,label='40c')
ax[0][1].plot(s,v40,'.b', markersize=5,label='40v')
ax[0][1].plot(x040,y040,'.r',markersize=8)
#ax[0][1].legend(loc='center right')
#ax[0][1].set_xlabel('s')
#ax[0][1].set_ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
#ax[0][1].savefig('comparacion/40.png')
#ax[0][1].clf()

ax[0][2].plot(svec,valv60,'--c')
ax[0][2].plot(svec,valc60,'--y')
ax[0][2].set_title('N=60')
ax[0][2].plot(s,c60,'.g', markersize=5,label='60c')
ax[0][2].plot(s,v60,'.b', markersize=5,label='60v')
ax[0][2].plot(x060,y060,'.r',markersize=8)
#ax[0][2].legend(loc='center right')
#ax[0][2].set_xlabel('s')
#ax[0][2].set_ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
#ax[0][2].savefig('comparacion/60.png')
#ax[0][2].clf()

s=df80[['s']].values
ax[1][0].plot(svec,valv80,'--c')
ax[1][0].plot(svec,valc80,'--y')
ax[1][0].set_title('N=80')
ax[1][0].plot(s,c80,'.g', markersize=5,label='80c')
ax[1][0].plot(s,v80,'.b', markersize=5,label='80v')
ax[1][0].plot(x080,y080,'.r',markersize=8)
ax[1][0].set_xlabel('s')
ax[1][0].set_ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
#ax[1][0].legend(loc='center right')
#plt.savefig('comparacion/80.png')
#plt.clf()
s=df100[['s']].values
ax[1][1].plot(svec,valv100,'--c')
ax[1][1].plot(svec,valc100,'--y')
ax[1][1].set_title('N=100')
ax[1][1].plot(s,c100,'.g', markersize=5,label='100c')
ax[1][1].plot(x0100,y0100,'.r',markersize=8)
ax[1][1].plot(s,v100,'.b', markersize=5,label='100v')
ax[1][1].set_xlabel('s')
#ax[1][1].set_ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
#ax[1][1].legend(loc='center right')
#plt.savefig('comparacion/100.png')
#plt.clf()

s=df120[['s']].values
ax[1][2].plot(svec,valv120,'--c')
ax[1][2].plot(svec,valc120,'--y')
ax[1][2].set_title('N=120')
ax[1][2].plot(s,c120,'.g', markersize=5,label='c')
ax[1][2].plot(x0120,y0120,'.r',markersize=8)
ax[1][2].plot(s,v120,'.b', markersize=5,label='v')
ax[1][2].set_xlabel('s')
ax[1][2].legend(bbox_to_anchor=(1, 1),
          bbox_transform=fig.transFigure)
#ax[1][2].set_ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
#ax[1][2].legend(loc='center right')
#plt.savefig('comparacion/120.png')
#plt.clf()
#handles, labels = ax.get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper center')
#fig.legend('')
plt.subplots_adjust(left=0.125,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.3,
                    hspace=0.35)
fig.savefig('comparacion/comparacion.png',dpi=270)
print(df40)

plt.clf()




s=df60[['s']].values

plt.plot(svec,valv20,'--c')
plt.plot(svec,valc20,'--y')
plt.plot(s,v20,'.b', markersize=5,label='20v')
plt.plot(s,c20,'.g',markersize=5,label='20c')
plt.plot(x020,y020,'.r',markersize=8)
plt.title('N=20')
plt.xlabel('s')
plt.savefig('comparacion/20.png')
plt.clf()

plt.plot(svec,valc40,'--y')
plt.plot(svec,valv40,'--c')
plt.plot(s,c40,'.g', markersize=5,label='40c')
plt.plot(s,v40,'.b', markersize=5,label='40v')
plt.plot(x040,y040,'.r',markersize=8)
plt.legend(loc='center right')
plt.xlabel('s')
plt.ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
plt.savefig('comparacion/40.png')
plt.clf()

plt.plot(svec,valv60,'--c')
plt.plot(svec,valc60,'--y')
plt.title('N=60')
plt.plot(s,c60,'.g', markersize=5,label='60c')
plt.plot(s,v60,'.b', markersize=5,label='60v')
plt.plot(x060,y060,'.r',markersize=8)
plt.legend(loc='center right')
plt.xlabel('s')
plt.ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
plt.savefig('comparacion/60.png')
plt.clf()


s=df80[['s']].values
plt.plot(svec,valv80,'--c')
plt.plot(svec,valc80,'--y')
plt.title('N=80')
plt.plot(s,c80,'.g', markersize=5,label='80c')
plt.plot(s,v80,'.b', markersize=5,label='80v')
plt.plot(x080,y080,'.r',markersize=8)
plt.ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
plt.legend(loc='center right')
plt.savefig('comparacion/80.png')
plt.clf()

s=df100[['s']].values

plt.plot(svec,valv100,'--c')
plt.plot(svec,valc100,'--y')
plt.title('N=100')
plt.plot(s,c100,'.g', markersize=5,label='100c')
plt.plot(s,v100,'.b', markersize=5,label='100v')
plt.xlabel('s')
plt.ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
plt.legend(loc='center right')
plt.savefig('comparacion/100.png')
plt.clf()

s=df120[['s']].values
plt.plot(svec,valv120,'--c')
plt.plot(svec,valc120,'--y')
plt.title('N=120')
plt.plot(s,c120,'.g', markersize=5,label='c')
plt.plot(x0120,y0120,'.r',markersize=8)
plt.plot(s,v120,'.b', markersize=5,label='v')
plt.xlabel('s')
#ax[1][2].legend(bbox_to_anchor=(1, 1),
#          bbox_transform=fig.transFigure)
#ax[1][2].set_ylabel(r'$A  \frac{\partial \beta P}{\partial N}$')
#ax[1][2].legend(loc='center right')
plt.savefig('comparacion/120.png')
plt.clf()



'''
g4035=np.loadtxt("40g35.txt",skiprows=1,usecols=1)
g4060=np.loadtxt("40g60.txt",skiprows=1,usecols=1)
#g4078=np.loadtxt("40gant78.txt",skiprows=1,usecols=0)
#g4080=np.loadtxt("40gant80.txt",skiprows=1,usecols=0)
x = np.linspace(-1, 1, num=20001)
ang=np.arccos(x)
grad=ang*180/(mt.pi)



plt.plot(grad,g4035,'-', markersize=8,label='w=0.35')
plt.plot(grad,g4060,'-', markersize=8,label='w=0.60')
#plt.plot(grad,g4078,'-', markersize=8,label='w=0.78')
#plt.plot(grad,g4080,'-', markersize=8,label='w=0.80')
plt.legend(loc='center right')
plt.title('n=40')
plt.show()



grad=np.loadtxt("grad.txt",usecols=0,delimiter=',')
gaf=np.loadtxt('20gaf35.txt',skiprows=1)
m=.35
b=np.log(m*np.exp(gaf)+(1-m)*(1+gaf))-gaf
plt.plot(grad,b)
plt.title('n=20, w=0.35')
plt.xlabel(r'$\theta$')
plt.ylabel('b')
plt.show()
'''
