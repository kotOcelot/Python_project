# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 18:53:44 2019

@author: Ocelot
"""

import numpy as np
import math
import matplotlib.pyplot as plt

G=0.45 #kpc^3*Myr^-2*(10^10*Ms)^-1
dt=0.01 #шаг сетки десять тысяч лет
Mb=1 #*10^10Ms - бар
ab=3.98 #kpc - бар
cb=1.25 #kpc - бар
eps=math.sqrt(ab**2-cb**2) #- бар
time=0
Omega=0.05 
Rh=3.6 #kpc
Vmax=0.22 #kpc/Myr
M=1 #*10^10Ms - балдж
b=0.32 #kpc - балдж
Time=1000
num=int(Time/dt)


def f0(psi):
	return ((math.sin(psi)/(((math.cos(psi))**2)))/2)

def f1(psi):
	return (math.log(np.abs((math.sin(psi)+1)/math.cos(psi))))

def w10(psi):
	return(2*f1(psi))

def w11 (psi):
	return (w10(psi) - 2*math.sin(psi))

def w12(psi):
	return (w11(psi) - (2/3)*((math.sin(psi))**3))

def w13(psi):
	return (w12(psi) - (2/5)*((math.sin(psi))**5))

def w20(psi):
	return (2*f0(psi) - f1(psi))

def w21(psi):
	return (2*f0(psi) - 3*f1(psi) + 2*(math.sin(psi)))

def w22(psi):
	return (2*f0(psi) - 5*f1(psi) + 4*(math.sin(psi)) + (2/3)*((math.sin(psi))**3))

def w30(psi):
	return((1/8)*(6*f1(psi)+(10*((math.sin(psi))**3)-6*(math.sin(psi)))/((math.cos(psi))**4)))

def w31(psi):
	return((1/8)*(30*f1(psi)+((50*((math.sin(psi))**3)-30*(math.sin(psi))-16*((math.sin(psi))**5))/((math.cos(psi))**4))))

def w40(psi):
	return((1/48)*(((66*((math.sin(psi))**5)+30*(math.sin(psi))-80*((math.sin(psi))**3))/((math.cos(psi))**6))-30*f1(psi)))


def mu(x,y):
	return(math.sqrt(((x**2)/(ab**2))+((y**2)/(cb**2))))

def psi1(x,y):
	return(math.acos(math.sqrt(((x**2)-(y**2)-(eps**2)+math.sqrt(4*(x**2)*(y**2)+(((y**2)-(x**2)+(eps**2))**2)))/(2*(x**2)))))

test=np.zeros((num))
x=np.zeros((num))
y=np.zeros((num))
vx=np.zeros((num))
vy=np.zeros((num))
test=np.zeros((num))
xout=np.zeros((num))
yout=np.zeros((num))
x[0]=8
y[0]=0
vx[0]=0.0
vy[0]=0.237

def force(x, y, teta):
    x1=x*math.cos(teta)+y*math.sin(teta)
    y1=-x*math.sin(teta)+y*math.cos(teta)
    if mu(x1,y1)>1:
        if np.abs(x1)>0.01:
            psi=psi1(x1,y1)
        else:
            psi=math.acos(y1/(math.sqrt((y1**2)+eps**2)))
    else:
        psi=math.acos(cb/ab)
    W11=w11(psi)
    W12=w12(psi)
    W13=w13(psi)
    W20=w20(psi)
    W21=w21(psi)
    W22=w22(psi)
    W30=w30(psi)
    W31=w31(psi)
    W40=w40(psi)
    abarxref=((105*G*Mb)/(32*eps))*(-((2*x1*W11)/(eps**2))+((4*(x1**3)*W12+4*x1*(y1**2)*W21)/(eps**4))-((6*W13*(x1**5)+12*W22*(x1**3)*(y1**2)+6*x1*W31*(y1**4))/(3*(eps**6))))
    abaryref=((105*G*Mb)/(32*eps))*(-((2*y1*W20)/(eps**2))+((4*y1*W21*(x1**2)+4*W30*(y1**3))/(eps**4))-((6*y1*W22*(x1**4)+12*W31*(x1**2)*(y1**3)+6*W40*(y1**5))/(3*(eps**6))))
    abarx=abarxref*math.cos(teta)-abaryref*math.sin(teta)
    abary=abarxref*math.sin(teta)+abaryref*math.cos(teta)
    ahalox=-x*((Vmax**2)/((Rh**2)+(x**2)+(y**2)))
    abuldgex=-((G*M*x)/((((x)**2)+((y)**2)+(b**2))**(3/2)))
    abuldgey=-((G*M*y)/((((x)**2)+((y)**2)+(b**2))**(3/2)))
    ahaloy=-(y)*((Vmax**2)/((Rh**2)+((x**2)+(y**2))))
    ax=abuldgex+ahalox+abuldgex
    ay=abuldgey+ahaloy+abuldgey
    return(ax, ay)

for i in range (1,2,1):
    time=time+dt
    teta=time*Omega
    ax=force(x[i-1],y[i-1],teta)[0]
    ay=force(x[i-1],y[i-1],teta)[1]
    x[i]=x[i-1]+dt*vx[i-1]+(1/2)*(dt**2)*ax
    y[i]=y[i-1]+dt*vy[i-1]+(1/2)*(dt**2)*ay
    vx[i]=vx[i-1]+(1/2)*dt*ax
    vy[i]=vy[i-1]+(1/2)*dt*ay
    xout[i]=x[i]*math.cos(teta)+y[i]*math.sin(teta)
    yout[i]=x[i]*math.sin(teta)+y[i]*math.cos(teta)
    
    
for i in range (2,num,1):
    time=time+dt
    teta=time*Omega
    ax=force(x[i-1],y[i-1],teta)[0]
    ay=force(x[i-1],y[i-1],teta)[1]
    x[i]=x[i-1]+dt*vx[i-1]
    y[i]=y[i-1]+dt*vy[i-1]
    vx[i]=vx[i-1]+dt*ax
    vy[i]=vy[i-1]+dt*ay
    xout[i]=x[i]*math.cos(teta)+y[i]*math.sin(teta)
    yout[i]=x[i]*math.sin(teta)+y[i]*math.cos(teta)


plt.plot(x[1:], y[1:])
#plt.show()
#np.savetxt('x_python.txt', x)
#np.savetxt('y_python.txt', y)
#np.savetxt('vx_python.txt', vx)
#np.savetxt('vy_python.txt', vy)
