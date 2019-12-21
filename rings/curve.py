from __future__ import unicode_literals
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import special as sp
from astropy.io import ascii
from astropy.table import Table

def rotcurve():
    """ Вычисление вкладов в кривую вращения отдельних компонетнов галактики и 
    суммарной кривой вращения, которая необходима для работы модели"""
    
    G = 6.67408*1.98892/(3.0856776**2)*(3.1556925**2)/3.0856776*0.1
    fkv = 3.1556925/3.0856776/1000
    theta = 0   
    nn = 400 
    dr = 12.0/nn
    nfi = 1000
    dfi = 2*math.pi/nfi
    
    """ Компоненты разложения потенциала бара """
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
   
    """ Собственно вклад бара. """
    def forcebar(x, y, theta):
        Mb = 0.13 
        ab = 4.2
        cb = 1.35
        eps = math.sqrt(ab**2-cb**2)
        x1 = x*math.cos(theta)+y*math.sin(theta)
        y1 = -x*math.sin(theta)+y*math.cos(theta)
        mu = (math.sqrt(((x**2)/(ab**2))+((y**2)/(cb**2))))
        psi1 = math.acos(math.sqrt(((x**2)-(y**2)-(eps**2)+math.sqrt(4*(x**2)*(y**2)+(((y**2)-(x**2)+(eps**2))**2)))/(2*(x**2))))
        if mu>1:
            if np.abs(x1)>0.01:
                psi = psi1
            else:
                psi = math.acos(y1/(math.sqrt((y1**2)+eps**2)))
        else:
            psi = math.acos(cb/ab)
        abarxref = ((105*G*Mb)/(32*eps))*(-((2*x1*w11(psi))/(eps**2))+((4*(x1**3)*w12(psi)+4*x1*(y1**2)*w21(psi))/(eps**4))-((6*w13(psi)*(x1**5)+12*w22(psi)*(x1**3)*(y1**2)+6*x1*w31(psi)*(y1**4))/(3*(eps**6))))
        abaryref = ((105*G*Mb)/(32*eps))*(-((2*y1*w20(psi))/(eps**2))+((4*y1*w21(psi)*(x1**2)+4*w30(psi)*(y1**3))/(eps**4))-((6*y1*w22(psi)*(x1**4)+12*w31(psi)*(x1**2)*(y1**3)+6*w40(psi)*(y1**5))/(3*(eps**6))))
        abarx = abarxref*math.cos(theta)-abaryref*math.sin(theta)
        abary = abarxref*math.sin(theta)+abaryref*math.cos(theta)
        return(abarx, abary)
        
    """ Вклады остальных компонентов. """     
    def forcebuldge(x, y):
        M = 0.05
        b = 0.30
        r = math.sqrt(x**2+y**2)
        abuldge = -M*G*r/((r**2+b**2)*(math.sqrt(r**2+b**2)))
        axbuldge = abuldge*x/r
        aybuldge = abuldge*y/r
        return(axbuldge, aybuldge)
        
    def forcehalo(x, y):
        Rh = 8
        Vmax = 206*fkv 
        r = math.sqrt(x**2+y**2);
        ahalo = -(Vmax**2)*r/(r**2+Rh**2)
        ax = ahalo*x/r
        ay = ahalo*y/r
        return(ax, ay)
    
    def forcedisk(x,y):
        Md = 0.350 
        rdisk = 2.5
        r = math.sqrt(x**2+y**2)
        ya = r/(2*rdisk)
        sigm = Md/(2*math.pi*rdisk*rdisk*(1-math.exp(-12.0/rdisk)*(1+(12.0/rdisk))))
        adisk = -(4*math.pi*G*sigm*rdisk*ya*ya*(sp.iv(0, ya)*sp.kv(0, ya)-sp.iv(1, ya)*sp.kv(1, ya)))/r
        ax = adisk*x/r
        ay = adisk*y/r 
        return(ax, ay)   
        
    """ Считаем саму кривую. """    
    rotcurve = np.zeros((nn, 6))
    for n in range (1, nn+1, 1):
        if (n//50)*50 == n:
            print('Итерация ' +str(n) + ' из 400.')
        rotcurve[n-1, 0] = n*dr
        asumbar = 0
        asum = 0
        asumdisk = 0
        asumhalo = 0
        asumbuldge = 0
        for j in range (0, nfi, 1):
            fi = (j-0.5)*dfi
            x1 = n*dr*math.cos(fi)
            y1 = n*dr*math.sin(fi)
            axb, ayb = forcebar(x1, y1, theta)
            asumbar = asumbar+axb*math.cos(fi)+ayb*math.sin(fi)
        axd, ayd = forcedisk(x1, y1)
        asumdisk = math.sqrt(axd**2+ayd**2)
        axbg, aybg = forcebuldge(x1, y1)
        asumbuldge = math.sqrt(axbg**2+aybg**2)
        axh, ayh = forcehalo(x1,y1)
        asumhalo = math.sqrt(axh**2+ayh**2)
        asum = abs(asumbar/nfi)+asumbuldge+asumhalo+asumdisk
        rotcurve[n-1, 1] = math.sqrt(rotcurve[n-1, 0]*abs(asum))/fkv
        rotcurve[n-1, 2] = math.sqrt(rotcurve[n-1, 0]*abs(asumbar)/nfi)/fkv
        rotcurve[n-1, 3] = math.sqrt(rotcurve[n-1, 0]*abs(asumdisk))/fkv
        rotcurve[n-1, 4] = math.sqrt(rotcurve[n-1, 0]*abs(asumbuldge))/fkv
        rotcurve[n-1, 5] = math.sqrt(rotcurve[n-1, 0]*abs(asumhalo))/fkv  
        
    """ Сохраняем картиночку. """    
    plt.figure(0)
    plt.plot(rotcurve[:,0], rotcurve[:,1], 'k-', linewidth = 2, label = 'Total')
    plt.plot(rotcurve[:,0], rotcurve[:,2], 'b-', label = 'Bar')
    plt.plot(rotcurve[:,0], rotcurve[:,3], 'g-', label = 'Disk')
    plt.plot(rotcurve[:,0], rotcurve[:,4], 'm-', label = 'Buldge')
    plt.plot(rotcurve[:,0], rotcurve[:,5], 'r-', label = 'Halo')
    plt.xlabel('Distance, kpc')
    plt.ylabel('Velocity, km/s')
    plt.title('Galaxy rotation curve')
    plt.legend(loc = 8, fontsize = 'small', framealpha=0.3)
    plt.grid()
    #plt.savefig('Rotation_curve.eps', dpi = 1000)
    plt.savefig('results/Rotation_curve.png', dpi = 1000)
    
    """ Сохраняем данные. """
    cols = ['Rc, kpc','V_total, km/s','V_bar, km/s','V_disk, km/s','V_halo, km/s','V_buldge, km/s']
    data = Table(rotcurve, names=cols)
    ascii.write(data, 'data/Rotation_curve.dat', overwrite = True)
    ascii.write(data, 'results/Rotation_curve.dat', overwrite = True)
    
