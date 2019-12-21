from __future__ import unicode_literals
import numpy as np
import math
from astropy.table import Table
import random as rd
import rings.wwd

""" Используя записанную кривую вращения, считаем движение газовых облаков в потенциале Галактики. """
G = 6.67408*1.98892/(3.0856776**2)*(3.1556925**2)/3.0856776*0.1
fkv = 3.1556925/3.0856776/1000
rmax = 12.0
rmin = 0.02
ncol = 0
dt = 0.005
dtcol = 0.05
acol = 0.001
ncf = round(dtcol/dt)
obtime = round(4.0/dt); 
Omega = 50*fkv
trot = 2*math.pi/Omega
tgrow = 4.0*trot
Rc, Vtot, Vdisk, Vbar = rings.wwd.read_data()

""" Ускорение бара - с учётом его вращения -> постепенное включение в течение 4 оборотов. """
def forcebar(x, y, theta, fstr, lent):
    Mb = 0.13
    ab = 4.2
    cb = 1.35
    eps = math.sqrt(ab**2-cb**2)
    x1 = x*math.cos(theta)+y*math.sin(theta)
    y1 = -x*math.sin(theta)+y*math.cos(theta)
    r = math.sqrt(x**2 + y**2)
    mu = (math.sqrt(((x**2)/(ab**2))+((y**2)/(cb**2))))
    psi1 = math.acos(math.sqrt(((x**2)-(y**2)-(eps**2)+math.sqrt(4*(x**2)*(y**2)+(((y**2)-(x**2)+(eps**2))**2)))/(2*(x**2))))
    if mu>1:
        if np.abs(x1)>0.01:
            psi = psi1
        else:
            psi = math.acos(y1/(math.sqrt((y1**2)+eps**2)))
    else:
        psi = math.acos(cb/ab)
    f0 = (math.sin(psi)/(((math.cos(psi))**2)))/2
    f1 = math.log(np.abs((math.sin(psi)+1)/math.cos(psi)))
    w10 = 2*f1
    w11 = w10 - 2*math.sin(psi)
    w12 = w11 - (2/3)*((math.sin(psi))**3)
    w13 = w12 - (2/5)*((math.sin(psi))**5)
    w20 = 2*f0 - f1
    w21 = 2*f0 - 3*f1 + 2*(math.sin(psi))
    w22 = 2*f0 - 5*f1 + 4*(math.sin(psi)) + (2/3)*((math.sin(psi))**3)
    w30 = (1/8)*(6*f1 + (10*((math.sin(psi))**3) - 6*(math.sin(psi)))/((math.cos(psi))**4))
    w31 = (1/8)*(30*f1 + ((50*((math.sin(psi))**3) - 30*(math.sin(psi))-16*((math.sin(psi))**5))/((math.cos(psi))**4)))
    w40 = (1/48)*(((66*((math.sin(psi))**5) + 30*(math.sin(psi)) - 80*((math.sin(psi))**3))/((math.cos(psi))**6))-30*f1)
    abarxref = ((105*G*Mb)/(32*eps))*(-((2*x1*w11)/(eps**2))+((4*(x1**3)*w12+4*x1*(y1**2)*w21)/(eps**4))-((6*w13*(x1**5)+12*w22*(x1**3)*(y1**2)+6*x1*w31*(y1**4))/(3*(eps**6))))
    abaryref = ((105*G*Mb)/(32*eps))*(-((2*y1*w20)/(eps**2))+((4*y1*w21*(x1**2)+4*w30*(y1**3))/(eps**4))-((6*y1*w22*(x1**4)+12*w31*(x1**2)*(y1**3)+6*w40*(y1**5))/(3*(eps**6))))
    abarx = abarxref*math.cos(theta)-abaryref*math.sin(theta)
    abary = abarxref*math.sin(theta)+abaryref*math.cos(theta)
    " Постепенное включение. "
    if (fstr < 1):
        irot = lent-1
        for krot in range (0, lent, 1): 
            if ((r < Rc[krot]) and (r >= Rc[krot-1])):
                irot = krot
        vcb = Vbar[irot-1] + (Vbar[irot] - Vbar[irot-1])/(Rc[irot] - Rc[irot-1])*(r - Rc[irot-1])
        abarx = fstr*abarx - (1 - fstr)*vcb**2/r*x1/r*fkv**2
        abary = fstr*abary - (1 - fstr)*vcb**2/r*y1/r*fkv**2
    return(abarx, abary)
    
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
""" Подключаем из кривой вращения, чтобы не считать функции Бесселя на каждой итерации. """
def forcedisk(x, y, i):
    r = math.sqrt(x**2 + y**2)
    adisk = -(Vdisk[i]*fkv)**2/r
    ay = adisk*y/r 
    ax = adisk*x/r
    return(ax, ay)
    
def forcetot(x, y, theta, fstr, lent, i):
    axbar, aybar = forcebar(x, y, theta, fstr, lent)
    axbuldge, aybuldge = forcebuldge(x, y)
    axdisk, aydisk = forcedisk(x, y, i)
    axhalo, ayhalo = forcehalo(x, y)
    ax = axbar + axbuldge + axhalo + axdisk
    ay = aybar + aybuldge + ayhalo + aydisk
    return(ax, ay)
    
""" Генерируем начальные параметры облаков в заданных пределах. """
def set_init(num, lent):
    init_i = np.zeros((num, 6))
    cols = ['x','y','vx','vy','stat','ind']
    init_i = Table(init_i, names=cols)
    n = 0
    while n < num:
        init_i['stat'][n] = 1
        init_i['ind'][n] = n
        x = rmax*rd.uniform(-1, 1) #Пределы по координатам задаём руками.
        y = rmax*rd.uniform(-1, 1)
        r = math.sqrt(x**2 + y**2)
        cos = x/r
        sin = y/r
        if (r > rmin) and (r < rmax):
            init_i['x'][n] = x
            init_i['y'][n] = y
            for k in range (0, lent, 1):
                if (r < Rc[k]) and (r >= Rc[k-1]):
                    l = k
            vt = Vtot[l-1] + (Vtot[l] - Vtot[l-1])/(Rc[l] - Rc[l-1])*(r - Rc[l-1])
            vt = vt + 5.0*rd.uniform(-1, 1) #Пределы по скорости берём из кривой вращения.
            vr = 5.0*rd.uniform(-1, 1)
            init_i['vx'][n] = (vr*cos - vt*sin)*fkv
            init_i['vy'][n] = (vr*sin + vt*cos)*fkv
            n += 1
    np.sort(init_i, order = 'x') #Сортируем, чтобы при столкновениях перебирать только соседние частицы.
    return(init_i)
    
""" Первый шаг моделирования. """   
def first_step(init_i, lent):
    arr_i = np.zeros((len(init_i), 8))
    cols = ['x','y','vx','vy','stat','ind','xb','yb']
    arr_i = Table(arr_i, names=cols)
    for n in range (0, len(init_i), 1):
        arr_i['stat'][n] = init_i['stat'][n]
        arr_i['ind'][n] = init_i['ind'][n]
        x0 = init_i['x'][n]
        y0 = init_i['y'][n]
        vx0 = init_i['vx'][n]
        vy0 = init_i['vy'][n]
        ax, ay = forcetot(x0, y0, 0, 1, lent, 0)
        x2 = x0 + vx0*dt + 0.5*ax*dt**2
        y2 = y0 + vy0*dt + 0.5*ay*dt**2
        vx2 = vx0 + ax *dt/2
        vy2 = vy0 + ay*dt/2
        arr_i['x'][n] = x2
        arr_i['y'][n] = y2
        arr_i['vx'][n] = vx2
        arr_i['vy'][n] = vy2
        r = math.sqrt(x2**2 + y2**2)
        if((r < rmin) or (r > rmax)):
            arr_i['stat'][n] = 0
    return(arr_i)

""" Основной цикл. Выводим для построения картинки только координаты. """    
def main_cycle(data_arr, ii, lent, ncol):
    time = dt*ii
    if (time < tgrow):
        fstr = time/tgrow #степень включения бара.
    theta = Omega*time
    iscol = 0
    if ((ii // ncf)*ncf == ii): 
        iscol = 1
    if (iscol == 1):
        np.sort(data_arr, order = 'x')
    for n in range (0, len(data_arr)-1, 1): #{4}
        isend = 1
        if ((data_arr['stat'][n] == 1) or (data_arr['stat'][n] == ii)):
            isend = 0
        k = n
        while (isend == 0):
            k += 1
            if (data_arr['x'][k] - data_arr['x'][n] > acol):
                isend = 1
            if (isend == 0 and abs(data_arr['y'][k] - data_arr['y'][n]) < acol and data_arr['stat'][k] == 1 or data_arr['stat'][k] == ii):
                ncol += 1
                vx1 = (data_arr['vx'][n] + data_arr['vx'][k])/2
                vy1 = (data_arr['vy'][n] + data_arr['vy'][k])/2
                data_arr['vx'][n] = vx1
                data_arr['vy'][n] = vy1 
                data_arr['vx'][k] = vx1
                data_arr['vy'][k] = vy1
                if (ncol // 10)*10 == ncol:
                    if (data_arr['stat'][n] == 1 and data_arr['stat'][k] == ii):
                        data_arr['stat'][n] = ii
                    if ((data_arr['stat'][n] == ii) and (data_arr['stat'][k] == 1)):
                       data_arr['stat'][k] == ii
                    if ((data_arr['stat'][n] == 1) and (data_arr['stat'][k] == 1)):
                        if (ncol // 20)*20 == ncol:
                            data_arr['stat'][n] = ii 
                        else:
                            data_arr['stat'][k] = ii
            if (k + 1 >= len(data_arr)):
                isend = 1 
    for n in range (0, len (data_arr), 1): 
        if (data_arr['stat'][n] > 1) and (ii == data_arr['stat'][n] + obtime):
            data_arr['stat'][n] = 1
    x2 = data_arr['x'] + data_arr['vx'] * dt
    for n in range (0, len(data_arr), 1):
        if (data_arr['stat'][n] > 0):
            x1 = data_arr['x'][n]
            y1 = data_arr['y'][n]
            vx1 = data_arr['vx'][n]
            vy1 = data_arr['vy'][n]
            x2 = x1 + vx1*dt
            y2 = y1 + vy1*dt
            ax, ay = forcetot(x2, y2, theta, fstr, lent, ii)
            vx2 = vx1 + ax*dt
            vy2 = vy1 + ay*dt
            data_arr['x'][n] = x2
            data_arr['y'][n] = y2
            data_arr['vx'][n] = vx2
            data_arr['vy'][n] = vy2
            xr = x2*math.cos(theta) + y2*math.sin(theta)
            yr = -x2*math.sin(theta) + y2*math.cos(theta)
            data_arr['xb'][n] = xr
            data_arr['yb'][n] = yr
            rr = math.sqrt(xr**2 + yr**2)
            if ((rr < rmin) or (rr > rmax)):
                data_arr['stat'][n] = 0	    
    return(data_arr)
