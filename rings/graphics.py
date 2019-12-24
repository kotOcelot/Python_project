from __future__ import unicode_literals
import matplotlib.pyplot as plt
from astropy.io import ascii
import matplotlib.ticker
from matplotlib.patches import Ellipse

def plot_distr(num, path) :
    """ Строим распределение облаков по диску галактики. """
    leni = 12
    res = 100
    dms = 2.54
    fs = leni / dms
    ab = 4.2*2
    cb = 1.35*2
    numpart = 10000
    data = ascii.read(str(path) + '/time_'+str(num)+'.dat', header_start=0, data_start = 1)
    if num == 0 or path == 'data':
        X = data['x']
        Y = data['y']
    else:
        X = data['xb']
        Y = data['yb']
    stat = data['stat']
    plt.figure(figsize=(fs, fs), dpi = res)
    ax = plt.subplot(1,1,1)
    for i in range (0, numpart, 1):
        if (i//500)*500 == i:
            print('Итерация ' +str(i) + ' из 10000.')
        if stat[i] == 1:
            plt.plot(X[i], Y[i], 'k.', markersize = 1)
        else:
            if stat[i] > 1:
             plt.plot(X[i], Y[i], 'r.', markersize = 1)              
    locator = matplotlib.ticker.MultipleLocator(base=6)
    ax.xaxis.set_major_locator(locator)
    ax.yaxis.set_major_locator(locator)
    bar = Ellipse((0, 0), ab, cb, linewidth = 1, fill = False)
    ax.add_patch(bar)
    plt.xlabel('x, kpc')
    plt.ylabel('y, kpc')
    plt.text(7, 11.5, 'T = ' + str(num) + ' Myr', fontsize = 8)
    #plt.savefig('data/distr_' + str(num) + '.eps')
    plt.savefig('results/distr_' + str(num) + '.png')


def plot_orbit(xr, path):
    """ Строим орбиту отдельного облака. """
    num = 400
    leni = 12
    res = 100
    dms = 2.54
    fs = leni / dms
    ab = 4.2*2
    cb = 1.35*2
    data = ascii.read(str(path) + '/time_0.dat', header_start=0, data_start = 1)
    X = data['x']
    I = data['ind']
    n=0
    for i in range(0, len(X)-1):
        if ((X[i] < xr) and (X[i+1] >= xr)):
            n = I[i]
            xj = X[i]
            yj = data['y'][i]
            yf = yj
            break
    plt.figure(figsize=(fs, fs), dpi = res)
    ax = plt.subplot(1,1,1)
    for k in range(1, num):
        if (k//50)*50 == k:
            print('Итерация ' +str(k) + ' из 400.')
        data = ascii.read(str(path) + '/time_' + str(k*5) +'.dat', header_start=0, data_start = 1)
        for l in range(0, len(X)):
            if (data['ind'][l] == n):
                if (data['stat'][l] != 0):
                    xi = data['x'][l]
                    yi = data['y'][l]
                    if (data['stat'][l] == 1):
                        plt.plot(xi/dms, yi/dms, 'k.', markersize = 3)
                        plt.plot([xj/dms, xi/dms], [yj/dms, yi/dms], 'k-', linewidth = 1)
                    else:
                        if data['stat'][l] > 1:
                            plt.plot(xi/dms, yi/dms, 'r.', markersize = 3)
                            plt.plot([xj/dms, xi/dms], [yj/dms, yi/dms], 'r-', linewidth = 1)
                    xj = data['x'][l]
                    yj = data['y'][l]
                else:
                    break
                break
    bar = Ellipse((0, 0), ab, cb, linewidth = 1, fill = False)
    ax.add_patch(bar)
    plt.text(-0.25, 0.25, 'x = ' + str(round(xr,3)) + ' kpc', fontsize = 7)
    plt.text(-0.25, -0.25, 'y = ' + str(round(yf, 3)) + ' kpc', fontsize = 7)
    plt.xlabel('x, kpc')
    plt.ylabel('y, kpc')
    #plt.savefig('results/orbit_' + str(n) + '.eps')
    plt.savefig('results/orbit_' + str(n) + '.png')

