from __future__ import unicode_literals
from astropy.io import ascii

def read_data():
    """ Чтение данных из кривой вращения. """
    data = ascii.read('data/Rotation_curve.dat', header_start = 0, data_start = 1)
    Rc1 = data['Rc, kpc']
    Vtot1 = data['V_total, km/s']
    Vdisk1 = data['V_disk, km/s']
    Vbar1 = data['V_bar, km/s']
    return(Rc1, Vtot1, Vdisk1, Vbar1)
    
def write_init(init_i, k):
    """ Записываем координаты, индексы и статусы газовых облаков после первой итерации. """
    ascii.write(init_i['x','y','stat', 'ind'], 'results/time_'+str(k)+'.dat', overwrite = True)

def write_res(darr_i, k):
    """ Записываем промежуточные данные. """
    ascii.write(darr_i['xb','yb','stat', 'ind'], 'results/time_'+str(k)+'.dat', overwrite = True)
    