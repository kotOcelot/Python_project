from __future__ import unicode_literals
import math
from multiprocessing import Process
from time import time
from scipy import special as sp
import scipy.integrate as integrate
import unittest
import rings.galaxy
import rings.wwd

""" Тестируем, правильно ли мы считаем интеграл в кривой вращения, функции Бесселя и время. """
class testing(unittest.TestCase):   
    def test_bessel(self): 
        n = 1000
        ds = math.pi/n
        ib0 = 0
        ib1 = 0
        x= 1.4
        for j in range (0, n):
            tt = (j - 0.5)*ds
            ib0 += math.exp(x*math.cos(tt))*ds 
            ib1 += +math.exp(x*math.cos(tt))*ds*math.cos(tt)
        ib0 = ib0/math.pi
        ib1 = ib1/math.pi
        kb0 = 0
        kb1 = 0
        ds = 2*math.pi/n
        for j in range(0, n):
            tt = (j - 0.5)*ds
            ch = (math.exp(tt) + math.exp(-tt))/2
            kb1 += math.exp(-x*ch)*ch*ds
            kb0 += math.exp(-x*ch)*ds
        i0 = sp.iv(0, x)
        i1 = sp.iv(1, x)
        k0 = sp.kv(0, x)
        k1 = sp.kv(1, x)
        self.assertTrue(((abs(ib0 - i0)<0.005)&(abs(ib1 - i1)<0.005)&(abs(kb0 - k0)<0.005)&(abs(kb1 - k1)<0.005)), msg = 'Deltas I0, I1, K0, K1: ' + str(abs(i0 - ib0)) + ', ' + str(abs(i1 - ib1)) + ', ' + str(abs(k0 - kb0)) + ', ' + str(abs(k1 - kb1)))

    def test_integrate(self):
            nfi = 1000
            asumbar = 0
            age = 250
            dfi = 2*math.pi / nfi
            for j in range (0, nfi, 1):
                fi = (j-0.5)*dfi
                x = math.cos(fi)
                y = math.sin(fi)
                axb, ayb = rings.galaxy.forcebar(x, y, 0, 1, age)
                asumbar = asumbar+axb*math.cos(fi)+ayb*math.sin(fi)
            abar1 = asumbar / nfi
            def abar(fi):
                ax, ay = rings.galaxy.forcebar(math.cos(fi), math.sin(fi), 0, 1, age)
                abart = math.cos(fi) + math.sin(fi)
                return(abart)
            abar2 = integrate.quad(lambda fi: abar(fi), 0, 2*math.pi)[0]
            self.assertTrue(abs(abar1 - abar2)<0.05 , msg = 'Delta a_Bar is: '+str(abar1 - abar2))

    def test_time(self):
        part = 10000
        age = 250
        ncol = 0
        Rc, Vtot, Vdisk, Vbar = rings.wwd.read_data()
        init = rings.galaxy.set_init(part, age)
        data = rings.galaxy.first_step(init, age)
        t2 = time()
        for k in range(0, 3):
            rings.galaxy.main_cycle(data, k, age, ncol)
        t3 = time()
        t_ord = t3 - t2
        t0 = time()
        for k in range(0, 3):
            Process(target = rings.galaxy.main_cycle, args = [data, k, age, ncol]).start()
        t1 = time()
        t_mult = t1 - t0
        self.assertTrue(t_mult < t_ord, msg = 'Время вычисления с использованием multiprocessing больше времени простого вычисления, эксперимент, очевидно, неудачный.')
