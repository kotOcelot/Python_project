from __future__ import unicode_literals
import rings.galaxy
import rings.graphics
import rings.curve
import rings.wwd


def main():
    print('Программа моделирует движение газовых облаков в потенциале Галактики и работает КРАЙНЕ медленно.')            
    print('Поскольку для преподавателя выполнять полный расчёт нецелесообразно, а с уменьшением шага или числа частиц теряется структура полученного изображения, готовые результаты расчёта находятся в папке /data/')
    print('Тем не менее, чтобы продемонстрировать работоспособность кода, проведём расчёт нескольких итераций.')
    print('Для этого сначала необходимо построить кривую вращения Галактики. Результат и картинка сохраняется в папку /results/. Начать? (y/n)')

    ans = input()
    if (ans == 'y'):
        rings.curve.rotcurve()
    if (ans == 'n'):
        exit()

    print('Теперь проведём три итерации основного цикла вычислений.')
    
    part = 10000
    ncol = 0
    dt = 0.005
    age = 400
    print('Для начала, зададим случайные начальные скорости и координаты облаков для T = 0 Myr.')
    init = rings.galaxy.set_init(part, age)
    print('Массив создан. Сохраняем его.')
    rings.wwd.write_init(init, 0)
    print('Первый шаг моделирования выполнен.')
    data = rings.galaxy.first_step(init, age)
    print('Переходим к основному циклу.')    
    for k in range (0, 3):
        print('Итерация ' + str(k) + '.')
        rings.galaxy.main_cycle(data, k, age, ncol)
        if ((k // 1000)*1000 == k or k == 1):
            kt = round(k*dt + 5)
            rings.wwd.write_res(data, kt)
            
    print('Данные записаны. Теперь построим изображения полученных распределений газовых облаков: при T = 0 и T = 5 Myr.')
    print('Это займёт некоторое время. Начать? (y/n)')
    
    ans1 = input()
    if (ans1 == 'y'):
        rings.graphics.plot_distr(0, 'results')
        print('Распределение для T = 0 Myr построено. Строим следующее.')
        rings.graphics.plot_distr(5, 'results')
    if (ans1 == 'n'):
        exit()
        
    print('Изображения сохранены в /results/.')
    print('Построим распределение частиц из готовых данных. Красным отметим OB - ассоциации, образованные в результате столкновений газовых облаков.')
    print('Введите время от начала моделирования в Myr от 0 до 2000 с шагом 5 (рекомендовано: 1950):')
    
    ans2 = input()
    rings.graphics.plot_distr(int(ans2), 'data')
    print('Изображения сохранены в /results/.')
    print('Наконец, построим орбиту одной точки. Введите начальную координату этой точки в кпк (от 0.03 до 12, рекомендовано: 8). Осторожно, долго!')
    
    ans3 = input()
    print('Начать? (y/n)')
    ans4 = input()
    if (ans4 == 'y'):
        rings.graphics.plot_orbit(float(ans3), 'data')
        print('Изображения сохранены в /results/.')  
    if (ans4 == 'n'):
        exit()  
    print('На этом всё.')

if __name__ == '__main__':
    main()
    

            