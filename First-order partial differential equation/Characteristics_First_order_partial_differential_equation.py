
# Импортируем необходимые библиотеки
from math import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d


# Определим функции характеристик через массивы, где итерирование
# будет идти по соответствующему неизвестному параметру: x0
# или t0. Здесь x0 и t0
# взяты с определенными шагами c помощью np.arange.

def family_1(x):
    return [((1 + cos(0.5 * pi * x0)) * (x0 - x)) for x0 in np.arange(-1, 0.1, 0.1)]


def family_2(x):
    return [((1 + exp(-t0)) * (-x) + t0) for t0 in np.arange(0, 1.1, 0.1)]


x_list = np.arange(-1, 0.1, 0.1)

# Создадим массив значений по x от −1 до 0
# с определенным шагом и соответсвущие массивы
# для функций с итерированием уже по x

x_list = np.arange(-1, 0.1, 0.1)
family_1_list = [family_1(x) for x in x_list]
family_2_list = [family_2(x) for x in x_list]

# Построим соответствующие графики.

plt.subplot(1, 2, 1)
plt.ylim(0, 4)
plt.xlim(-1, 0)
plt.plot(x_list, family_1_list)
plt.title('Характеристики', loc='right')
plt.ylabel('t')
plt.xlabel('x')
plt.subplot(1, 2, 2)
plt.ylim(0, 2)
plt.xlim(-1, 0)
plt.plot(x_list, family_2_list)
plt.ylabel('t')
plt.xlabel('x')

plt.show()
