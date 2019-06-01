#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 11:00:54 2019

@author:  Mikhailov Mikhail
"""

import matplotlib.pyplot as plt
import numpy as np
import math as mt
from mpl_toolkits.mplot3d import axes3d

# Зададим: ϵ - точность в методе Ньютона, N - количество шагов по x,
# М - количество шагов по y , а также границы нашей сетки.

epsilon = 0.00001
N = 100
M = 100
T_begin = 0
T_end = 1
X_begin = 0
X_end = -1

# Соответственно, элементарные шаги.

h_x = (X_end - X_begin) / (N - 1)
h_t = (T_end - T_begin) / (M - 1)

# Создадим двумерный массив размерами с нашу сетку (N×M),
# в ячейках которого будут храниться соответствующие искомые значения.

y = np.zeros((M, N))

# Начнем заполнять его начальным и граничным значениями.

for n in range(N):
    y[0][n] = (mt.cos(mt.pi * h_x * n * 0.5))

for m in range(M):
    y[m][0] = mt.exp(-h_t * m)


# Определим вспомогательные функции.

def F(m, n):
    return mt.log(y[m][n] + 1)


def df(mp1, np1):
    return (1 / (2 * h_t) - 0.5 / (h_x * (y[mp1][np1] + 1)))


# Разностная схема будет иметь вид.

def f(mp1, np1):
    n = np1 - 1
    m = mp1 - 1
    de = (y[mp1][n] - y[m][n] + y[mp1][np1] - y[m][np1]) / (2. * h_t) - (
            F(mp1, np1) - F(mp1, n) + F(m, np1) - F(m, n)) / (2. * h_x)
    return (de)


# Перейдем к методу Ньютона, пробегая по всей сетке.

eps = epsilon + 1
while eps > epsilon:
    eps = 0
    for m in range(M)[0:M - 1]:
        for n in range(N)[0:N - 1]:
            ep = f(m + 1, n + 1) / df(m + 1, n + 1)
            y[m + 1][n + 1] = y[m + 1][n + 1] - ep
            if abs(ep) > eps:
                eps = abs(ep)

            # Построим график решения.

tm = np.linspace(T_begin, T_end, num=M)
xn = np.linspace(X_begin, X_end, num=N)
X, T = np.meshgrid(xn, tm)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf_1 = ax.plot_wireframe(X, T, y, rstride=10, cstride=1)
plt.title('Решение')
plt.xlabel('X')
plt.ylabel('T')
plt.show()

# Погрешности

print('y=', y[5][5])
print('yтеор=', mt.cos(mt.pi * xn[5] / 2) + mt.exp(-tm[5]) - 1)
print(y[5][5] - (mt.cos(mt.pi * xn[5] / 2) + mt.exp(-tm[5]) - 1))
