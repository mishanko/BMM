import matplotlib.pyplot as plt
import numpy as np
import math as mt
from mpl_toolkits.mplot3d import axes3d
from tqdm import tqdm
import sys

# * Объявим некоторые 
L, T = mt.pi, mt.pi # * x in (0, L) и длина расчетного промежутка по времени 
H = 3 # * y in (0, H) 
N = 100 # * число узлов по x 
M = 100 # * число узлов по y 

hx = L / (N-1) # * шаг по x 
hy = H / (M-1) # * шаг по y 

# J = 2 + round(2 * T * mt.sqrt(1 / hx ** 2 + 1 / hy ** 2))
J = 100

tau =  T / (J-1) # * шаг по t 

# * Задание сетки---------------------------------------------------- 

x = np.linspace(0, L, N-1)
y = np.linspace(0, H, M-1)
t = np.linspace(0, T, J-1)

# * ------------------------------------------------------------------ 
# *  Аналитическое решение-------------------------------------------- 

u = np.zeros((N, M, J))
for j in range(J-1):
    for n in range(N-1):
        for m in range(M-1):
            u[n, m, j] = (np.exp(-4 * t[j]) + 1/17 * 
                          (4 * np.sin(t[j]) - np.cos(t[j]) + np.exp(-4 * t[j]))) * \
                            np.cos(x[n]) * np.sin(y[m] * mt.pi / 3)

# * ----------------------------------------------------------------- 
# * Численное решение------------------------------------------------ 

v, err = np.zeros((N, M, J)), np.zeros((N, M, J)) # * численное решение и погрешность 

# * Начальные условия:
for n in range(N-1):
    for m in range(M-1):
        v[n, m, 0] = np.cos(x[n])

gamma_x = 4 * tau / hx ** 2
gamma_y = 4 * tau / hy ** 2

def F_1(n, m, j):
    return 0.5 * gamma_y * (v[n, m-1, j-1] + v[n, m+1, j-1]) + (1 - gamma_y) * v[n, m, j-1] + 2 * tau * np.cos(x[n]) * np.sin(tau*((j + 1) / 2))

def F_2(n, m, j):
    return 0.5 * gamma_x * (v[n-1, m, j-1] + v[n+1, m, j-1]) + (1 - gamma_x) * v[n, m, j-1] + 2 * tau * np.cos(x[n]) * np.sin(tau*((j + 1) / 2))

def progonka_x(m, j):
    d = np.zeros(N)
    sigma = np.zeros(N)
    d[1] = 0
    sigma[1] = 1
    A = 0.5 * gamma_x
    B = 1 + gamma_x
    C = 0.5 * gamma_x

    for n in range(1, N-1):
        Fm = -F_1(n, m, j)
        d[n+1] = C / (B - A * d[n])
        sigma[n+1] = (Fm + A * sigma[n]) / (B - A * d[n])

    v[N-1, m, j] = sigma[N-2] / (1 - d[N-2])
    
    for n in range(N-1, 0, -1):
        v[n-1, m, j] = d[n] * v[n, m, j] + sigma[n]
        err[n-1, m, j] = v[n-1, m, j] - u[n-1, m, j]
        
def progonka_y(n, j):
    d = np.zeros(M)
    sigma = np.zeros(M)
    d[1] = 0
    sigma[1] = 1
    A = 0.5 * gamma_y
    B = 1 + gamma_y
    C = 0.5 * gamma_y

    for m in range(1, M-1):
        Fm = -F_2(n, m, j)
        d[m+1] = C / (B - A * d[m])
        sigma[m+1] = (Fm + A * sigma[m]) / (B - A * d[m])
        
    v[n, M-1, j] = sigma[M-2] / (1 - d[M-2])
    err[n, M-1, j] = v[n, M-1, j] - u[n, M-1, j]
    for m in range(M-1, 0, -1):
        v[n, m-1, j] = d[m] * v[n, m, j] + sigma[m]
        err[n, m-1, j] = v[n, m-1, j] - u[n, m-1, j]


for j in tqdm(range(1, J-1)): #выполнение расчетов
    for m in range(1, M-1):
        progonka_x(m, j)
    for n in range(1, N-1):
        progonka_y(n, j+1)
        
    # * граничные условия:
    for n in range(N):
        v[n, 1, j] = v[n, 0, j]
        err[n, 0, j] = v[n, 0, j] - u[n, 0, j]
        v[n, M-2, j] = v[n, M-1, j]
        err[n, M-1, j] = v[n, M-1, j] - u[n, M-1, j]

    for m in range(M):
        v[1, m, j] = v[0, m, j]
        err[0, m, j] = v[0, m, j] - u[0, m, j]
        v[N-2, m, j] = v[N-1, m, j] 
        err[N-1, m, j] = v[N-1, m, j] - u[N-1, m, j]

j_ = 10

ym = np.linspace(0, H, num=M)
xn = np.linspace(0, L, num=N)
X, Y = np.meshgrid(xn, ym)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# print(v.shape, X.shape, Y.shape, v[:, :, j_].shape)
surf_1 = ax.plot_wireframe(Y, X, v[:, :, j_], rstride=10, cstride=1)
surf_2 = ax.plot_wireframe(Y, X, u[:, :, j_], rstride=10, cstride=1, color='r')

plt.title('Решение')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

