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

hx = L / N # * шаг по x 
hy = H / M # * шаг по y 

# J = 2 + round(2 * T * mt.sqrt(1 / hx ** 2 + 1 / hy ** 2))
J = 100

tau =  T / J # * шаг по t 
x = np.zeros((N, )) # * массив узлов по x 
y = np.zeros((M, )) # * массив узлов по y 
t = np.zeros((J, )) # * массив узлов по t 

# * Задание сетки---------------------------------------------------- 

for n in range(N):
    x[n] = n * hx / 2 - (n + 1) * hx
for m in range(M):
    y[m] = m * hy / 2 - (m + 1) * hy
for j in range(J):
    t[j] = j * tau

# * ------------------------------------------------------------------ 
# *  Аналитическое решение-------------------------------------------- 

u = np.zeros((M, N, J))
for j in range(J):
    for n in range(N):
        for m in range(M):
            u[m, n, j] = (np.exp(-4 * t[j]) + 1/17 * (4 * np.sin(t[j]) - np.cos(t[j]) + np.exp(-4 * t[j]))) * np.cos(x[n]) * np.cos(y[m]/3)


# * ----------------------------------------------------------------- 
# * Численное решение------------------------------------------------ 

v, err, Fx, Fy = np.zeros((M, N, J)), np.zeros((M, N, J)), np.zeros((M, N, J)), np.zeros((M, N, J)) # * численное решение и погрешность 

# * Начальные условия:
for n in range(N):
    for m in range(M):
        v[m, n, 0] = np.cos(x[n])

alpha_x=np.zeros((M, N)) 
beta_x=np.zeros((M, N)) 
alpha_y=np.zeros((M, N)) 
beta_y=np.zeros((M, N)) 
Ay = 2 * tau / hy**2
By = 2 * tau / hy**2
Ax = 2 * tau / hx**2
Bx = 2 * tau / hx**2
Cx = (1 + 4 * tau / hx**2) 
Cy = (1 + 4 * tau / hy**2)

for j in tqdm(range(J-1)):
    # * граничные условия:
    for m in range(M):
        v[m, 0, j] = v[m, 1, j]
        err[m, 0, j] = v[m, 0, j] - u[m, 0, j]
        v[m, N-1, j] = v[m, N-2, j]
        err[m, N-1, j] = v[m, N-1, j] - u[m, N-1, j]

    for n in range(N):
        v[0, n, j] = v[1, n, j]
        err[0, n, j] = v[0, n, j] - u[0, n, j]
        v[M-1, n, j] = v[M-2, n, j] 
        err[M-1, n, j] = v[M-1, n, j] - u[M-1, n, j]
        
    # * прямой ход прогонки по x---------------------------------     
    for n in range(1, N-1):
        for m in range(1, M-1):
            Fx[m, n, j] = 1 / 2 * (v[m-1, n, j] - 2 * v[m, n, j] + v[m+1, n, j]) * 4 * tau / hx**2 + v[m, n, j] * 4 + 2 * tau * np.cos(x[n]) * np.sin(t[j] + tau / 2)

    for m in range(1, M-1):   
        alpha_x[m, 1] = 1
        beta_x[m, 1] = 0   
        for n in range(1, N-1):
            alpha_x[m, n+1] = Bx / (Cx - Ax * alpha_x[m, n])
            beta_x[m, n+1] = (Fx[m, n, j] + Ax * beta_x[m, n]) / (Cx - Ax * alpha_x[m, n])

    # * обратный ход прогонки по x-------------------------------
    for m in range(1, M-1):
        v[m, N-1, j] = v[m, N-2, j]
        for n in range(N-2, -1, -1):
            v[m, n, j] = alpha_x[m, n+1] * v[m, n+1, j] + beta_x[m, n+1]
            err[m, n, j] = v[m, n, j] - u[m, n, j]

    # ------------------------------------------------------------------------------

    # * прямой ход прогонки по y---------------------------------     
    for m in range(1, M-1):
        for n in range(1, N-1):
            Fy[m, n, j] = 1 / 2 * (v[m, n-1, j] - 2 * v[m, n, j] + v[m, n+1, j]) * 4 * tau / hy**2 + v[m, n, j] * 4 + 2 * tau * np.cos(x[n]) * np.sin(t[j] + tau / 2)

    for n in range(1, N-1):   
        alpha_y[1, n] = 1
        beta_y[1, n] = 0   
        for m in range(1, M-1):
            alpha_y[m+1, n] = By / (Cy - Ay * alpha_y[m, n])
            beta_y[m+1, n] = (Fy[m, n, j] + Ay * beta_y[m, n]) / (Cy - Ay * alpha_y[m, n])

    # * обратный ход прогонки по y-------------------------------
    for n in range(1, N-1):
        v[M-1, n, j] = v[M-2, n, j]
        for m in range(M-2, -1, -1):
            v[m, n, j] = alpha_y[m+1, n] * v[m+1, n, j] + beta_y[m+1, n]
            err[m, n, j] = v[m, n, j] - u[m, n, j]

    

    

j_ = 30

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

