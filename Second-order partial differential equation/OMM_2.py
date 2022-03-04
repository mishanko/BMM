import matplotlib.pyplot as plt
import numpy as np
import math as mt
from mpl_toolkits.mplot3d import axes3d
from tqdm import tqdm
import sys

# * Объявим некоторые 
L, T = mt.pi, mt.pi # * x in (0, L) и длина расчетного промежутка по времени 
H = 3 # * y in (0, H) 
N = 101 # * число узлов по x 
M = 51 # * число узлов по y 

hx = L / (N - 1) # * шаг по x 
hy = H / (M - 1) # * шаг по y 

# J = 2 + round(2 * T * mt.sqrt(1 / hx ** 2 + 1 / hy ** 2))
J = 101

tau =  T / (J - 1) # * шаг по t 
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
            u[m, n, j] = (np.exp(-4 * t[j]) + (1 / 17) * (4 * np.sin(t[j]) - np.cos(t[j]) + np.exp(-4 * t[j]))) * np.cos(x[n]) * np.cos(y[m]/3)


# * ----------------------------------------------------------------- 
# * Численное решение------------------------------------------------ 

v, err = np.zeros((M, N, J)), np.zeros((M, N, J)) # * численное решение и погрешность 

# * Начальные условия:
for n in range(N):
    for m in range(M):
        v[m, n, 0] = np.cos(x[n])

w = np.zeros((M, N)) # * вспомогательная функция при переходе со слоя на слой 
alpha_x=np.zeros((N-1, )) 
beta_x=np.zeros((N-1, )) 
alpha_y=np.zeros((M-1, )) 
beta_y=np.zeros((M-1, )) 
Ay = 2 * tau / hy**2
By = 2 * tau / hy**2
Ax = 2 * tau / hx**2
Bx = 2 * tau / hx**2
Cx = (1 + 4 * tau / hx**2) 
Cy = (1 + 4 * tau / hy**2)

for j in tqdm(range(J-1)):
    # * вычисление вспомогательной функции w ?----------------------- 
    for m in range(1, M-1):
        alpha_x[0] = 1
        # * прямой ход прогонки по x--------------------------------- 
        for n in range(1, N-1):
            F = 1 / 2 * (v[m-1, n, j] - 2 * v[m, n, j] + v[m+1, n, j]) * 4 * tau / hx**2 + (1 - 4 * tau / hx**2) * v[m, n, j] + 2 * tau * np.cos(x[n]) * np.sin(t[j] + tau / 2)
            alpha_x[n] = Bx / (Cx - Ax * alpha_x[n-1])
            beta_x[n] = (F + Ax * beta_x[n-1]) / (Cx - Ax * alpha_x[n-1])

        # * обратный ход прогонки по x-------------------------------
        w[m, N-1] = beta_x[N-2] / (1 - alpha_x[N-2])
        for n in range(N-2, -1, -1):
            w[m, n-1] = alpha_x[n] * w[m, n] + beta_x[n]

        # w[m, N-1] = w[m, N-2]
        # w[m, 0] = w[m, 1]
        # * вычисление функции v на новом слое по времени------------ 
        alpha_y[0] = 1

    for n in range(1, N-1):
        
        # * прямой ход прогонки по y--------------------------------- 
        for m in range(1, M-1):
            F = 1/2 * (w[m, n-1] - 2 * w[m, n] + w[m, n+1]) * 4 * tau / hy**2 + (1 - 4 * tau / hy**2) * w[m, n] + 2 * tau * np.cos(x[n]) * np.sin(t[j] + tau / 2) 
            alpha_y[m] = By / (Cy - Ay * alpha_y[m-1])
            beta_y[m] = (F + Ay * beta_y[m-1]) / (Cy - Ay * alpha_y[m-1])
            
        # * обратный ход прогонки по y-------------------------------
        v[M-1, n, j+1] = beta_y[M-2] / (1 - alpha_y[M-2])
        err[M-1, n, j+1] = v[M-1, n, j+1] - u[M-1, n, j+1]
        for m in range(M-2, -1, -1):
            v[m-1, n, j+1] = alpha_y[m] * v[m, n, j+1] + beta_y[m]
            err[m, n, j+1] = v[m, n, j+1] - u[m, n, j+1]

    # * граничные условия:
    for m in range(M):
        v[m, 0, j+1] = v[m, 1, j+1]
        err[m, 0, j+1] = v[m, 0, j+1] - u[m, 0, j+1]
        v[m, N-1, j+1] = v[m, N-2, j+1]
        err[m, N-1, j+1] = v[m, N-1, j+1] - u[m, N-1, j+1]

    for n in range(N):
        v[0, n, j+1] = v[1, n, j+1]
        err[0, n, j+1] = v[0, n, j+1] - u[0, n, j+1]
        v[M-1, n, j+1] = v[M-2, n, j+1] 
        err[M-1, n, j+1] = v[M-1, n, j+1] - u[M-1, n, j+1]

j_ = 20

ym = np.linspace(0, H, num=M)
xn = np.linspace(0, L, num=N)
X, Y = np.meshgrid(xn, ym)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
print(v.shape, X.shape, Y.shape, v[:, :, j_].shape)
surf_1 = ax.plot_wireframe(X, Y, v[:, :, j_], rstride=10, cstride=1)
surf_2 = ax.plot_wireframe(X, Y, u[:, :, j_], rstride=10, cstride=1, color='r')

plt.title('Решение')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

