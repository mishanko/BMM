import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


Nx = 160 # число шагов по x
Ny = 160 # число шагов по y
M = 100 #число шагов по t

x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
t = np.linspace(0, 1, M)

h_x = x[1] - x[0]
h_y = y[1] - y[0]
tau = t[1] - t[0]

gamma_x = tau / h_x**2
gamma_y = tau / h_y**2

def lam(k):
    return (np.pi ** 2) * (1 + (2 * k + 1) ** 2)

def exact(x, y, t):
    z = 0
    for k in range(10):
        z = z + (np.sin((np.pi * x) * (2 * k + 1)) * (np.exp(t) - np.exp(-lam (k) * t))) / ((2 * k + 1) * np.pi * (lam(k) + 1))
    return 4 * np.cos(np.pi * y) * z


u = np.zeros((Nx, Ny, 2 * M + 1))

def F_1(i1, i2, j):
    return 0.5 * gamma_y * (u[i1, i2-1, j-1] + u[i1, i2+1, j-1]) + (1 - gamma_y) * u[i1, i2, j-1] + 0.5 * tau * np.exp(tau * (j + 1 / 2)) * np.cos(np.pi * y[i2])

def F_2(i1, i2, j):
    return 0.5 * gamma_x * (u[i1-1, i2, j-1] + u[i1+1, i2, j-1]) + (1 - gamma_x) * u[i1, i2, j-1] + 0.5 * tau * np.exp(tau * (j + 1 / 2)) * np.cos(np.pi * y[i2])

def progonka_x(i2, j):

    d = np.zeros(Nx)
    sigma = np.zeros(Nx)
    d[1] = 0
    sigma[1] = 0
    A = 0.5 * gamma_x
    B = 1 + gamma_x
    C = 0.5 * gamma_x

    for m in range(1, Nx-1):
        Fm = -F_1(m, i2, j)
        d[m+1] = C / (B - A * d[m])
        sigma[m+1] = (Fm - A * sigma[m]) / (A * d[m] - B)
    u[Nx-1, i2, j] = 0
    u[0, i2, j] = 0

    for m in range(Nx-1, 0, -1):
        u[m-1, i2, j] = d[m] * u[m, i2, j] + sigma[m]

def progonka_y(i1, j):

    d = np.zeros(Ny)
    sigma = np.zeros(Ny)
    d[1] = 1
    sigma[1] = 0
    A = 0.5 * gamma_y
    B = 1 + gamma_y
    C = 0.5 * gamma_y
    
    for m in range(1, Ny-1):
        Fm = -F_2(i1, m, j)
        d[m+1] = C / (B - A * d[m])
        sigma[m+1] = (Fm - A * sigma[m]) / (A * d[m] - B)
    u[i1, Ny - 1, j] = sigma[-1] / (1 - d[-1])

    for m in range(Ny-1, 0, -1):
        u[i1, m-1, j] = d[m] * u[i1, m, j] + sigma[m]

u[:, :, 0] = 0 #начальное условие

for j in tqdm(range(1, 2 * M, 2)): #выполнение расчетов
    for i2 in range(1, Ny - 1):
        progonka_x(i2, j)
    for i1 in range(1, Nx - 1):
        progonka_y(i1, j+1)

j_ = 80

p = np.zeros((Nx, Ny))
def p_(j_):
    for i in tqdm(range(Nx)):
        for j in range(Ny):
            p[i][j] = exact(i * h_x, j * h_y, j_ * tau)
    return p

ym = np.linspace(0, 1, num=Ny)
xn = np.linspace(0, 1, num=Nx)
X, Y = np.meshgrid(xn, ym)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
print(p.shape, X.shape, Y.shape)
nnn = p_(j_)[:, :]
surf_1 = ax.plot_wireframe(X, Y, nnn, rstride=10, cstride=1, color='r')
surf_2 = ax.plot_wireframe(X, Y, u[:, :, j_], rstride=10, cstride=1)
surf_3 = ax.plot_wireframe(X, Y, np.abs(nnn - u[:, :, j_]), rstride=10, cstride=1, color='g')
print(np.mean(np.abs(nnn - u[:, :, j_])))

plt.title('Решение')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

