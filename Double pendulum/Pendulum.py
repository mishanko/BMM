# Двойной маятник. Анализ проводится с точки зрения лагранжевой механики.
# Лагранжевы переменные - это угол верхнего стержня, угол нижнего стержня, измеренный от вертикали.

# Импортируем библиотеку для визуализации
from vpython import *

# Задаем параметры сцены
scene.width = scene.height = 600
scene.range = 1.8
scene.title = "A double pendulum"

# Начальные ускорение и скорость первого матника
add1 = 0
ad1 = 0

# Начальный угол отклонения первого маятника
a1 = pi

# Начальные ускорение и скорость второго матника
add2 = 0
ad2 = 0

# Начальный угол отклонения второго маятника
a2 = pi / 2

# Начальные координаты первого и второго матника
x1 = 0
x2 = 0
y1 = 0
y2 = 0

# Массы маятников 1 и 2 соответственно
M1 = 10
M2 = 10

# Длины стержней 1 и 2 соответственно
L1 = 0.5
L2 = 0.5

# Ускорение свободного падения
g = 10.4

# Величина изменения времени
dt = 0.0002
t = 0
d = 0.05  # толщина палочки
gap = 2 * d  # расстояние между двумя верхними стержнями

L1display = L1 + d  # чтобы палочка не выпирал из краев верхних стержней
L2display = L2 + d / 2  # чтобы палочка не выпирала из краев нижнего стержня

hpedestal = 1.3 * (L1 + L2)  # высота пьедестала
wpedestal = 0.1  # ширина пьедестала
tbase = 0.05  # толщина основания
wbase = 100 * gap  # ширина основания
offset = 10 * gap  # длина верхней палочки
pedestal_top = vec(0, hpedestal / 2, 0)  # чтобы картинка была по центру и в центре пьедестала

pedestal = box(pos=pedestal_top - vec(0, hpedestal / 2, offset),
               size=vec(wpedestal, 1.1 * hpedestal, wpedestal),
               color=vec(0.4, 0.4, 0.5))  # форма пьедестала, размер и цвет
base = box(pos=pedestal_top - vec(0, hpedestal + tbase / 2, offset),
           size=vec(wbase, tbase, wbase),
           color=pedestal.color)  # форма основания, размер и цвет
axle1 = cylinder(pos=pedestal_top - vec(0, 0, gap / 2 - d / 4), axis=vec(0, 0, -1),
                 size=vec(offset, d / 4, d / 4), color=color.yellow)  # форма палочки1, размер и цвет

bar1 = box(pos=pedestal_top + vec(L1display / 2 - d / 2, 0, -(gap + d) / 2),
           size=vec(L1display, d, d), color=color.blue)  # форма стержня1, размер и цвет

bar1b = box(pos=pedestal_top + vec(L1display / 2 - d / 2, 0, (gap + d) / 2),
            size=vec(L1display, d, d), color=bar1.color)  # форма стержня11, размер и цвет


# Реализация вращения
bar1.rotate(angle=-pi / 2, axis=vec(0, 0, 1), origin=vec(axle1.pos.x, axle1.pos.y, bar1.pos.z))
bar1.rotate(angle=a1, axis=vec(0, 0, 1), origin=vec(axle1.pos.x, axle1.pos.y, bar1.pos.z))

bar1b.rotate(angle=-pi / 2, axis=vec(0, 0, 1), origin=vec(axle1.pos.x, axle1.pos.y, bar1b.pos.z))
bar1b.rotate(angle=a1, axis=vec(0, 0, 1), origin=vec(axle1.pos.x, axle1.pos.y, bar1b.pos.z))

pivot1 = vec(axle1.pos.x, axle1.pos.y, 0)  # начальное положение палочки 1

axle2 = cylinder(pos=pedestal_top + vec(L1, 0, -(gap + d) / 2), axis=vec(0, 0, 1),
                 size=vec(gap + d, axle1.size.y / 2, axle1.size.y / 2), color=axle1.color)  # форма палочки2, размер и цвет

axle2.rotate(angle=-pi / 2, axis=vec(0, 0, 1), origin=vec(axle1.pos.x, axle1.pos.y, axle2.pos.z))
axle2.rotate(angle=a1, axis=vec(0, 0, 1), origin=vec(axle1.pos.x, axle1.pos.y, axle2.pos.z))

bar2 = box(pos=axle2.pos + vec(L2display / 2 - d / 2, 0, (gap + d) / 2),
           size=vec(L2display, d, d), color=color.blue) # форма стержня2, размер и цвет

bar2.rotate(angle=-pi / 2, axis=vec(0, 0, 1), origin=vec(axle2.pos.x, axle2.pos.y, bar2.pos.z))
bar2.rotate(angle=a2, axis=vec(0, 0, 1), origin=vec(axle2.pos.x, axle2.pos.y, bar2.pos.z))

while True:
    rate(1 / dt)
    # Решение уравнения Лагранжа
    add1 = (M2 * L2 * add2 * cos(a1 - a2) + M2 * L2 * ad2 * ad2 * sin(a1 - a2) + (M1 + M2) * g * sin(a1)) / (
                -(M1 + M2) * L1)
    add2 = (L1 * add1 * cos(a1 - a2) - L1 * ad1 * ad1 * sin(a1 - a2) + g * sin(a2)) / (-L2)
    ad1 = ad1 + add1 * dt
    a1 = a1 + ad1 * dt
    da1 = ad1 * dt
    ad2 = ad2 + add2 * dt
    da2 = ad2 * dt
    a2 = a2 + ad2 * dt
    x1 = L1 * sin(a1)
    x2 = L1 * sin(a1) + L2 * sin(a2)
    y1 = -L1 * cos(a1)
    y2 = -L1 * cos(a1) - L2 * cos(a2)
    t = t + dt

    bar1.rotate(angle=da1, axis=vec(0, 0, 1), origin=pivot1)
    bar1b.rotate(angle=da1, axis=vec(0, 0, 1), origin=pivot1)
    pivot2 = vec(axle2.pos.x, axle2.pos.y, pivot1.z)  # начальное положение палочки 2
    axle2.rotate(angle=da1, axis=vec(0, 0, 1), origin=pivot1)
    bar2.rotate(angle=da2, axis=vec(0, 0, 1), origin=pivot2)
    bar2.pos = pivot2 + bar2.axis / 2  # чтобы нижная палочка крепился чуть ниже топа
