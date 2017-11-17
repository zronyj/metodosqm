# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:17:34 2017

@author: Rony J. Letona
@email: zronyj@gmail.com

Descripcion
-----------
Pequeno script para realizar el calculo de la energia y forma de la funcion de
onda dela particula en una caja.

Referencias
-----------

(1) Mueller, M. (2002). Fundamentals of Quantum Chemistry (1ra ed). New York:
    Kluwer Academic Publishers.

(2) Piela, L. (2007). Ideas of Quantum Chemistry (1ra ed). Oxford: Elsevier.

"""

import numpy as np
import scipy.constants as cnts
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# *******************************************
# Funcion de onda de la particula en una caja
# Referencia (2) p. 146 eq. (4.6)
# *******************************************
def psi(x, n = 1, L = 1, m = cnts.m_e):
    N = (2 / L)**0.5
    f = np.sin(n * np.pi * x / L)
    return N * f

# *************************************
# Energia del oscilador armonico para n
# Referencia (2) p. 146 eq. (4.5)
# *************************************
def energia(n = 1, L = 1, m = cnts.m_e):
    return n**2 * cnts.h**2 / (8 * m * L**2)

# ****************************************
# Funcion para graficar la funcion de onda
# de la particula en una caja
# ****************************************
def orbital(n = 1, L = 1, m = cnts.m_e):
    x = np.arange(0, L, 0.001)
    vf = np.vectorize(psi)
    y = vf(x, n=n, L=L, m=m)
    y2 = vf(x, n=n, L=L, m=m)**2
    plt.title("Funcion de onda de la particula en una caja")
    plt.plot(x, y, "r--", label="$\Psi$")
    plt.plot(x, y2, "b-", label="$\Psi^2$")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.grid(True)
    plt.show()

# ****************************************
# Funcion para graficar la funcion de onda
# del oscilador armonico cuantico en 2
# dimensiones nx, ny = 0, 1, ... 10
# ****************************************
def orbital2D(nx = 2, ny = 2, Lx = 1, Ly = 1, mx = cnts.m_e,
              my = cnts.m_e):
    vf = np.vectorize(psi)
    x = np.arange(0, Lx, 0.005)
    y = np.arange(0, Ly, 0.005)
    X, Y = np.meshgrid(x, y)
    Z = vf(X, n=nx, L=Lx, m=mx) * vf(Y, n=ny, L=Ly, m=my)
    plt.contourf(x, y, Z, 50, cmap=cm.bone)
    plt.show()