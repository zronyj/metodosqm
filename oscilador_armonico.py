# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 20:38:52 2017

@author: Rony J. Letona
@email: zronyj@gmail.com

Descripcion
-----------
Pequeno script para realizar el calculo de la energia y forma de la funcion de
onda del oscilador armonico.

Referencias
-----------

(1) Lowe, J. P., Peterson, K. A. (2006). Quantum Chemistry (3ra ed). London:
    Elsevier.

(2) Piela, L. (2007). Ideas of Quantum Chemistry (1ra ed). Oxford: Elsevier.

"""
import numpy as np
import scipy.special as spe
import scipy.constants as cnts
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Funcion para calcular factorial n!
factorial = lambda n: np.prod( np.array( [i for i in range(1,n)] ) )

# **************************************
# Funcion de onda del oscilador armonico
# Referencia (2) p. 165 eq. (4.16)
# **************************************
def psi(x, n = 1, f = 0.015, m = cnts.m_e):
    omega = 2 * np.pi * f
    a = m * omega / cnts.hbar
    y = a**0.5 * x
    cts = (a / np.pi)**0.25 * 1 / (2**n * factorial(n))**0.5
    poly_herm = spe.eval_hermite(n, y)
    return cts * poly_herm * np.exp(-y**2 / 2)

# *************************************
# Energia del oscilador armonico para n
# Referencia (2) p. 165 eq. (4.17)
# *************************************
def energia(n = 1, f = 0.015):
    omega = 2 * np.pi * f
    return cnts.hbar * omega * (n + 0.5)

# ****************************************
# Funcion para graficar la funcion de onda
# del oscilador armonico cuantico n = 0-10
# ****************************************
def orbital(n = 0, f = 0.015, m = cnts.m_e, d = [-0.2, 0.2, 0.0001]):
    x = np.arange(d[0], d[1], d[2])
    vf = np.vectorize(psi)
    y = vf(x, n=n, f=f, m=m)
    y2 = vf(x, n=n, f=f, m=m)**2
    plt.title("Funcion de onda del oscilador armonico")
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
def orbital2D(nx = 1, ny=1, fx = 0.015, fy = 0.015, mx = cnts.m_e,
              my = cnts.m_e, d = [-0.15, 0.15, 0.005]):
    vf = np.vectorize(psi)
    x = np.arange(d[0], d[1], d[2])
    y = np.arange(d[0], d[1], d[2])
    X, Y = np.meshgrid(x, y)
    Z = vf(X, n=nx, f=fx, m=mx) * vf(Y, n=ny, f=fy, m=my)
    plt.contourf(x, y, Z, 50, cmap=cm.bone)
    plt.show()

# *****************************************
# Funcion para graficar la funcion de onda
# del oscilador armonico cuantico en varios
# niveles n = 3, 4, ... 11
# *****************************************
def orbitalN(n_max = 6, f = 0.015, m = cnts.m_e, d = [-0.2, 0.2, 0.0005]):
    vf = np.vectorize(psi)
    niveles = [0] * n_max
    for i in range(n_max):
        x = np.arange(d[0], d[1], d[2])
        niveles[i] = vf(x, n=i, f=f, m=m)**2
    niveles = np.array(niveles[::-1])
    plt.imshow(niveles, cmap=cm.bone)
    plt.axis("auto")
    plt.axis("off")
    plt.show()
