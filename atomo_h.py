# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 10:20:15 2017

@author: Rony J. Letona
@email: zronyj@gmail.com

Descripcion
-----------
Pequeno script para realizar el calculo de la energia y forma de los orbitales
del atomo de hidrogeno partiendo de sus funciones.

"""

import numpy as np
import scipy.special as spe
import scipy.constants as cnts
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Cambio de coordenadas de cartesianas a esfericas
rho = lambda x, y, z: (x**2 + y**2 + z**2)**0.5
theta = lambda x, y: np.arctan(y/x)
phi = lambda x, y, z: np.arctan((x**2 + y**2)**0.5 / z)

# Funcion para calcular factorial n!
factorial = lambda n: np.prod( np.array( [i for i in range(1,n)] ) )

# *******************************
# Funcion de coordenadas radiales
# *******************************
def Rho(n, l, r, Z=1):
    rho = 2 * Z * r / n
    el = rho**l
    N_ln = (Z**3 * factorial(n - l - 1) / (n**4 * factorial(l + n)))**0.5
    L_nl = spe.assoc_laguerre(rho, l, n) # Polinomio asociado de Laguerre
    return 2 * el * N_ln * np.exp(-rho/2) * L_nl

# *********************************
# Funcion de coordenadas azimutales
# *********************************
def Theta(x, y, m):
    t = theta(x, y)
    return np.exp(m*t*1j) / (2*np.pi)**0.5

# ******************************
# Funcion de coordenadas polares
# ******************************
def Phi(x, y, z, l, m):
    N_lm = ( (2*l+1)/2 * factorial(l - m)/factorial(l + m))**0.5
    p = phi(x, y, z)
    P_lm = spe.lpmv(m, l, np.cos(p)) # Polinomio asociado de Legendre
    return N_lm * P_lm

# *******************
# Armonicos Esfericos
# *******************
def Y_ae(x, y, z, l, m):
    ytp = Theta(x, y, m) * Phi(x, y, z, l, m)
    #ytp = spe.sph_harm(m, l, phi(x, y, z), theta(x, y))
    return np.real(ytp)

# **************************************
# Funcion de onda del atomo de hidrogeno
# **************************************
def Psi2(x, y, z, n, l, m, Z=1):
    r = rho(x, y, z)
    return (Rho(n, l, r, Z) * Y_ae(x, y, z, l, m))**2

# *************************************
# Energia del atomo de hidrogeno para n
# *************************************
def E_h(n, Z=1):
    a = -Z**2 * cnts.m_e**4 * np.exp(4)
    d = ( 8 * cnts.h**2 * cnts.epsilon_0**2 * n**2)
    return a / d
                                    

# ***************************************
# Funcion para graficar el orbital radial
# en 2 dimensiones
# ***************************************
def orbital_r(n, l, Z=1, d=[-1,5,0.1]):
    x = np.arange(d[0], d[1], d[2])
    vr = np.vectorize(Rho)
    y = vr(n, l, x, Z)
    y2 = vr(n, l, x, Z)**2
    plt.title("Funcion de onda del atomo de hidrogeno $H \cdot$")
    plt.plot(x, y, "r--", label="$\Psi$")
    plt.plot(x, y2, "b-", label="$\Psi^2$")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axes = plt.gca()
    axes.set_ylim([-0.2, 1])
    axes.set_xlim([-0.5, d[1] - 0.5])
    plt.grid(True)
    plt.show()

# *****************************************
# Funcion para graficar los orbitales del
# atomo de hidrogeno en 3 diferentes planos
# *****************************************
def orbital2D(n=1, l=0, m=0, Z=1, d=[-4,4,40]):
    x = np.linspace(d[0], d[1], d[2])
    y = np.linspace(d[0], d[1], d[2])
    z = np.linspace(d[0], d[1], d[2])
    Xi, Yi, Zi = np.meshgrid(x, y, z)
    vf = np.vectorize(Psi2)
    orb = vf(x=Xi, y=Yi, z=Zi, n=n, l=l, m=m, Z=Z)
    plano_xz = orb[:,:,int(d[2]/2)]
    plano_yz = orb[:,int(d[2]/2),:]
    plano_xy = orb[int(d[2]/2),:,:]
    fig = plt.figure()
    
    ax = fig.add_subplot(221)
    ax.title.set_text("Eje X - plano yz")
    plt.contourf(y, z, plano_yz, 20, cmap=cm.bone)
    plt.colorbar()
    
    ay = fig.add_subplot(222)
    ay.title.set_text("Eje Y - plano xz")
    plt.contourf(x, z, plano_xz, 20, cmap=cm.bone)
    plt.colorbar()
    
    az = fig.add_subplot(223)
    az.title.set_text("Eje Z - plano xy")
    plt.contourf(x, y, plano_xy, 20, cmap=cm.bone)
    plt.colorbar()
    
    fig.tight_layout()
    fig.set_size_inches(w=5.4,h=4.4)
    