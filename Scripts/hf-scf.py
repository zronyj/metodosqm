# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 14:31:00 2017

@author: Rony J. Letona
@email: zronyj@gmail.com

Descripcion
-----------
Pequeno script para realizar el calculo de la energia y forma de los orbitales
de la molecula de hidrogeno H2 o HeH+ mediante el metodo de Hartree-Fock con
3 sets de bases diferentes: STO-2G, STO-3G y STO-6G.

Referencias
-----------

(1) Cruzeiro, V. W., Roitberg, A., & Polfer, N. C. (2016). Interactively
    Applying the Variational Method to the Dihydrogen Molecule: Exploring
    Bonding and Antibonding. Journal of Chemical Education, 93(9), 1578-1585.
    doi:10.1021/acs.jchemed.6b00017

(2) EMSL Basis Set Exchange. (n.d.). Recuperado noviembre 08, 2017,
    de https://bse.pnl.gov/bse/portal

(3) Feller, D. (1996). The role of databases in support of computational
    chemistry calculations. Journal of Computational Chemistry, 17(13),
    1571-1586.
    doi:10.1002/(sici)1096-987x(199610)17:13<1571::aid-jcc9>3.0.co;2-p

(4) Page, T. R., Boots, C. A., & Freitag, M. A. (2008). Quantum Chemistry:
    Restricted Hartree-Fock SCF Calculations Using Microsoft Excel. Journal of
    Chemical Education, 85(1), 159. doi:10.1021/ed085p159

(5) Schuchardt, K. L., Didier, B. T., Elsethagen, T., Sun, L., Gurumoorthi, V.,
    Chase, J., . . . Windus, T. L. (2007). Basis Set Exchange:  A Community
    Database for Computational Sciences. Journal of Chemical Information and
    Modeling, 47(3), 1045-1052. doi:10.1021/ci600510j

(6) Stewart, B., Hylton, D. J., & Ravi, N. (2013). A Systematic Approach for
    Understanding Slater–Gaussian Functions in Computational Chemistry. Journal
    of Chemical Education, 90(5), 609-612. doi:10.1021/ed300807y

(7) Szabo, A., & Ostlund, N. S. (1996). Modern quantum chemistry: introduction
    to advanced electronic structure theory. Mineola, NY: Dover Publications.

"""

from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# *************************
# Sets de bases
# Referencia (2)
# *************************
    
GTO = {"H":{}, "He":{}}

# Coeficiente exponencial Alfa para Hidrogeno
GTO["H"]["a"] = {2:[1.309756377, 0.233135974],
                 3:[3.42525091, 0.62391373, 0.16885540], 
                 6:[35.52322122, 6.513143725, 1.822142904,
                    0.625955266, 0.243076747, 0.100112428]}

# Coeficiente C para Hidrogeno
GTO["H"]["c"] = {2:[0.430128498, 0.678913531],
                 3:[0.15432897, 0.53532814, 0.44463454],
                 6:[0.00916359628, 0.04936149294, 0.16853830490,
                    0.37056279970, 0.41649152980, 0.13033408410]}

# Coeficiente exponencial Alfa para Helio
GTO["He"]["a"] = {2:[2.4328790, 0.4330510],
                  3:[6.36242139, 1.15892300, 0.31364979],
                  6:[65.98456824, 12.09819836, 3.384639924,
                     1.162715163, 0.451516322, 0.185959356]}

# Coeficiente C para Helio
GTO["He"]["c"] = {2:[0.430128498, 0.678913531],
                  3:[0.15432897, 0.53532814, 0.44463454],
                  6:[0.00916359628, 0.04936149294, 0.16853830490,
                     0.37056279970, 0.41649152980, 0.13033408410]}

# *********************************
# Funciones de onda (GTOs)
# Referencia (7) p. 153 eq. (3.203)
# *********************************
a, r, R = symbols("alpha r R")

phi = (2 * a / pi)**(3./4) * exp(-a * (r - R)**2)

def Phi (rn, Ru, base = GTO["H"], b = 3):
    wf = 0
    for i in range(b):
        wf += base["c"][b][i] * \
        phi.subs([(a, base["a"][b][i]), (r, rn), (R, Ru)])
    return wf

def Psi (r, R, base = GTO["H"], b = 3):
    wf = 0
    for i in range(b):
        a = base["a"][b][i]
        wf += base["c"][b][i] * (2 * a / np.pi)**0.75 * np.exp(-a * (r - R)**2)
    return wf

# ***********************************
# Funcion para graficar los orbitales
# con respecto a la distancia
# Referencia (1)
# Referencia (4)
# ***********************************
def orbital (c, r = 1.4632, b1 = GTO["H"], b2 = GTO["H"], b = 3):
    dom = np.arange(-3.5,r+3.55,0.05)
    psi1 = [c[0,0] * Psi(x, 0, base=b1, b=b) + c[1,0] * \
            Psi(x, r, base=b2, b=b) for x in dom]
    psi2 = [c[0,1] * Psi(x, 0, base=b1, b=b) + c[1,1] * \
            Psi(x, r, base=b2, b=b) for x in dom]
    psi3 = [y**2 for y in psi1]
    psi4 = [y**2 for y in psi2]
    plt.title("Funciones de onda resultantes del calculo HF-SCF")
    plt.plot(dom, psi1, "r--", label="$\Psi_1$")
    plt.plot(dom, psi2, "b--", label="$\Psi_2$")
    plt.plot(dom, psi3, "g-", label="$\Psi_{1}^{2}$")
    plt.plot(dom, psi4, "k-", label="$\Psi_{2}^{2}$")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.grid(True)
    plt.show()

def orbital2D (C, X, F, r = 1.4632, b1 = GTO["H"], b2 = GTO["H"], b = 3, delta=0.02):
    domx = np.arange(-3.5,r+3.5+delta,delta)
    domy = np.arange(-1.5-r/2,1.5+r/2+delta,delta)
    densmap1 = [[0]*len(domx)]*len(domy)
    densmap2 = [[0]*len(domx)]*len(domy)
    for y in range(len(domy)):
        q1 = Psi(domy[y], 0, base=b1, b=b)
        q2 = Psi(domy[y], 0, base=b2, b=b)
        for x in range(len(domx)):
            p1 = Psi(domx[x], 0, base=b1, b=b)
            p2 = Psi(domx[x], r, base=b2, b=b)
            densmap1[y][x] = (C[0,0]*p1*q1 + C[1,0]*p2*q2)**2
            densmap2[y][x] = (C[0,1]*p1*q1 + C[1,1]*p2*q2)**2
        densmap1[y] = tuple(densmap1[y])
        densmap2[y] = tuple(densmap2[y])
    dm1 = np.array(densmap1, dtype=np.float32)
    dm2 = np.array(densmap2, dtype=np.float32)
    energias = ener_orbs(X, F)
    fig = plt.figure()
    e1 = energias[0]
    a1 = fig.add_subplot(221)
    a1.title.set_text("Orbital 1\nE = " + str(e1))
    a1.imshow(dm1,cmap=cm.hot)
    e2 = energias[1]
    a2 = fig.add_subplot(224)
    a2.title.set_text("Orbital 2\nE = " + str(e2))
    a2.imshow(dm2,cmap=cm.hot)
    plt.show()

# **********************************
# Funcion para diagonalizar matrices
# Referencia (7) p. 144 eq. (3.169)
# **********************************
def diagon(m=np.matrix([[1,-0.5],[-2,1.5]])):
    m_eigenval, m_eigenvec = np.linalg.eig(m)
    mvp = np.matrix(np.diag([1/i**0.5 for i in m_eigenval.tolist()]),
                    dtype=np.float64)
    return (m_eigenvec * mvp)[::-1]

# **************************************
# Funcion para diagonalizar la matriz F'
# Referencia (7) p. 427
# **************************************
def diagon2(m=np.matrix([[1,-0.5],[-2,1.5]])):
    if abs(m[0,0] - m[1,1]) < 1E-8:
        temp = np.pi/4
    elif (abs(m[0,1]) < 1E-8) or (abs(m[1,0]) < 1E-8):
        temp = np.pi/4
    else:
        temp = 0.5 * np.arctan( 2*m[0,1] / (m[0,0] - m[1,1]) )
    ed = np.cos(temp)
    ec = np.sin(temp)
    return np.matrix( [[ed,ec],[ec,-ed]] )
    
# *******************************
# Matriz de traslape S
# Referencia (4)
# Referencia (7) p. 412 eq. (A.9)
# *******************************
def s (b1 = GTO["H"], b2 = GTO["H"], r = 0, b = 3):
    suv = 0
    for i in range(b-1,-1,-1):
        for j in range(b-1,-1,-1):
            du = (2*b1["a"][b][i]/np.pi)**0.75 * b1["c"][b][i]
            dv = (2*b2["a"][b][j]/np.pi)**0.75 * b2["c"][b][j]
            Sa = b1["a"][b][i] + b2["a"][b][j]
            suv += du * dv * (np.pi/Sa)**1.5 * \
            np.exp( -((b1["a"][b][i] * b2["a"][b][j])/Sa) * r**2 )
    return suv

def S(R = [0,1.4632], b1 = GTO["H"], b2 = GTO["H"], b = 3):
    MS = [[R[0]-R[0],R[1]-R[0]],
          [R[0]-R[1],R[1]-R[1]]]
    for i in range(2):
        for j in range(2):
            MS[i][j] = s(b1=b1, b2=b2, r=MS[i][j], b=b)
    return np.matrix(MS, dtype=np.float64)

# ********************************
# Matriz de energia cinetica T
# Referencia (7) p. 412 eq. (A.11)
# ********************************
def t(b1 = GTO["H"], b2 = GTO["H"], r = 0, b = 3):
    tuv = 0
    for i in range(b-1,-1,-1):
        for j in range(b-1,-1,-1):
            du = (2*b1["a"][b][i]/np.pi)**0.75 * b1["c"][b][i]
            dv = (2*b2["a"][b][j]/np.pi)**0.75 * b2["c"][b][j]
            Sa = b1["a"][b][i] + b2["a"][b][j]
            Ma = b1["a"][b][i] * b2["a"][b][j]
            tuv += du * dv * (Ma/Sa) * (3 - 2*Ma/Sa * r**2) * \
            (np.pi/Sa)**1.5 * np.exp( -((Ma/Sa) * r**2 ))
    return tuv

def T(R=[0,1.4632], b1 = GTO["H"], b2 = GTO["H"], b = 3):
    MT = [[R[0]-R[0],R[1]-R[0]],
          [R[0]-R[1],R[1]-R[1]]]
    for i in range(2):
        for j in range(2):
            MT[i][j] = t(b1=b1, b2=b2, r=MT[i][j], b=b)
    return np.matrix(MT, dtype=np.float64)

# ********************************
# Matriz de energia potencial V
# Referencia (7) p. 415 eq. (A.33)
# ********************************
def v (b1 = GTO["H"], b2 = GTO["H"], r=1.4632, Rp=[0,0,1,1], Z=[1,1], b = 3):
    vuv = 0
    for k in range(len(Z)):
        R = [l*r for l in Rp[2*k:2*k+2]]
        for i in range(b-1,-1,-1):
            for j in range(b-1,-1,-1):
                du = (2*b1["a"][b][i]/np.pi)**0.75 * b1["c"][b][i]
                dv = (2*b2["a"][b][j]/np.pi)**0.75 * b2["c"][b][j]
                Sa = b1["a"][b][i] + b2["a"][b][j]
                Rr = (b1["a"][b][i] * R[0] + b2["a"][b][j] * R[1]) / Sa
                Ruv = abs(R[1]-R[0])
                Rur = (Rr - Ruv)**2
                t = Sa*Rur
                if t < 1E-5:
                    F = 1 - t
                else:
                    F = 0.5*(np.pi/t)**0.5 * erf(t**0.5)
                temp = (-2*du*dv*np.pi*Z[k]/Sa) * F * \
                np.exp( -((b1["a"][b][i] * b2["a"][b][j])/Sa * Ruv**2 ))
                vuv += N(temp)
    return vuv

def V (r = 1.4632, Z=[1,1], b1 = GTO["H"], b2 = GTO["H"], b = 3):
    MV = [[[0,0,1,1],[0,1,1,0]],
          [[1,0,0,1],[1,1,0,0]]]
    for i in range(2):
        for j in range(2):
            MV[i][j] = v(b1=b1, b2=b2, r=r, Rp=MV[i][j], Z=Z, b=b)
    return np.matrix(MV, dtype=np.float64)

# *********************************
# Matriz de Hamiltonano H
# Referencia (7) p. 141 eq. (3.153)
# *********************************
def H (R=[0,1.4632], Z=[1,1], b1 = GTO["H"], b2 = GTO["H"], b = 3):
    th = T(R = R, b1 = b1, b2 = b2, b = 3)
    vh = V(r = R[1] - R[0], Z = Z, b1 = b1, b2 = b2, b = b)
    h = th + vh
    return h

# **********************************
# Integral de dos electrones <uv|ls>
# Referencia (4)
# Referencia (7) p. 416 eq. (A.41)
# **********************************
def ide (o = [0,0,0,0], R=[0,1], b1 = GTO["H"], b2 = GTO["H"], b = 3):
    res = 0
    base = [b1,b2]
    for i in range(b-1,-1,-1):
        for j in range(b-1,-1,-1):
            for k in range(b-1,-1,-1):
                for l in range(b-1,-1,-1):
                    s1 = base[o[0]]["a"][b][i] + base[o[1]]["a"][b][j]
                    s2 = base[o[2]]["a"][b][k] + base[o[3]]["a"][b][l]
                    s3 = s1 + s2
                    m1 = base[o[0]]["a"][b][i] * base[o[1]]["a"][b][j]
                    m2 = base[o[2]]["a"][b][k] * base[o[3]]["a"][b][l]
                    index = [i,j,k,l]
                    d = [(2*base[o[m]]["a"][b][index[m]]/np.pi)**0.75 * \
                         base[o[m]]["c"][b][index[m]] for m in range(4)]
                    contrac = np.prod(d)
                    R1 = R[o[1]] - R[o[0]]
                    R2 = R[o[3]] - R[o[2]]
                    Rr = (base[o[0]]["a"][b][i] * R[o[0]] + \
                          base[o[1]]["a"][b][j] * R[o[1]]) / s1
                    Rs = (base[o[2]]["a"][b][k] * R[o[2]] + \
                          base[o[3]]["a"][b][l] * R[o[3]]) / s2
                    Rsr = (Rs - Rr)**2
                    piterm = 2*np.pi**2.5 / (s1*s2*s3**0.5)
                    expterm = np.exp(-m1/s1 * R1**2 - m2/s2 * R2**2)
                    t = s1 * s2 / s3 * Rsr
                    if t < 1E-5:
                        F = 1 - t
                    else:
                        F = 0.5 * (np.pi/t)**0.5 * erf(t**0.5)
                    res += contrac * piterm * expterm * F
    return N(res)

# *************************************
# Matriz de dos electrones G
# Referencia (7) p. 141 eq. (3.154)
# Referencia (7) p. 427
# *************************************
def G (r=1.4632, p=np.matrix([[0,0],[0,0]], dtype=np.float64), b1=GTO["H"],
       b2=GTO["H"], b=3):
    g = [[0,0],[0,0]]
    for u in range(2):
        for v in range(2):
            for w in range(2):
                for x in range(2):
                    j = ide(o=[u,v,w,x], R=[0,r], b1=b1, b2=b2, b=b)
                    k = ide(o=[u,x,w,v], R=[0,r], b1=b1, b2=b2, b=b)
                    g[u][v] += p[w,x] * ( j - 0.5 * k )
    return np.matrix(g, dtype=np.float64)

# *********************************
# Matriz de densidad electronica P
# Referencia (7) p. 139 eq. (3.145)
# *********************************
def P (C=np.matrix([[0,0],[0,0]], dtype=np.float64)):
    p = [[C[0,0]*C[0,0], C[0,0]*C[1,0]],
         [C[1,0]*C[0,0], C[1,0]*C[1,0]]]
    return 2*np.matrix(p, dtype=np.float64)

# *************************************
# Metodo de Campo Autoconsistente (SCF)
# Referencia (7) p. 145
# *************************************
def SCF (r = 1.4632, Z=[1,1], b1 = GTO["H"], b2 = GTO["H"], b = 3, vbs=False):
    R = [0, r]
    if vbs: print("*) Generando matriz de traslape S.")
    s_scf = S(R=R, b1=b1, b2=b2, b=b)       # Referencia (7) p. 412 eq. (A.9)
    if vbs: print("\n*) Generando hamiltoniano H.")
    h_scf = H(R=R, Z=Z, b1=b1, b2=b2, b=b)  # Referencia (7) p. 141 eq. (3.153)
        
    # Diagonalizar matriz S y hallar matriz X
    if vbs: print("\n*) Diagonalizando matriz S y hallando matriz diagonal X.")
    X = diagon(m=s_scf)
    Xa = X.getH()
    
    # Estimar matriz de densidad P
    if vbs: print("\n*) Creando matriz de densidad P.")
    p_scf = np.matrix([[0,0],[0,0]], dtype=np.float64)  # Referencia (7) p. 148
    
    # Comenzar proceso iterativo
    if vbs: print("\n*) Comenzando con el SCF.")
    for iteracion in range(50):
        # Construir matriz de Fock F
        # F = H + G
        if vbs: print("\n**) Generando la matriz de Fock: calculando \
integrales de dos electrones.")
        g_scf = G(r=r, p=p_scf, b1=b1, b2=b2, b=b)
        f_scf = h_scf + g_scf     # Referencia (7) p. 141 eq. (3.154)
        
        # Construir matriz F'
        # F' = X_adj * F * X
        if vbs: print("**) Cambiando la base de F.")
        f_tra = Xa * f_scf * X
        
        # Diagonalizar matriz F y constuir matriz C'
        if vbs: print("**) Diagonalizando F' y generando C'.")
        c_tra = diagon2(m=f_tra)
        
        # Construir matriz C
        # C = X * C'
        if vbs: print("**) Construyendo matriz de coeficientes C.")
        c_scf = X * c_tra
        
        # Construir matriz P a partir de matriz C
        if vbs: print("**) Recalculando matriz de densidad P.")
        p_temp = P(C=c_scf)
        
        print("\nConcluida la " + str(iteracion + 1) + ". iteracion.\n")
        
        # Revisar convergencia
        if np.linalg.norm(p_temp - p_scf) < 1E-4: # Referencia (7) p. 148
            print("\n\n-->El campo autoconsistente SI ha convergido!")
            return {"S":s_scf,"H":h_scf,"X": X,"F":f_scf,"C":c_scf,"P":p_temp}
        else:
            p_scf = p_temp
    print("\n\n-->El campo autoconsistente NO ha convergido!\nRevisar supuestos.")
    return {"S":s_scf,"H":h_scf,"X": X,"F":f_scf,"C":c_scf,"P":p_temp}

# **********************************
# Calculo de energia de cada orbital
# Referencia (4)
# **********************************
def ener_orbs(x, f):
    return np.diag(x.getH() * f * x)

# *********************************
# Calculo de energia electronica
# Referencia (7) p. 150 eq. (3.184)
# *********************************
def ener_elec(p, h, f):
    m = 0.5 * p.reshape(1,4) * (h + f).reshape(1,4).T
    return m[0,0]

# *********************************
# Calculo de energia total
# Referencia (7) p. 150 eq. (3.185)
# *********************************
def ener_tot(r=1.4632, Z=[1,1], elec=0):
    nuc = Z[0] * Z[1] / r
    return nuc + elec