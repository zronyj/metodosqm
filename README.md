# Métodos para Mecánica Cuántica

A continuación se hallan varios métodos para realizar cálculos cuánticos sencillos.
La idea principal de estos es calcular la energía y densidad electrónica de los fenómenos:
* Partícula en una Caja
* Oscilador Armónico `oscilador_armonico.py`
* Átomo de Hidrógeno (H) `atomo_h.py`
* Molécula de Hidrógeno (H<sub>2</sub>) `hf-scf.py`

Cabe resaltar que lo que se halla en cada archivo solo son las funciones para realizar los cálculos; el cálculo como tal y los parámetros deseados deben ser ensamblados por el usuario. A continuación se proponen algunos ejemplos.

## Partícula en una caja

## Oscilador Armónico

## Átomo de Hidrógeno [H]
```python
# Funcion para el calculo de la funcion de onda radial (Rho)
rho = Rho(n=1, l=0, r=0.5, Z=1)

# Funcion para el calculo de la funcion de onda azimutal (Theta)
theta = Theta(x = 0.3, y = 0.5, m = 0)

# Funcion para el calculo de la funcion de onda polar (Phi)
phi = Phi(x = 0, y = 0, z = 1, l=1, m=-1)

# Calculo de la energia del atomo de hidrogeno
e_hidrogeno = E_h(n = 1, Z = 1)

# Grafica de la funcion de onda a lo largo del radio atomico
orbital_r(n = 1, l = 0, Z = 1, d = [-0.5, 5, 0.1])

# Grafica de los orbitales del atomo de hidrogeno
orbital2D(n = 1, l = 0, m = 0, Z = 1, d = [-4, 4, 40])
```

## Molécula de Hidrógeno [H<sub>2</sub>] (utilizando el método de Hartree-Fock)
```python
# Calculo de las constantes para la molecula de H2
orb = SCF(r = 1.4632, Z=[1,1], b1 = GTO["H"], b2 = GTO["H"], b = 2, vbs=True)

# Calculo de las energias de enlace (negativa) y antienlace (positiva)
E_anti_enlace = ener_orbs(x = orb["X"], f = orb["F"])

# Calculo de la energia electronica y total de la molecula
E_elec = ener_elec(p = orb["P"], h = orb["H"], f = orb["F"])
E_tot = ener_tot(r = 1.4632, Z=[1,1], elec=E_elec)

# Graficas de las funciones de onda a lo largo del eje
orbital(orb["C"], r = 1.4632, b1 = GTO["H"], b2 = GTO["H"], b = 2)
orbital2D(orb["C"], r = 1.4632, b1 = GTO["H"], b2 = GTO["H"], b = 2)
```