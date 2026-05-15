import numpy as np
import matplotlib.pyplot as plt

# 1.- Gráficas de Equilibración (Ejercicio 1)

# Cargamos el archivo saltando el encabezado (skiprows=1)
datos_eq = np.loadtxt('datos_eq.txt', skiprows=1)

iteracion = datos_eq[:, 0]
mag_eq = datos_eq[:, 1]
ene_eq = datos_eq[:, 2]

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig1.suptitle('Ejercicio 1: Tiempo de Equilibración (T=2.0)', fontsize=14)

ax1.plot(iteracion, mag_eq, color='blue')
ax1.set_xlabel('Iteraciones')
ax1.set_ylabel('Magnetización <M>')
ax1.grid(True)

ax2.plot(iteracion, ene_eq, color='red')
ax2.set_xlabel('Iteraciones')
ax2.set_ylabel('Energía <E>')
ax2.grid(True)


# 2.- Gráficas de Transición de Fase (Ejercicio 2)

datos_fase = np.loadtxt('datos_fase.txt', skiprows=1)

T = datos_fase[:, 0]
M = datos_fase[:, 1]
E = datos_fase[:, 2]
Cv = datos_fase[:, 3]
Chi = datos_fase[:, 4]

fig2, axs = plt.subplots(2, 2, figsize=(12, 8))
fig2.suptitle('Ejercicio 2: Transición de Fase', fontsize=14)

axs[0, 0].plot(T, M, marker='o', color='blue')
axs[0, 0].set_ylabel('Magnetización Media |<M>|')
axs[0, 0].grid(True)

axs[0, 1].plot(T, E, marker='o', color='red')
axs[0, 1].set_ylabel('Energía Media <E>')
axs[0, 1].grid(True)

axs[1, 0].plot(T, Cv, marker='o', color='green')
axs[1, 0].set_xlabel('Temperatura T')
axs[1, 0].set_ylabel('Calor Específico Cv')
axs[1, 0].grid(True)

axs[1, 1].plot(T, Chi, marker='o', color='purple')
axs[1, 1].set_xlabel('Temperatura T')
axs[1, 1].set_ylabel('Susceptibilidad Chi')
axs[1, 1].grid(True)


# 3.- Gráficas de Correlación (Ejercicio 3)

datos_corr = np.loadtxt('datos_corr.txt', skiprows=1)

R = datos_corr[:, 0]
C_Tc = datos_corr[:, 1]
C_Hot = datos_corr[:, 2]

fig3, ax3 = plt.subplots(figsize=(8, 5))
fig3.suptitle('Ejercicio 3: Longitud de Correlación', fontsize=14)

ax3.plot(R, C_Tc, marker='o', label='T = 2.3 (Punto Crítico)', color='orange')
ax3.plot(R, C_Hot, marker='s', label='T = 3.0 (Fase Desordenada)', color='teal')
ax3.set_xlabel('Distancia R')
ax3.set_ylabel('Correlación C(R)')
ax3.legend()
ax3.grid(True)

# Mostrar todas las figuras
plt.show()