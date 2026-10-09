# -*- coding: utf-8 -*-
"""Monte Carlo simulation of the 2D Ising model using Metropolis updates.

Tracks total energy and magnetization, estimates a burn-in point, and plots configurations.
"""

import numpy as np
import matplotlib.pyplot as plt

J = 1.0
N2 = N1 = 30
pasos = 100000
T = 1.0
semilla = 69
guardar_cada = pasos // 6

def energia(espines, J=1.0):
    N1, N2 = espines.shape
    E = 0.0
    for i in range(N1):
        for j in range(N2):
            der1 = espines[i, (j + 1) % N2]
            der2 = espines[(i + 1) % N1, j]
            E += -J * espines[i, j] * (der1 + der2)
    return E

def deltaE(espines, i, j, J=1.0):
    N1, N2 = espines.shape
    sij = espines[i, j]
    interacción = (espines[(i - 1) % N1, j] + espines[(i + 1) % N1, j] + espines[i, (j - 1) % N2] + espines[i, (j + 1) % N2])
    return 2.0 * J * sij * interacción #se consigue de Edesp -Eant con 2 int

def magnetizacion_por_spin(espines):
    return espines.sum() / espines.size

def ejecutar_metropolis(N1=N1, N2= N2, J=1.0, T=1.0, pasos=100000, semilla=semilla, guardar_cada=None, pasos_snap=None, init_aleatorio=True):

    rng = np.random.default_rng(semilla)
    espines = rng.choice([1, -1], size=(N1, N2))
    E = energia(espines, J)
    suma_M = espines.sum()
    energias = [E]
    magnetizaciones = [suma_M / (N1 * N2)]

    captura_spin = [] #método para tomar los slices de la cadena a graficar
    pasos_snap_set = set(pasos_snap) if pasos_snap is not None else set()

    for paso in range(1, pasos + 1):
        i = rng.integers(0, N1)
        j = rng.integers(0, N2)

        s_ant = espines[i, j]

        dE = deltaE(espines, i, j, J)

        if dE < 0.0 or rng.random() < np.exp(-dE / T):#metrópolis cond.

            espines[i, j] = -s_ant
            E += dE #optimización
            suma_M += -2 * s_ant #optm.

        energias.append(E)
        magnetizaciones.append(suma_M / (N1 * N2))

        ###
        if (guardar_cada is not None and paso % guardar_cada == 0) or (paso in pasos_snap_set):
            captura_spin.append((paso, espines.copy()))
        ###


    return np.array(energias), np.array(magnetizaciones), captura_spin, espines.copy()

def burn_inn(energias, tol=1e-3, window=200, consecutivos=5):
    n = len(energias)
    if n <= window or window < 2:
        return 0

    k = np.ones(window) / window
    j = np.convolve(energias, k, mode='valid')#Método pro de suavizado para revisar estabilidad
    diffs = np.abs(np.diff(j))


    L = len(diffs)
    for p in range(0, L - consecutivos + 1):#revisón de estabilidad
        if np.all(diffs[p:p + consecutivos] < tol):
            burn_in = p + window
            return min(burn_in, n - 1)
    return 0


######## opcional
def graficar_toro(espines, R=2.0, r=0.8, elev=20, azim=30, cmap='coolwarm'):
    filas, columnas = espines.shape
    u = np.linspace(0, 2 * np.pi, filas, endpoint=False)
    v = np.linspace(0, 2 * np.pi, columnas, endpoint=False)
    U, V = np.meshgrid(u, v, indexing='ij')
    X = (R + r * np.cos(V)) * np.cos(U)
    Y = (R + r * np.cos(V)) * np.sin(U)
    Z = r * np.sin(V)
    min_val = espines.min()
    max_val = espines.max()
    if max_val != min_val:
        norm = (espines - min_val) / (max_val - min_val)
    else:
        norm = np.full_like(espines, 0.5, dtype=float)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, facecolors=plt.get_cmap(cmap)(norm), rcount=filas, ccount=columnas,
                    shade=False, linewidth=0, antialiased=False)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    plt.show()

pasos_snap = None
Toro = True

energias, magnetizaciones, instantaneas, espines_finales = ejecutar_metropolis(
    N1=N1, N2=N2, J=J, T=T, pasos=pasos, semilla=semilla,
    guardar_cada=guardar_cada, pasos_snap=pasos_snap, init_aleatorio=True
)

burn_in = burn_inn(energias, tol=1e-2, window=400, consecutivos=3)

if burn_in + 1 < len(magnetizaciones):
    avg_M = np.mean(magnetizaciones[burn_in + 1:])
    std_err = np.std(magnetizaciones[burn_in + 1:], ddof=1) / np.sqrt(len(magnetizaciones[burn_in + 1:]))
else:
    avg_M = np.mean(magnetizaciones)
    std_err = np.std(magnetizaciones, ddof=1) / np.sqrt(len(magnetizaciones))


plt.figure(figsize=(10, 4))
plt.plot(energias, linewidth=0.6)
plt.axvline(burn_in, color='C3', linestyle='--', label=f'burn_in={burn_in}')
plt.xlabel("Paso")
plt.ylabel("Energía total")
plt.title(f"Energía vs paso (T={T}, N={N1}x{N2})")
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(magnetizaciones, linewidth=0.6)
plt.axvline(burn_in, color='C3', linestyle='--', label=f'burn_in={burn_in}')
plt.xlabel("Paso")
plt.ylabel("Magnetización por espín")
plt.title(f"Magnetización vs paso (T={T})")
plt.legend()
plt.tight_layout()
plt.show()


for paso, conf in instantaneas:
    plt.figure(figsize=(4, 4))
    plt.imshow(conf, cmap='binary', interpolation='nearest', aspect='equal')
    plt.axis('off')
    plt.title(f'Paso {paso}')
    plt.show()

plt.figure(figsize=(5, 5))
plt.imshow(espines_finales, cmap='binary', interpolation='nearest', aspect='equal')
plt.axis('off')
plt.title('Configuración final')
plt.show()

graficar_toro(espines_finales, R=2.5, r=0.9)
print("burn_in:", burn_in)
print("M_medio:", avg_M)
print("Error:", std_err)