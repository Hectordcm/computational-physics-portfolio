# -*- coding: utf-8 -*-
"""Ondas_gravitacionales_original.ipynb


    https://colab.research.google.com/drive/1hqTl1-OY8SZemfAlp2YX5VOnT1bc7bws

# Plantillas desde SXS
Primero instalamos la libreria sxs `!pip install sxs -q`
"""


# Commented out IPython magic to ensure Python compatibility.
import sxs # import sxs
import matplotlib.pyplot as plt
import numpy as np
# %config InlineBackend.figure_format = 'retina'

sxs_bhns_0001 = sxs.load("SXS:BHNS:0001")
w = sxs_bhns_0001.h

"""`w.metadata` los datos de los datos"""

print('tipo de objetos: '+w.metadata.object_types)

t0 = w.metadata.reference_time # tiempo de referencia, al inicio hay "junk radiation"
print(f't0 = {t0}')

w.data # todos los multipolos

ell, m = 2 , 2 # seleccionar el multipolo
i22 = w.index(ell, m) # seleccionar el indice del multipolo correspondiente
print(i22)
# guardar la onda 2 2 en w_2_2
w_2_2 = w[:, i22]

plt.figure()
plt.plot(w.t, np.real(w_2_2))
plt.plot(w.t, np.imag(w_2_2))
#plt.xlim(2000,2500)
plt.xlabel(r'$t$')
plt.ylabel(r'$h_{22}$')
plt.show()

"""Definimos las variables"""

phi = np.angle(w_2_2)
time = w_2_2.t

np.unwrap(-phi)[16460]

"""Grafica de los datos"""

#plt.plot(time, -phi)
plt.plot(time, np.unwrap(-phi))
plt.xlabel(r'$t/M$')
plt.xlim(0,2500)
plt.ylim(0,300)
plt.ylabel(r'$\varphi$')
plt.show()

"""Número de datos"""

len(time)

phi[16460]

"""Calcular los maximos y minimos y verificar si el tiempo esta equiespaciado



"""

indices_after_t0 = np.where(time >= t0)

unwrapped_phi_after_t0 = np.unwrap(-phi)[indices_after_t0]
time_after_t0 = time[indices_after_t0]

Min_after_t0 = np.min(unwrapped_phi_after_t0)
Max = np.max(np.unwrap(-phi))

print(f"Mínimo después de t={t0}: {Min_after_t0}")
print(f"Máximo general: {Max}")

índice_min_after_t0_filtered = np.where(unwrapped_phi_after_t0 == Min_after_t0)
índice_min_after_t0 = indices_after_t0[0][índice_min_after_t0_filtered[0]][0]

print(f"Índice del mínimo después de t={t0} (en array original): {índice_min_after_t0}")

índice_max = np.where(np.unwrap(-phi) == Max)[0][0]
print(f"Índice del máximo general (en array original): {índice_max}")

num_elements_between_min_max = índice_max - índice_min_after_t0 + 1
print(f"Número de elementos entre el índice {índice_min_after_t0} y {índice_max}: {num_elements_between_min_max}")

dt = np.diff(time)
print(f"Media: {dt.mean():.6f}, Min: {dt.min():.6f}, Max: {dt.max():.6f}")
print("Paso constante" if np.allclose(dt, dt.mean(), rtol=1e-3) else "Paso no constante")

"""Intervalo de interpolación"""

print(time[4586],np.unwrap(-phi)[4586])
print(time[19558],np.unwrap(-phi)[19558])

"""Encontrar los tiempos igualmente esapciados y los indices de los phi, y verifico que la pendiente es creciente




"""

Tiempos_lagrange = np.linspace(time[4586], time[19558], 20)
indices_lagrange = [np.abs(time - t).argmin() for t in Tiempos_lagrange]
Tiempos_lagrange = time[indices_lagrange]
selected_phi_lagrange = np.unwrap(-phi)[indices_lagrange]

print(Tiempos_lagrange)
print(selected_phi_lagrange)

Tiempos = np.linspace(time[4586], time[19558], 101)
indices = [np.abs(time - t).argmin() for t in Tiempos]
Tiempos = time[indices]
selected_phi = np.unwrap(-phi)[indices]

diff_selected_phi_lagrange = np.diff(selected_phi_lagrange)
print("Diferencia entre pasos sucesivos de selected_phi_lagrange:")
print(diff_selected_phi_lagrange)

diff_selected_phi = np.diff(selected_phi)
print("Diferencia entre pasos sucesivos de selected_phi:")
print(diff_selected_phi)

if np.any(diff_selected_phi < 0):
  print("Al menos una entrada en diff_selected_phi es negativa.")
else:
  print("Ninguna entrada en diff_selected_phi es negativa.")


plt.figure()
plt.plot(time, np.unwrap(-phi), label='Curva Original')
plt.plot(Tiempos_lagrange, selected_phi_lagrange, 'o', color='red', label='Puntos seleccionados')
plt.xlabel(r'$t/M$')
plt.ylabel(r'$\varphi$')
plt.xlim(0, 2500)
plt.ylim(0,300)
plt.legend()
plt.show()

"""Interpolaciones"""

from scipy.interpolate import interp1d, CubicSpline, PchipInterpolator

def lagrange_polynomial(x_data, y_data):
    n = len(x_data)
    def P(x):
        total = 0
        for i in range(n):
            Li = 1
            for j in range(n):
                if i != j:
                    Li *= (x - x_data[j]) / (x_data[i] - x_data[j])
            total += y_data[i] * Li
        return total
    return P

Lagrange = lagrange_polynomial(Tiempos_lagrange, selected_phi_lagrange)
cs_default = CubicSpline(Tiempos, selected_phi)
cs_natural = CubicSpline(Tiempos, selected_phi, bc_type='natural')
cs_1_0 = CubicSpline(Tiempos, selected_phi, bc_type=((1, 0.0), (1, 0.0)))
pchip = PchipInterpolator(Tiempos, selected_phi)

"""Grafica de las interpolaciones"""

t = np.linspace(time[4586], time[19558], 300)

plt.figure()
plt.plot(time[4586:19558], np.unwrap(-phi)[4586:19558], label='Curva original', alpha = 0.3)
#plt.plot(Tiempos_lagrange, selected_phi_lagrange, 'o', color='red', label='Puntos para la interpolacion de lagrange')
#plt.plot(t, Lagrange(t), label='Lagrange')
plt.plot(Tiempos[:-1], selected_phi[:-1], 'o', color='red', label='Puntos seleccionados', markersize=2, markeredgewidth=0.5)
plt.plot(t, cs_default(t), label='Cs', linestyle="--") # Not a knot
#plt.plot(t, cs_natural(t), label='Cs Natural') # Segunda derivada = 0 en los extremos
#plt.plot(t, cs_1_0(t), label='Cs 1_0') # Primera derivada = 0 en los extremos
#plt.plot(t, pchip(t), label='Pchip')
plt.xlabel(r'$t/M$')
plt.ylabel(r'$\varphi$')
plt.ylim(0,200)
plt.xlim(0,2500)
plt.legend()
plt.title('Resultado de la interpolación')
plt.show()

plt.figure()
plt.plot(time[4586:19558], np.unwrap(-phi)[4586:19558], label='Curva original', alpha = 0.3)
#plt.plot(Tiempos, selected_phi, 'o', color='red', label='Selected Points')
#plt.plot(t, Lagrange(t), label='Lagrange')
plt.plot(Tiempos[:-1], selected_phi[:-1], 'o', color='red', label='Puntos seleccionados', markersize=5, markeredgewidth=0.5)
plt.plot(t, cs_default(t), label='Cubic Spline', linestyle="--")
#plt.plot(t, cs_natural(t), label='Cubic Spline Natural')
#plt.plot(t, cs_1_0(t), label='Cs 1_0')
#plt.plot(t, pchip(t), label='Pchip')
plt.xlabel(r'$t/M$')
plt.ylabel(r'$\varphi$')
plt.ylim(130,190)
plt.xlim(2050,2360)
plt.legend()
plt.title('Resultado de la interpolación')
plt.show()

"""Crear grupo de entrenamiento y grupo de prueba"""

x_train_lagrange = Tiempos_lagrange[0::2]
x_test_lagrange = Tiempos_lagrange[1::2]
y_train_lagrange = selected_phi_lagrange[0::2]
y_test_lagrange_0 = selected_phi_lagrange[1::2]

x_train = Tiempos[0::2]
x_test = Tiempos[1::2]
y_train = selected_phi[0::2]
y_test = selected_phi[1::2]

Lagrange_interp = lagrange_polynomial(x_train, y_train)
Cs_interp = CubicSpline(x_train, y_train)
Cs_natural_interp = CubicSpline(x_train, y_train, bc_type='natural')
Cs_1_0_train = CubicSpline(x_train, y_train, bc_type=((1, 0.0), (1, 0.0)))
Pchip_interp = PchipInterpolator(x_train, y_train)

y_test_lagrange_1 = Lagrange_interp(x_test_lagrange)
y_test_Cs = Cs_interp(x_test)
y_test_Cs_natural = Cs_natural_interp(x_test)
y_test_Cs_1_0 = Cs_1_0_train(x_test)
y_test_Pchip = Pchip_interp(x_test)

SE_lagrange = (y_test_lagrange_0 - y_test_lagrange_1)**2
SE_Cs = (y_test - y_test_Cs)**2
SE_Cs_natural = (y_test - y_test_Cs_natural)**2
SE_Cs_1_0 = (y_test - y_test_Cs_1_0)**2
SE_Pchip = (y_test - y_test_Pchip)**2

AE_lagrange = np.abs(y_test_lagrange_0 - y_test_lagrange_1)
AE_Cs = np.abs(y_test - y_test_Cs)
AE_Cs_natural = np.abs(y_test - y_test_Cs_natural)
AE_Cs_1_0 = np.abs(y_test - y_test_Cs_1_0)
AE_Pchip = np.abs(y_test - y_test_Pchip)

MSE_lagrange = np.mean(SE_lagrange)
MSE_CS = np.mean(SE_Cs)
MSE_Cs_natural = np.mean(SE_Cs_natural)
MSE_Cs_1_0 = np.mean(SE_Cs_1_0)
MSE_Pchip = np.mean(SE_Pchip)

MAE_lagrange = np.mean(AE_lagrange)
MAE_CS = np.mean(AE_Cs)
MAE_Cs_natural = np.mean(AE_Cs_natural)
MAE_Cs_1_0 = np.mean(AE_Cs_1_0)
MAE_Pchip = np.mean(AE_Pchip)

print(f"MSE lagrange = {MSE_lagrange:.3f}")
print(f"MSE CS = {MSE_CS:.3f}")
print(f"MSE Cs natural = {MSE_Cs_natural:.3f}")
print(f"MSE Cs 1_0 = {MSE_Cs_1_0:.3f}")
print(f"MSE Pchip = {MSE_Pchip:.3f}")

print(f"MAE lagrange = {MAE_lagrange:.3f}")
print(f"MAE CS = {MAE_CS:.3f}")
print(f"MAE Cs natural = {MAE_Cs_natural:.3f}")
print(f"MAE Cs 1_0 = {MAE_Cs_1_0:.3f}")
print(f"MAE Pchip = {MAE_Pchip:.3f}")

plt.figure()
plt.scatter(x_test_lagrange, SE_lagrange, label='Lagrange')
plt.scatter(x_test, SE_Cs, label='CS')
plt.scatter(x_test, SE_Cs_natural, label='CS natural')
plt.scatter(x_test, SE_Cs_1_0, label='CS 1_0')
plt.scatter(x_test, SE_Pchip, label='Pchip')
plt.yscale('log')
plt.xlabel('x')
plt.ylabel('SE')
plt.legend(loc='upper center')
plt.show()

plt.figure()
plt.scatter(x_test_lagrange, AE_lagrange, label='Lagrange')
plt.scatter(x_test, AE_Cs, label='CS')
plt.scatter(x_test, AE_Cs_natural, label='CS natural')
plt.scatter(x_test, AE_Cs_1_0, label='CS 1_0')
plt.scatter(x_test, AE_Pchip, label='Pchip')
plt.yscale('log')
plt.xlabel('x')
plt.ylabel('AE')
plt.legend(loc='upper center')
plt.show()

"""Derivada de la fase"""

Tiempo = np.linspace(time[4586], time[19558], 500)
ind = [np.abs(time - t).argmin() for t in Tiempo]
Tiempo = time[ind]
selected_phi = np.unwrap(-phi)[ind]
print(f"Tiempo inicial: {Tiempo[0]}")
print(f"Tiempo final: {Tiempo[-1]}")
h = np.diff(Tiempo)[0]
print(f"paso es de {h}")

v_forward = (selected_phi[1:] - selected_phi[:-1]) / h
v_central = (selected_phi[2:] - selected_phi[:-2]) / (2 * h)
v_backward = (selected_phi[1:] - selected_phi[:-1]) / h

cs = CubicSpline(Tiempo, selected_phi)
Velocidad = cs.derivative(1)(Tiempo)
D_WM = cs.derivative(2)(Tiempo)
w_isco = 2*6**-1.5

"""Error"""

Vel  = Velocidad[:-1]
promedio = (v_forward + v_backward) / 2

err_rel  = np.mean((promedio - Vel) / Vel)
err_abs  = np.mean(np.abs(promedio - Vel))
err_quad = np.mean((promedio - Vel)**2)

err_rel1  = np.mean((v_forward - v_backward) / v_backward)
err_abs1  = np.mean(np.abs(v_forward - v_backward))
err_quad1 = np.mean((v_forward - v_backward)**2)

#La comparación con el central tiene problemas de dimensionalidad.

print(err_abs, err_quad, err_rel, err_rel1, err_abs1, err_abs1)
#print(f'La frecuencia angular es{Vel/2*np.pi}')

import matplotlib.pyplot as plt
Tiempo = np.linspace(time[4586], time[19558], 500)
ind = [np.abs(time - t).argmin() for t in Tiempo]
Tiempo = time[ind]
selected_phi = np.unwrap(-phi)[ind]
# --- Gráfica phi ---
plt.figure(figsize=(15, 6))
plt.plot(Tiempo, selected_phi, 'ko-', label='Discretizado', markersize=3)
plt.ylabel("phi")
plt.xlabel("t/M")
plt.legend()
plt.title("Evolución de phi")
plt.show()


# --- Gráfica velocidades ---
plt.figure()
#plt.plot(Tiempo[:-1], v_forward, '.--', label='Forward (1 Orden)', markersize=3)
#plt.plot(Tiempo[1:], v_backward, '.--', label='Backward (1 Orden)', markersize=3)
plt.plot(Tiempo, Velocidad, '-', label='Curva original')
#plt.plot(Tiempo[1:-1], v_central, linestyle='--', label='Central (2 Orden)', markersize=3)
plt.axhline(y=w_isco, label=f"w_isco = {w_isco}", color = 'green')
#plt.xlim(2250,2350)
plt.ylabel("WM")
plt.xlabel("t/M")
plt.legend()
plt.title("Resultados derivación numérica")
plt.show()

"""Ceros con Biseccion y Newton-Raphson"""

f = dict(zip(Tiempo, Velocidad - w_isco))
df = dict(zip(Tiempo, D_WM))

def biseccion(f, x0, x1, kmax=200, tol=1.e-3):
    x0  = Tiempo[np.abs(x0 - Tiempo).argmin()]
    x1  = Tiempo[np.abs(x1 - Tiempo).argmin()]
    f0 = f[x0]

    for k in range(1,kmax):
        x2 = (x0 + x1) / 2
        x2  = Tiempo[np.abs(x2 - Tiempo).argmin()]
        f2 = f[x2]
        if np.abs(f2) < tol:
          x2new = x2
          break
        if f0 * f2 < 0:
            x1 = x2
        else:
            x0, f0 = x2, f2
        x2new = (x0 + x1) / 2
        xdiff = abs(x2new - x2)
        if abs(xdiff / x2new) < tol:
            break
        else:
            x2new = None
    print(f"{k} iteraciones, {x2new}")
    return x2new

def newton_raphson(f, df, x0, tol=1e-3, max_iter=100):
    x0  = Tiempo[np.abs(x0 - Tiempo).argmin()]
    x = x0
    for i in range(max_iter):
        fx = f[x]
        if abs(fx) < tol:
            return x, i+1
            break
        dfx = df[x]
        if abs(dfx) < 1e-12:
            print ("Derivada pequeña")
        x = x - fx/ dfx
        x = Tiempo[np.abs(x - Tiempo).argmin()]
    print(f"{i} iteraciones, {x}")

ini_guess = 2140
fin_guess = 2180
midd_guess = 2160

t_wisco_bi = biseccion(f, ini_guess, fin_guess)
t_wisco_ne_tuple = newton_raphson(f, df, midd_guess)
t_wisco_ne = t_wisco_ne_tuple[0]
iterations = t_wisco_ne_tuple[1]


print(f'Bisección (Orden 1): El tiempo en el que WM = W_ISCO es {t_wisco_bi} Con una diferencia de {f[t_wisco_bi]}')
print(f'Newton (Orden 2): El tiempo en el que WM = W_ISCO es {t_wisco_ne} Con una diferencia de {f[t_wisco_ne]} en {iterations} iteraciones')

index_t_wisco = np.where(time == t_wisco_bi)[0][0]
print(f"El índice de t_wisco es: {index_t_wisco}")



"""Edos

"""

from scipy.integrate import solve_ivp

m1 = sxs_bhns_0001.metadata.get("reference_mass1", None)
m2 = sxs_bhns_0001.metadata.get("reference_mass2", None)

M = m1 + m2
m1 = m1 / M  # normalizadas a M=1
m2 = m2 / M
M = 1.0
nu = (m1 * m2) / (m1 + m2)**2

t_junk = time[4586]
t_isco  = t_wisco_bi
F_junk = Velocidad[0] / (2 * np.pi)
v0 = (np.pi * M * F_junk)**(1/3)

def dvdt(t, v):
    gamma = 0.577216
    C1 = 1
    C2 = -(743/336 + 11/4 * nu) * v**2
    C3 = 4*np.pi * v**3
    C4 = (34103/18144 + 13661/2016*nu + 59/18*nu**2) * v**4
    C5 = -(4159/672*np.pi + 189/8*np.pi*nu) * v**5
    C6 = (16447322263/139708800
             - 1712/105*np.pi**2
             + (451/48*np.pi**2 - 56198689/217728)*nu
             + 541/896*nu**2
             - 5605/2592*nu**3
             - 856/105*np.log(16*v**2)) * v**6
    C7 = -(4415/4032 - 358675/6048*nu + 91/45*nu**2) * np.pi * v**7

    dv = (32/5) * nu / M * v**9 * (C1 + C2 + C3 + C4 + C5 + C6 + C7)
    return dv

# 1) RK4 manual

def rk4(f, y0, t_array):
    n = len(t_array)
    y = np.zeros(n)
    y[0] = y0
    h = t_array[1] - t_array[0]

    for i in range(n - 1):
        k1 = f(t_array[i], y[i])
        k2 = f(t_array[i] + 0.5*h, y[i] + 0.5*h*k1)
        k3 = f(t_array[i] + 0.5*h, y[i] + 0.5*h*k2)
        k4 = f(t_array[i] + h, y[i] + h*k3)
        y[i+1] = y[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    return y

# 1) RK45 y LSODA

sol_rk45 = solve_ivp(dvdt, t_span=[t_junk, t_isco], y0=[v0], method='RK45', dense_output = True, rtol=1e-5, atol=1e-7)
sol_lsoda = solve_ivp(dvdt,t_span=[t_junk, t_isco], y0=[v0], method='LSODA', dense_output = True, rtol=1e-5, atol=1e-7)

t_eval = np.linspace(t_junk, t_isco, 2000)
v_rk45_dense  = sol_rk45.sol(t_eval)[0]
v_lsoda_dense = sol_lsoda.sol(t_eval)[0]


dt = 0.1
t_rk4 = np.arange(t_junk, t_isco + dt, dt)
v_rk4 = rk4(dvdt, v0, t_rk4)


plt.figure(figsize=(12,6))
plt.plot(t_rk4, v_rk4, label="RK4 manual (dt=0.1)", color="blue", linewidth=2, alpha=0.6)
plt.plot(t_eval, v_rk45_dense, label="RK45", color="orange", linestyle="--", linewidth=2)
plt.plot(t_eval, v_lsoda_dense, label="LSODA", color="green", linestyle=":", linewidth=2)
plt.xlabel("t")
plt.ylabel("v(t)")
plt.title("Comparación de resultados - TaylorT4")
plt.legend()
plt.grid(True)

plt.show()

i1, i2 = 4586, 16461 #t_junk hasta t_isco
t_sxs = time[i1:i2]
phi_sxs = np.unwrap(-phi)[i1:i2]
cs = CubicSpline(t_sxs, phi_sxs)
omega_sxs = cs.derivative(1)(t_sxs)
f_gw_sxs = omega_sxs / (2*np.pi)
v_sxs = (np.pi * M * f_gw_sxs)**(1/3)

interp_v_sxs = interp1d(t_sxs, v_sxs, kind="cubic", fill_value="extrapolate")
v_sxs_interp = interp_v_sxs(t_rk4)


plt.figure(figsize=(12,6))
plt.plot(t_rk4, v_rk4, label="RK4 clásico (dt = 0.1)", linewidth=2)
plt.plot(t_eval, v_rk45_dense, label="RK45 (dense output)", linestyle="--", linewidth=2)
plt.plot(t_eval, v_lsoda_dense,label="LSODA (dense output)", linestyle=":", linewidth=2)
plt.plot(t_rk4, v_sxs_interp,'--', label="SXS interpolado", linewidth=2)
plt.xlabel("t")
plt.ylabel("v(t)")
plt.title("Comparación: RK4 – RK45 – LSODA – SXS")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#Zoom

plt.figure(figsize=(12,6))
plt.plot(t_rk4, v_rk4, label="RK4 clásico", linewidth=2)
plt.plot(t_eval, v_rk45_dense, label="RK45 (dense)", linestyle="--", linewidth=2)
plt.plot(t_eval, v_lsoda_dense,label="LSODA (dense)", linestyle=":", linewidth=2)
plt.plot(t_rk4, v_sxs_interp, '--', label="SXS interpolado", linewidth=2)
plt.xlim(t_rk4[0], t_rk4[8000])
plt.ylim(0.29,0.32)
plt.xlabel("t")
plt.ylabel("v(t)")
plt.title("Zoom al comienzo: RK4 – RK45 – LSODA – SXS")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

sol_ref_rk45 = solve_ivp(dvdt, [t_junk, t_isco], [v0], method="RK45", rtol=1e-14, atol=1e-16, dense_output=True)
v_ref_rk45_final = sol_ref_rk45.sol(t_isco)[0]
t_dense = np.linspace(t_junk, t_isco, 2000)
v_dense = sol_ref_rk45.sol(t_dense)[0]

N_array = np.array([10, 20, 50, 100, 200, 500])
h_array = (t_isco - t_junk) / N_array
error_rk4 = []
for N in N_array:
    t_array = np.linspace(t_junk, t_isco, N + 1)
    v_rk4 = rk4(dvdt, v0, t_array)
    err = abs(v_rk4[-1] - v_ref_rk45_final)
    error_rk4.append(err)
error_rk4 = np.array(error_rk4)

plt.figure(figsize=(9,6))
plt.loglog(h_array, error_rk4, 's-', label='Error RK4')
plt.loglog(h_array, 1e-18 * h_array**4, 'k--', label='Pendiente $O(h^4)$')
plt.title("Convergencia RK4 usando RK45 como referencia")
plt.xlabel("Paso h")
plt.ylabel("Error")
plt.grid(True)
plt.legend()
plt.show()
