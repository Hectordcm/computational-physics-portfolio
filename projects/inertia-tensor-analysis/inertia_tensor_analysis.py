# -*- coding: utf-8 -*-
"""Compute center of mass and inertia tensors for 2D and 3D mass distributions."""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data_path = Path(__file__).with_name("datosmasas.csv")
if not data_path.exists():
    raise FileNotFoundError(f"Required data file not found: {data_path}")

df = pd.read_csv(data_path)

"""##Análisis 2D"""

import numpy as np
import matplotlib.pyplot as plt

#Datos 2D
masas=df['masas'].values
coords_2d=df[['x','y']].values

#Momento de orden cero
M_total=masas.sum()
print('El momento de orden cero es:', M_total)

#Momento de orden uno
x_cm=(coords_2d[:,0]*masas).sum()/M_total
y_cm=(coords_2d[:,1]*masas).sum()/M_total
coords_cm_2d=coords_2d-np.array([x_cm,y_cm])
print('El momento de orden uno es:', x_cm,y_cm)

#Momento de orden dos
I_xx=(masas*coords_cm_2d[:,1]**2).sum()
I_yy=(masas*coords_cm_2d[:,0]**2).sum()
I_xy=-(masas*coords_cm_2d[:,0]*coords_cm_2d[:,1]).sum()
I2D=np.array([[I_xx,I_xy], [I_xy,I_yy]])
print('El momento de orden dos es:', I2D)

#Autovalores y autovectores
autoval2D,autovec2D=np.linalg.eigh(I2D)
autovec2D=autovec2D.T

#Gráfico partículas y centro de masa
plt.figure(figsize=(7,7))
plt.scatter(df['x'],df['y'],s=masas*10,alpha=0.7,label='Partículas')
plt.scatter([x_cm],[y_cm],color='red',s=100,label='Centro de masa')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Distribución de masa en 2D')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

"""**¿Los vectores base del sistema cartesiano constituyen una base propia para esta distribución de masa? Esto es: ¿Los vectores cartesianos son autovectores del tensor momento de inercia?**

No, la matriz de autovectores del tensor momento de inercia es:
"""

print('Los autovalores son:', autoval2D)
print('Los autovectores son:', autovec2D)

"""Cuyos autovectores no corresponden a los vectores base del sistema cartesiano.

**Encuentre los ejes principales de inercia para esta distribución de masas. Esto es aquellos vectores propios del tensor de inercia, que forma una base ortogonal respecto a la cual la distribución de las masas se organiza de forma mas simple.**
"""

#Gráfico incluyendo ejes principales
plt.figure(figsize=(7,7))
plt.scatter(df['x'],df['y'],s=masas*10,alpha=0.7,label='Partículas')
plt.scatter([x_cm],[y_cm],color='red',s=100,label='Centro de masa')
for val,vec in zip(autoval2D, autovec2D.T):
    eje=vec*np.sqrt(val)*0.05
    plt.arrow(x_cm,y_cm,eje[0],eje[1],head_width=2,head_length=2,linewidth=2,label=f'Eje (λ={val:.1f})')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Distribución + Ejes Principales en 2D')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

"""**Encuentre la matriz de transformación de la base cartesiana a la base de autovectores conformada por los ejes principales.**

La matriz de transforamción es la matriz creada por los autvoectores del tensor de inercia.
"""

T = autovec2D
print(T)

"""## Análisis 3D

"""

from mpl_toolkits.mplot3d import Axes3D

#Datos 3D
coords_3d = df[['x', 'y', 'z']].values

#Momento de orden cero
print('El momento de orden cero es:', M_total)

#Momento de orden uno
z_cm=(coords_3d[:,2]*masas).sum()/M_total
coords_cm_3d=coords_3d-np.array([x_cm,y_cm,z_cm])
print('El momento de orden uno es:', x_cm,y_cm,z_cm)

#Momento de orden dos
I_xx=(masas*(coords_cm_3d[:,1]**2+coords_cm_3d[:,2]**2)).sum()
I_yy=(masas*(coords_cm_3d[:,0]**2+coords_cm_3d[:,2]**2)).sum()
I_zz=(masas*(coords_cm_3d[:,0]**2+coords_cm_3d[:,1]**2)).sum()
I_xy=-(masas*coords_cm_3d[:,0]*coords_cm_3d[:,1]).sum()
I_xz=-(masas*coords_cm_3d[:,0]*coords_cm_3d[:,2]).sum()
I_yz=-(masas*coords_cm_3d[:,1]*coords_cm_3d[:,2]).sum()
I3D=np.array([[I_xx, I_xy, I_xz],
                [I_xy, I_yy, I_yz],
                [I_xz, I_yz, I_zz]])
print('El momento de orden dos es:',I3D)

#Autovalores y autovectores
autoval3D, autovec3D = np.linalg.eigh(I3D)

#Gráfico partículas y centro de masa
fig=plt.figure(figsize=(10, 8))
ax=fig.add_subplot(111, projection='3d')
ax.scatter(df['x'],df['y'],df['z'],s=masas*10,alpha=0.6,label='Partículas')
ax.scatter([x_cm],[y_cm],[z_cm],color='red',s=100,label='Centro de masa')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Distribución de masa en 3D')
ax.legend()
plt.show()

"""**¿Los vectores base del sistema cartesiano constituyen una base propia para esta distribución de masa? Esto es: ¿Los vectores cartesianos son autovectores del tensor momento de inercia?**

No, la matriz de autovectores del tensor momento de inercia es:
"""

print('Los autovalores son:', autoval3D)
print('Los autovectores son:', autovec3D)

"""Cuyos autovectores no corresponden a los vectores base del sistema cartesiano.

**Encuentre los ejes principales de inercia para esta distribución de masas. Esto es aquellos vectores propios del tensor de inercia, que forma una base ortogonal respecto a la cual la distribución de las masas se organiza de forma mas simple.**
"""

#Gráfico incluyendo ejes principales
fig=plt.figure(figsize=(10, 8))
ax=fig.add_subplot(111, projection='3d')
ax.scatter(df['x'],df['y'],df['z'],s=masas*10,alpha=0.7,label='Partículas')
ax.scatter([x_cm],[y_cm],[z_cm],color='red',s=100,label='Centro de masa')
for val,vec in zip(autoval3D,autovec3D.T):
    escala=np.sqrt(val)*0.1
    eje=vec*escala
    ax.quiver(x_cm,y_cm,z_cm,eje[0],eje[1],eje[2],color='black',linewidth=2,arrow_length_ratio=0.02,label=f'Eje (λ={val:.1f})')
ax.set_box_aspect([1, 1, 1])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Distribución + Ejes Principales en 3D')
ax.legend()
plt.show()

"""**Encuentre la matriz de transformación de la base cartesiana a la base de autovectores conformada por los ejes principales.**

La matriz de transforamción es la matriz creada por los autvoectores del tensor de inercia.
"""

T=autovec3D
print(T)