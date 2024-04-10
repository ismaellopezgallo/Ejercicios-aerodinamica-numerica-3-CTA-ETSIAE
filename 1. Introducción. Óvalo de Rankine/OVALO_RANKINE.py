# %% 1. Importación de módulos
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %% 2. Definición de las funciones de cálculo del campo de velocidades
def velocity_field(x, z, Q, x_s, z_s):
    """
    Calcula el campo de velocidades inducido por una serie finita de manantiales y sumideros.
    
    Parámetros:
    x, z -- coordenadas en el plano xz
    Q -- magnitud de los manantiales/sumideros
    x_s, z_s -- coordenadas de los manantiales/sumideros
    
    Devuelve:
    u, w -- componentes del campo de velocidades
    """
    u = np.zeros_like(x, dtype='f8')
    w = np.zeros_like(z, dtype='f8')
    N = len(Q)
    for i in range(N):
        r = np.sqrt((x - x_s[i])**2 + (z - z_s[i])**2)
        u += Q[i] * (x - x_s[i]) / (2 * np.pi * r**2)
        w += Q[i] * (z - z_s[i]) / (2 * np.pi * r**2)
    return u, w



# %% 3. Cálculo del campo de velocidades, función potencial y función de corriente

# Definimos la magnitud y ángulo de la velocidad incidente
U_inf = 5
alpha = 0

# Definimos las coordenadas de los manantiales y sumideros
Q = np.array([50, -50])
x_s = np.array([-5, 5])
z_s = np.array([0, 0])

# Se definen también sus coordenadas en el plano complejo, para el cálculo del potencial.
t_s = x_s + 1j*z_s
num_sources = len(Q)

# Definimos una malla para las coordenadas x e z
nx, nz = 40, 40
x = np.linspace(-10, 10, nx)
z = np.linspace(-10, 10, nz)
Xm, Zm = np.meshgrid(x, z)
Tm = Xm + 1j*Zm

# Calculamos las componentes u y v del campo de velocidades inducido por los manantiales o sumideros
u, w = velocity_field(Xm, Zm, Q, x_s, z_s)

# Añadimos el componente de la velocidad incidente
u += U_inf * np.cos(alpha*np.pi/180)
w += U_inf * np.sin(alpha*np.pi/180)

# Construimos la funcion potencial
pot = 0
for s in range(num_sources):
    pot += Q[s]*np.log(Tm - t_s[s]) / (2*np.pi)
pot += U_inf * np.exp(-1j*alpha*np.pi/180) * Tm

# De la función potencial extraemos tanto el potencial de velocidades (Phi) como la función de corriente (Psi)
phi = pot.real
psi = pot.imag

# %% 4. Visualizamos la función de corriente
fig, ax = plt.subplots()                                              # Creamos la figura
divider = make_axes_locatable(ax)                                     # Ajustamos la posición de la barra de niveles
cax = divider.append_axes('right', size='5%', pad=0.05)

cont = ax.contourf(Xm, Zm, psi, 10, cmap='bwr')                       # Pintamos los contornos
ax.contour(Xm, Zm, psi, levels=[0], colors=['black'], linewidths=1.2) # Representamos la línea de corriente psi=0    
fig.colorbar(cont, cax=cax, orientation='vertical')                   # Pintamos la barra de niveles

N = np.sqrt(u**2 + w**2)                                                
ax.quiver(Xm, Zm, u/N, w/N, scale=50)                                           # Representamos el campo vectorial normalizado

ax.set_xlabel('x')                                                    # Etiquetamos de forma apropiada, y mostramos en pantalla
ax.set_ylabel('z')
ax.set_title(r'Función de corriente $\Psi$')
#ax.set(xlim=(-3, -1), ylim=(-1, 1))
ax.set_aspect('equal','box')
# plt.show()

# %% 5. Visualizamos el campo de velocidades
fig, axes = plt.subplots(1,2, figsize=(10,8)) 
ax1, ax2 = axes.flatten()                                             # Creamos la figura, con dos plots integrados

dividerU = make_axes_locatable(ax1)                                   # Ajustamos la posición de la barra de niveles en U
caxU = dividerU.append_axes('right', size='5%', pad=0.05)
dividerV = make_axes_locatable(ax2)                                   # Ajustamos la posición de la barra de niveles en W
caxV = dividerV.append_axes('right', size='5%', pad=0.05)

ax1.set_title(r'Velocidad horizontal $\Phi_x$')
cont = ax1.contourf(Xm, Zm, u, 10, cmap='hot')                        # Pintamos los contornos de U
ax1.contour(Xm, Zm, u, levels=10, colors=['black'], linewidths=0.8)   # Representamos lineas de isocontornos
fig.colorbar(cont, cax=caxU, orientation='vertical')                  # Pintamos la barra de niveles

ax2.set_title(r'Velocidad vertical $-\Phi_z$')
cont = ax2.contourf(Xm, Zm, w, 10, cmap='hot')                        # Pintamos los contornos
ax2.contour(Xm, Zm, w, levels=10, colors=['black'], linewidths=0.8)   # Representamos la línea de corriente psi=0    
fig.colorbar(cont, cax=caxV, orientation='vertical')                  # Pintamos la barra de niveles


for ax in axes:
    ax.set_xlabel('x')                                                # Etiquetamos de forma apropiada, y mostramos en pantalla
    ax.set_ylabel('z')
    #ax.set(xlim=(-3, -1), ylim=(-1, 1))
    ax.set_aspect('equal','box')
plt.tight_layout()
plt.show()

# %% Calculos del problema


# cargamos los datos de nuestro manantial

# si queremos introducirlo por teclado
# print('Coordenadas del punto')
# xx = float(input('x= '))
# zz = float(input('z= '))

# para el ejercicio lo vamos introducir desde el codigo para no repetirlo 
xx = -2.4325
zz = 3.4305

print(f'Coordenadas del punto: ({xx}, {zz})')

# dif_x = np.abs(x - xx)
# dif_z = np.abs(x - zz)
# x_arr, z_arr = np.argmin(dif_x), np.argmin(dif_z)       # con esto conseguimos las posiciones del array que contienen la informacion de nuestro punto
# print('x_arr, z_arr = ', Xm[x_arr, z_arr], Zm[x_arr, z_arr]) # para asegurarnos de que estamos en el punto que queremos

# resultado de la velocidad
u_res, w_res = velocity_field(xx, zz, Q, x_s, z_s)
u_res += U_inf      # le añadimos la velocidad incidente desde el infinito
print(f'Velocidad horizontal (u): {u_res} m/s')
print(f'Velocidad vertical (w): {w_res} m/s')

# calculo de la presion
p_ref = 1e5 
rho = 1
# lo calculamos con la ecuacion de bernouilli
p_est = p_ref + 0.5 * rho * (U_inf**2 - (u_res**2 + w_res**2))

print(f'Presión en el punto: {p_est} Pa')

# %% para pruebas 
# print('x=', x[x_arr], 'z=', x[z_arr])
# print(x[x_arr : x_arr + 2])
# print(z[z_arr :z_arr + 2])

# print(np.shape(u))

uu = 24.8732
ww = -30.5972
