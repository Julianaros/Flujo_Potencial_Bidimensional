import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Necesario para gráficos 3D

# Función para cargar la función de corriente
def load_streamfunction(filename):
    data = np.loadtxt(filename)
    x = data[:, 0].astype(int)
    y = data[:, 1].astype(int)
    z = data[:, 2]

    Nx = x.max() + 1
    Ny = y.max() + 1

    U = np.zeros((Nx, Ny))
    for xi, yi, zi in zip(x, y, z):
        U[xi, yi] = zi

    return U.T  # Transponer para que coincidan los ejes (y, x)

# Función para cargar el campo de velocidades
def load_velocity_field(filename):
    x, y, u, v = [], [], [], []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() == "":
                continue
            xi, yi, ui, vi = map(float, line.strip().split())
            x.append(xi)
            y.append(yi)
            u.append(ui)
            v.append(vi)
    return np.array(x), np.array(y), np.array(u), np.array(v)

# Cargar función de corriente
U = load_streamfunction("streamfunction.dat")

# Gráfica 2D - contornos
plt.figure(figsize=(10, 5))
cp = plt.contour(U, levels=30, colors='k', linewidths=0.8)  # curvas negras, finas
plt.clabel(cp, inline=True, fontsize=9, fmt="%.1f")  # etiquetas suaves

plt.title("Líneas de corriente (ψ)", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel("y", fontsize=12)

# Opcional: mejorar estilo de ejes
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.4)

plt.tight_layout()
plt.savefig("streamfunction_2D_lines.png", dpi=300)
plt.show()


# Gráfica 3D - superficie
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
Xg, Yg = np.meshgrid(np.arange(U.shape[1]), np.arange(U.shape[0]))
ax.plot_surface(Xg, Yg, U, cmap='viridis')
ax.set_title("Función de corriente (streamfunction) - Superficie 3D")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("ψ (streamfunction)")
plt.tight_layout()
plt.savefig("streamfunction_3D.png")
plt.show()

# Cargar campo de velocidades
x, y, u, v = load_velocity_field("velocity_field.dat")

# Gráfica quiver (campo de velocidades)
plt.figure(figsize=(8, 5))
plt.quiver(x, y, u, v, scale=30, color='blue')
plt.title("Campo de velocidades (normalizado)")
plt.xlabel("x")
plt.ylabel("y")
plt.xlim(0, 30)   # Cambia estos valores según tu dominio
plt.ylim(0, 15)
#plt.axis('equal')
plt.tight_layout()
plt.savefig("campo_velocidades.png")
plt.show()
