##  SIMULACIÓN DEL FLUJO BIDIMENSIONAL ALREDEDOR DE UNA VIGA MEDIANTE EL MÉTODO DE VORTICIDAD-FUNCIÓN DE CORRIENTE

Este proyecto estudia el comportamiento de un **fluido viscoso e incompresible en régimen estacionario** fluyendo alrededor de una **viga sólida sumergida**. Se considera un flujo bidimensional en el plano (x, y), ignorando la dirección z por simetría.

###  Hipótesis Físicas

- **Incompresibilidad**: El fluido tiene densidad constante \( \rho \), lo que implica:
  \[
  \nabla \cdot \vec{v} = 0
  \]

- **Régimen estacionario**: El campo de velocidades no depende del tiempo:
  \[
  \frac{\partial}{\partial t} = 0
  \]

- **Fluido viscoso**: Se considera la viscosidad cinemática \( \nu \), responsable de la disipación en el flujo.

---

###  Ecuaciones Fundamentales

Partimos de las ecuaciones de Navier–Stokes para flujo incompresible:

\[
(\vec{v} \cdot \nabla)\vec{v} = -\frac{1}{\rho} \nabla P + \nu \nabla^2 \vec{v}
\]

junto con la condición de incompresibilidad:

\[
\nabla \cdot \vec{v} = 0
\]

donde \( \vec{v} = (v_x, v_y) \) es el campo de velocidades y \( P \) la presión.

---

###  Reformulación: Función de Corriente y Vorticidad

Para garantizar automáticamente la incompresibilidad y simplificar el sistema:

- Se define la **función de corriente** \( u(x, y) \):
  \[
  v_x = \frac{\partial u}{\partial y}, \quad v_y = -\frac{\partial u}{\partial x}
  \]

- La **vorticidad escalar** \( \omega \) se define como:
  \[
  \omega = \nabla \times \vec{v} = \frac{\partial v_y}{\partial x} - \frac{\partial v_x}{\partial y}
  \]

- Lo anterior lleva a la ecuación de Poisson para la función de corriente:
  \[
  \omega = -\nabla^2 u
  \]

- Además, al tomar el rotacional de Navier–Stokes se obtiene una ecuación para la evolución de \( \omega \):
  \[
  \nu \nabla^2 \omega = \frac{\partial u}{\partial y} \frac{\partial \omega}{\partial x} - \frac{\partial u}{\partial x} \frac{\partial \omega}{\partial y}
  \]

---

###  Interpretación Física

- Las **líneas de nivel** de \( u \) representan las **líneas de corriente**, es decir, trayectorias tangentes al flujo.
- La vorticidad \( \omega \) mide la **rotación local del fluido**.
- Si \( \omega = 0 \), el flujo es **irrotacional**.
- Las ecuaciones anteriores se resuelven de forma acoplada para obtener el patrón de flujo alrededor de la viga.

---

###  Número de Reynolds

El número de Reynolds caracteriza el régimen del flujo:

\[
Re = \frac{V_0 L}{\nu}
\]

| Régimen del flujo        | Intervalo de Re        |
|--------------------------|------------------------|
| Flujo viscoso (Stokes)   | \( Re \ll 1 \)         |
| Transición               | \( Re \sim 1 \)        |
| Flujo inercial           | \( Re \gg 1 \)         |
| Flujo adherido           | \( Re < 5 \)           |
| Formación de vórtices    | \( 5 < Re < 40 \)      |
| Calle de von Kármán      | \( Re > 40 \)          |

También se define el **número de Reynolds de malla** como:

\[
R = \frac{V_0 h}{\nu}
\]

Este controla la estabilidad numérica del esquema y debe mantenerse dentro de límites adecuados.

---

