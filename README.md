## SIMULACIÓN DEL FLUJO BIDIMENSIONAL ALREDEDOR DE UNA VIGA MEDIANTE EL MÉTODO DE VORTICIDAD-FUNCIÓN DE CORRIENTE

Este proyecto estudia el comportamiento de un **fluido viscoso e incompresible en régimen estacionario** fluyendo alrededor de una **viga sólida sumergida**. Se considera un flujo bidimensional en el plano (x, y), ignorando la dirección z por simetría.

###  Hipótesis Físicas

- **Incompresibilidad**: El fluido tiene densidad constante ρ. Entonces:
  ∇ · v = 0

- **Régimen estacionario**: El campo de velocidades no depende del tiempo:
  ∂/∂t = 0

- **Fluido viscoso**: Se incluye la viscosidad cinemática ν, que introduce disipación.

---

###  Ecuaciones Fundamentales

Partimos de las ecuaciones de Navier–Stokes para flujo incompresible:

- (v · ∇)v = -(1/ρ) ∇P + ν ∇²v  
- ∇ · v = 0

donde v = (vx, vy) y P es la presión.

---

###  Reformulación: Función de Corriente y Vorticidad

Para simplificar y garantizar incompresibilidad:

- Función de corriente u(x, y):
  - vx = ∂u/∂y
  - vy = -∂u/∂x

- Vorticidad escalar ω:
  - ω = ∂vy/∂x - ∂vx/∂y = -∇²u

- Ecuación para la evolución de ω:
  - ν ∇²ω = (∂u/∂y)(∂ω/∂x) - (∂u/∂x)(∂ω/∂y)

---

###  Interpretación Física

- Las líneas de nivel de u son las líneas de corriente.
- ω mide la rotación local del fluido. Si ω = 0, el flujo es irrotacional.
- Las ecuaciones de u y ω se resuelven acopladas para obtener el campo de flujo.

---

###  Número de Reynolds

Define el régimen del flujo:

- Re = V0 * L / ν

| Régimen del flujo        | Intervalo de Re        |
|--------------------------|------------------------|
| Flujo viscoso (Stokes)   | Re ≪ 1                 |
| Transición               | Re ≈ 1                 |
| Flujo inercial           | Re ≫ 1                 |
| Flujo adherido           | Re < 5                 |
| Formación de vórtices    | 5 < Re < 40            |
| Calle de von Kármán      | Re > 40                |

También se define el número de Reynolds de malla:

- R = V0 * h / ν

Controla la estabilidad numérica del esquema en diferencias finitas.
