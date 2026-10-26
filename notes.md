# MPC for Trajectory Tracking: Steering and Throttle/Brake

Notation follows the course slides *Low-Level Motion Planning* (Intelligent Vehicles and Autonomous Driving, Università di Trento). Slide numbers are given in brackets, e.g. [s.27].

---

## 1. Purpose

A high-level planner provides a **reference trajectory**, made of a **path** and a **speed profile** [s.4, s.28]. The low-level controller has to produce the actuator commands that make the vehicle track it:

- the **steering angle** $\delta$ at the front wheel (or the driver steering angle $\delta_D = \tau_D\,\delta$ [s.13]);
- a **single longitudinal command** $p \in [-1, 1]$, where $p > 0$ is throttle (motor torque) and $p < 0$ is brake.

Path tracking and speed tracking could be handled by two separate controllers (e.g. Stanley + PID [s.5]). That separation is only valid far from the vehicle's maximum performance [s.4]. Here we use a **combined** Model Predictive Controller (MPC) [s.27–29], which computes both commands together. It can therefore account for the coupling between lateral and longitudinal dynamics, for example slowing down to stay on the path in a tight corner, or braking before a curve while steering.

At each control step the MPC:

1. reads the current vehicle state and projects it onto the path (curvilinear coordinates);
2. predicts the vehicle motion over a horizon of length $\zeta_f - \zeta_i$ along the path [s.27];
3. solves an optimal control problem that keeps the predicted lateral position close to $n_{ref}$ and the speed close to $v_{x,ref}$, while respecting actuator and tire limits;
4. applies only the first control and repeats at the next step (*receding horizon*).

---

## 2. Notation

### 2.1 Path and curvilinear coordinates [s.7–9]

| Symbol | Meaning |
|---|---|
| $s$ | curvilinear abscissa along the reference path |
| $n$ | lateral position of the center of mass $G$ w.r.t. the path |
| $\xi$ | relative yaw angle of the vehicle w.r.t. the path tangent |
| $\theta$ | orientation of the local path tangent |
| $\psi$ | vehicle attitude (yaw), $\psi = \theta + \xi$ |
| $k(s)$ | path curvature, $k(s) = 1/R(s)$ |
| $\zeta$ | independent variable of the MPC (curvilinear abscissa) [s.27] |

### 2.2 Vehicle (single-track model) [s.12–17, s.39–41]

| Symbol | Meaning |
|---|---|
| $u$ | longitudinal speed of $G$ (vehicle frame) |
| $v$ | lateral speed of $G$ (vehicle frame) |
| $\Omega = \dot\psi$ | yaw rate |
| $\beta$ | chassis side-slip angle |
| $\delta$ | front-wheel steering angle; $\delta_D = \tau_D\delta$ at the steering wheel, $\tau_D \in [10, 20]$ |
| $L_f,\ L_r$ | distances of $G$ from the front and rear axle; $L = L_f + L_r$ (wheelbase) |
| $m,\ I_{zz}$ | vehicle mass and yaw moment of inertia |
| $\alpha_f,\ \alpha_r$ | front and rear tire side-slip angles |
| ${}^tF_{xf},\ {}^tF_{yf}$ | front tire longitudinal/lateral forces, in the **tire** frame |
| $F_{xr},\ F_{yr}$ | rear axle longitudinal/lateral forces |
| $C_f,\ C_r$ | front and rear cornering stiffness (the slides use a single $C_{F\alpha}$; see §3.3) |
| $F_{Ax}$ | aerodynamic drag force |
| $R_f,\ R_r$ | wheel radii |
| $K_{us}$ | understeering gradient [s.22–24] |
| $a_y$ | lateral acceleration |

### 2.3 Controls and references [s.27–28]

| Symbol | Meaning |
|---|---|
| $\vec{u}$ | control vector (note: **not** the longitudinal speed $u$) |
| $T_{mot},\ T_{br}$ | motor and brake torques |
| $p$ | single throttle/brake command, $p \in [-1,1]$ |
| $n_{ref}(\zeta)$ | reference lateral position (usually $0$: follow the path itself) |
| $v_{x,ref}(\zeta)$ | reference speed profile |
| $\mathcal{J},\ \mathcal{M},\ \mathcal{L}$ | cost functional, Mayer (terminal) term, Lagrange (running) term |

> **Notation clash.** The slides use $u$ both for the longitudinal speed and, as $\vec{u}$, for the control vector. Here the arrow always means the control vector. The Lagrange term is written $\mathcal{L}$ to distinguish it from the wheelbase $L$.

---

## 3. Vehicle Model

The vehicle is represented by a **dynamic single-track (bicycle) model with linear tires** [s.39–41]. The two wheels of each axle are lumped into one, which is acceptable since $\delta$ is small and the side-slip angles are nearly equal on the same axle [s.13–14].

### 3.1 Curvilinear kinematics [s.9]

$$
\dot s = \frac{u\cos\xi - v\sin\xi}{1 - n\,k(s)}, \qquad
\dot n = u\sin\xi + v\cos\xi, \qquad
\dot\xi = \Omega - k(s)\,\dot s
$$

### 3.2 Dynamics [s.40]

$$
\begin{aligned}
m(\dot u - \Omega v) &= F_{xr} + {}^tF_{xf}\cos\delta - {}^tF_{yf}\sin\delta - F_{Ax} \\
m(\dot v + \Omega u) &= F_{yr} + {}^tF_{yf}\cos\delta + {}^tF_{xf}\sin\delta \\
I_{zz}\,\dot\Omega &= L_f\left({}^tF_{yf}\cos\delta + {}^tF_{xf}\sin\delta\right) - L_r F_{yr}
\end{aligned}
$$

The term $-{}^tF_{yf}\sin\delta$ (omitted on s.40) is kept on purpose. It is the speed loss caused by cornering, and the controller should anticipate it. Drag is modeled as $F_{Ax} = \tfrac12 \rho\, C_d A\, u^2$, plus rolling resistance if needed. Self-aligning torques $M_{z}$ are neglected, as on s.41.

### 3.3 Tires [s.40–41]

Lateral forces are linear in the slip angles, which is valid for small slips:

$$
{}^tF_{yf} = C_f\,\alpha_f, \qquad F_{yr} = C_r\,\alpha_r
$$

$$
\alpha_f = \delta - \arctan\!\left(\frac{v + \Omega L_f}{u}\right), \qquad
\alpha_r = -\arctan\!\left(\frac{v - \Omega L_r}{u}\right)
$$

**Two stiffnesses instead of one.** The slides use a single $C_{F\alpha}$, but the ratio between $C_f$ and $C_r$ is what sets the steering characteristic. In the linear range:

$$
\delta - \frac{L}{R} = K_{us}\,a_y, \qquad K_{us} = \frac{m}{L}\left(\frac{L_r}{C_f} - \frac{L_f}{C_r}\right)
$$

with $K_{us} > 0$ for understeer and $K_{us} < 0$ for oversteer [s.22–24]. In practice, $C_f$ and $C_r$ can be identified by fitting $K_{us}$ from a steady-state (handling diagram) test.

### 3.4 Longitudinal actuation: from $p$ to forces

The slides propose $\{T_{mot}, T_{br}\}$ as controls [s.28]. Here they are merged into a single command $p$:

$$
T_{mot} = \max(p, 0)\; T_{mot}^{max}(u), \qquad
T_{br} = \max(-p, 0)\; T_{br}^{max}
$$

The torques are converted to forces at the ground through the wheel radius, neglecting wheel dynamics and longitudinal slip (s.16 is not used):

- the **motor** force $T_{mot}/R$ is applied to the driven axle;
- the **brake** force $T_{br}/R$ is split between the axles with a fixed ratio $\gamma$ (e.g. 0.6 front / 0.4 rear).

For a rear-wheel-drive car, for example:

$$
{}^tF_{xf} = -\gamma\,\frac{T_{br}}{R_f}, \qquad F_{xr} = \frac{T_{mot}}{R_r} - (1-\gamma)\frac{T_{br}}{R_r}
$$

The map has a kink at $p = 0$. For a gradient-based solver it can be smoothed, or the MPC can use the total longitudinal force $F_x$ as its input, with $p$ recovered afterwards from an identified inverse map $p = g^{-1}(F_x, u)$. A small deadband around $p = 0$ prevents chattering between throttle and brake.

### 3.5 State and input vectors

$$
\vec{x} = [\,n,\ \xi,\ u,\ v,\ \Omega,\ \delta\,]^T, \qquad
\vec{u} = [\,\dot\delta,\ p\,]^T
$$

Using the steering rate $\dot\delta$ as the input, with $\delta$ as a state, produces smooth commands and models the steering actuator's rate limit. Time $t$ can be added as a state when integrating over $\zeta$.

---

## 4. Optimal Control Problem [s.27–28]

### 4.1 Spatial formulation

The curvilinear abscissa $\zeta$ is the independent variable, and the horizon is a distance $\zeta_f - \zeta_i$ [s.27]. Every time derivative is converted using $\dot s$:

$$
\frac{d\vec{x}}{d\zeta} = \frac{\dot{\vec{x}}}{\dot s} = f(\vec{x}(\zeta), \vec{u}(\zeta)), \qquad \frac{dt}{d\zeta} = \frac{1}{\dot s}
$$

With this formulation the references $n_{ref}(\zeta)$ and $v_{x,ref}(\zeta)$ come directly from the planner, without any time parametrization. It requires $\dot s > 0$, i.e. forward driving.

### 4.2 Cost

$$
\min_{\vec{u}\in\mathcal{U}} \mathcal{J}, \qquad
\mathcal{J} = \mathcal{M}\big(\vec{x}(\zeta_i), \vec{x}(\zeta_f)\big) + \int_{\zeta_i}^{\zeta_f} \mathcal{L}\big(\vec{x}(\zeta), \vec{u}(\zeta)\big)\, d\zeta
$$

$$
\mathcal{L} = w_N\,(n - n_{ref})^2 + w_U\,(u - v_{x,ref})^2 + w_\xi\,\xi^2 + w_{\dot\delta}\,\dot\delta^2 + w_{p}\,\Delta p^2
$$

The first two terms are the ones on s.27. The others are regularizations: $w_\xi$ damps heading oscillations, and $w_{\dot\delta}$, $w_p$ smooth the commands.

### 4.3 Constraints $\vec{c}(\vec{x}, \vec{u}) \ge 0$ [s.28]

- **Steering:** $|\delta| \le \delta_{max}$ and $|\dot\delta| \le \dot\delta_{max}$
- **Throttle/brake:** $-1 \le p \le 1$
- **Tire validity:** $|\alpha_f|,\ |\alpha_r| \le \alpha_{max}$ (about 4–6°). Linear tires never saturate, so without this constraint the optimizer would assume unlimited lateral grip.
- **Grip (combined slip):** a friction ellipse $\left(\frac{a_x}{a_{x,max}}\right)^2 + \left(\frac{a_y}{a_{y,max}}\right)^2 \le 1$, a simple way to capture the lateral/longitudinal coupling recalled on s.4
- **Corridor (optional):** $n_{min}(\zeta) \le n \le n_{max}(\zeta)$ for road borders [s.7]

Constraints are better implemented as **soft** constraints (with slack variables) to avoid infeasibility.

### 4.4 Boundary conditions $\vec{b} = 0$ [s.28]

The initial state is fixed to the measured one, $\vec{x}(\zeta_i) = \vec{x}_0$. Terminal conditions are imposed softly through $\mathcal{M}$, e.g. $\mathcal{M} = w_{N,f}\,n(\zeta_f)^2 + w_{\xi,f}\,\xi(\zeta_f)^2$.

---

## 5. Implementation Notes

What `mpc.py` does about each practical issue:

- **Low speed.** The slip angles divide by $u$, and the spatial formulation needs $\dot s > 0$. The tire equations use $\max(u, 1\ \text{m/s})$ and the MPC keeps $u \ge 2$ m/s. Starting from standstill would need the kinematic single-track model [s.10–11] at low speed: here the lap starts already at the racing line speed.
- **Integration.** One RK4 step per node in the MPC; the simulated car uses 10 RK4 steps per control period. At low speed the $v$ and $\Omega$ dynamics are fast (stiff) and need more sub-steps or an implicit integrator.
- **Solver.** The full nonlinear problem is solved with CasADi + IPOPT at every control step (about 20 ms). Faster alternatives are the real-time iteration scheme (e.g. acados), or an LTV-MPC linearized along the previous solution, which yields a QP.
- **Horizon.** 2.5 s at the reference speed with 40 nodes, i.e. 30–150 m on this track. The horizon does not need to cover the whole braking distance, because the reference speed already slows down before the corners.
- **Delays.** Not modeled: the simulated car applies the command at once. On a real car, predict the initial state forward by the actuator and computation delay before solving, to avoid the instability mentioned on s.42.
- **Warm start.** The solver starts from the previous solution, moved to the new nodes. At the first step it starts from the curvature feedforward $\delta_{ff} = k\,(L + K_{us}u^2)$ [s.49].
- **Output.** $\dot\delta$ and $p$ of the first node are applied for one control period (0.05 s); $\delta$ (or $\delta_D = \tau_D\delta$) is what the steering actuator receives.

---

## 6. Model Choice and Validity

| Model | Valid range | Notes |
|---|---|---|
| Kinematic single track [s.10–11] | low speed, small $a_y$ (≲ 4 m/s²) | no tire parameters, well-behaved at low speed |
| **Dynamic single track, linear tires** (this document) | up to ≈ 0.5 g of lateral acceleration | needs $m$, $I_{zz}$, $L_f$, $L_r$, $C_f$, $C_r$; captures under/oversteer |
| Dynamic single track, Pacejka / combined slip | up to the grip limit | needed only near maximum performance [s.4, s.29] |

**Suggested workflow:** first get the full MPC pipeline (references, projection, $p$ mapping, constraints) working with the kinematic model. Then replace $f$ with the dynamic model described here. If the controller misbehaves after the swap, the cause is in the model rather than in the MPC setup.
