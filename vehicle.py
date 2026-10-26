"""Vehicle parameters and the dynamic single-track model with linear tyres (see notes.md, section 3).

State and input names used everywhere:
    u, v, r   longitudinal speed, lateral speed, yaw rate (vehicle frame)
    delta     front wheel steering angle
    ddelta    steering rate (input)
    p         throttle / brake command in [-1, 1], > 0 throttle, < 0 brake (input)
"""

import tomllib

import casadi as ca

G = 9.81  # m/s^2
RHO = 1.2  # kg/m^3, air density


def load_vehicle(file="configs/vehicle.toml"):
    with open(file, "rb") as f:
        return tomllib.load(f)


def understeer_gradient(veh):
    L = veh["lf"] + veh["lr"]
    return veh["m"] / L * (veh["lr"] / veh["Cf"] - veh["lf"] / veh["Cr"])


def pos(x, eps=1e-2):
    """Smooth max(x, 0), so that the solver sees no kink at p = 0."""
    return 0.5 * (x + ca.sqrt(x**2 + eps**2))


def body(veh):
    """[u, v, r, delta], [ddelta, p] -> time derivatives, and [alpha_f, alpha_r, ax, ay]."""
    z, w = ca.SX.sym("z", 4), ca.SX.sym("w", 2)
    u, v, r, d = z[0], z[1], z[2], z[3]
    ddelta, p = w[0], w[1]
    m, R, lf, lr = veh["m"], veh["R"], veh["lf"], veh["lr"]

    # slip angles, with the speed kept away from 0
    uc = ca.fmax(u, 1.0)
    alpha_f = d - ca.atan((v + lf * r) / uc)
    alpha_r = -ca.atan((v - lr * r) / uc)

    # p -> motor torque (rear axle, limited by power) and brake torque (split with gamma)
    T_max = veh["T_mot"] - pos(veh["T_mot"] - veh["P_max"] * R / uc, 10.0)
    T_mot, T_br = pos(p) * T_max, pos(-p) * veh["T_br"]
    Fxf = -veh["gamma"] * T_br / R
    Fxr = (T_mot - (1 - veh["gamma"]) * T_br) / R
    Fyf, Fyr = veh["Cf"] * alpha_f, veh["Cr"] * alpha_r
    F_drag = 0.5 * RHO * veh["CdA"] * u**2

    # forces and yaw moment in the vehicle frame
    Fx = Fxr + Fxf * ca.cos(d) - Fyf * ca.sin(d)
    Fy = Fyr + Fyf * ca.cos(d) + Fxf * ca.sin(d)
    Mz = lf * (Fyf * ca.cos(d) + Fxf * ca.sin(d)) - lr * Fyr

    dz = ca.vertcat(r * v + (Fx - F_drag) / m, -r * u + Fy / m, Mz / veh["Izz"], ddelta)
    return ca.Function("body", [z, w], [dz, ca.vertcat(alpha_f, alpha_r, Fx / m, Fy / m)])


def spatial(veh):
    """[n, xi, u, v, r, delta] derivatives w.r.t. the path abscissa s; k is the path curvature.
    n = lateral offset from the path, xi = heading relative to the path."""
    x, w, k = ca.SX.sym("x", 6), ca.SX.sym("w", 2), ca.SX.sym("k")
    n, xi, u, v, r = x[0], x[1], x[2], x[3], x[4]
    dz, aux = body(veh)(x[2:], w)
    sdot = (u * ca.cos(xi) - v * ca.sin(xi)) / (1 - n * k)
    dx = ca.vertcat(u * ca.sin(xi) + v * ca.cos(xi), r - k * sdot, dz) / sdot
    return ca.Function("spatial", [x, w, k], [dx, aux])


def cartesian(veh):
    """[X, Y, psi, u, v, r, delta] derivatives w.r.t. time: the simulated car."""
    y, w = ca.SX.sym("y", 7), ca.SX.sym("w", 2)
    psi, u, v = y[2], y[3], y[4]
    dz, aux = body(veh)(y[3:], w)
    dy = ca.vertcat(u * ca.cos(psi) - v * ca.sin(psi), u * ca.sin(psi) + v * ca.cos(psi), y[5], dz)
    return ca.Function("cartesian", [y, w], [dy, aux])


def rk4(f, x, h, steps=1):
    """Integrate dx = f(x) over h with fixed-step Runge-Kutta 4."""
    h = h / steps
    for _ in range(steps):
        k1 = f(x)
        k2 = f(x + h / 2 * k1)
        k3 = f(x + h / 2 * k2)
        k4 = f(x + h * k3)
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x
