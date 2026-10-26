"""Step 1: minimum lap time OCP, saved as the racing line for the MPC.

The car is a point mass driven by its tyre accelerations ax, ay, inside the friction circle and the engine
limits. Its state is written along the track centre line (abscissa s, see notes.md, section 4.1):
    n   lateral position w.r.t. the centre line (> 0 to the left)
    xi  heading relative to the centre line
    V   speed
    t   time
and the lap time t(L) is minimised over one lap that ends in the state it starts from.
"""

import os
import tomllib

import casadi as ca
import matplotlib.pyplot as plt
import numpy as np

import plots
from track import RACING_LINE, read_csv
from vehicle import G, RHO, load_vehicle, rk4


def min_lap_time(track, veh, ds, grip, margin):
    L = track["s"][-1]
    N = int(np.ceil(L / ds))
    s = np.linspace(0, L, N + 1)
    at = lambda key, z: np.interp(z, track["s"], track[key])
    kappa = at("kappa", (s[:-1] + s[1:]) / 2)  # curvature in the middle of each interval
    n_max = at("w_left", s) - veh["width"] / 2 - margin
    n_min = -at("w_right", s) + veh["width"] / 2 + margin

    # point mass dynamics w.r.t. s, integrated over one interval
    x, a, k = ca.SX.sym("x", 4), ca.SX.sym("a", 2), ca.SX.sym("k")

    def f(x):
        n, xi, V = x[0], x[1], x[2]
        sdot = V * ca.cos(xi) / (1 - n * k)
        dV = a[0] - 0.5 * RHO * veh["CdA"] * V**2 / veh["m"]
        return ca.vertcat(V * ca.sin(xi), a[1] / V - k * sdot, dV, 1) / sdot

    F = ca.Function("F", [x, a, k], [rk4(f, x, s[1] - s[0])])

    opti = ca.Opti()
    X = opti.variable(4, N + 1)  # [n, xi, V, t] at the nodes
    A = opti.variable(2, N)  # [ax, ay] in each interval
    n, xi, V, t = X[0, :], X[1, :], X[2, :], X[3, :]
    ax, ay = A[0, :], A[1, :]

    for i in range(N):
        opti.subject_to(X[:, i + 1] == F(X[:, i], A[:, i], kappa[i]))
    opti.subject_to(X[:3, N] == X[:3, 0])  # periodic lap
    opti.subject_to(t[0] == 0)

    m, R = veh["m"], veh["R"]
    opti.subject_to(ax**2 + ay**2 <= (grip * veh["mu"] * G) ** 2)  # friction circle
    opti.subject_to(opti.bounded(-veh["T_br"] / R, m * ax, veh["T_mot"] / R))  # brake and drive torque
    opti.subject_to(m * ax * V[:-1] <= veh["P_max"])  # engine power
    opti.subject_to(opti.bounded(n_min, n.T, n_max))  # track borders
    opti.subject_to(opti.bounded(-1, xi, 1))  # driving forward along the track
    opti.subject_to(V >= 5)  # away from V = 0, where the dynamics divide by V

    # lap time, plus a small penalty on the changes of acceleration (smooth controls)
    dA = A - ca.horzcat(A[:, -1], A[:, :-1])
    opti.minimize(t[N] + 1e-3 * ca.sumsqr(dA))

    # initial guess: centre line at constant speed
    v0 = 20.0
    opti.set_initial(V, v0)
    opti.set_initial(t, s / v0)
    opti.set_initial(ay, v0**2 * kappa)

    opti.solver("ipopt", {"expand": True}, {"max_iter": 1000, "sb": "yes"})
    print(f"Solving it with IPOPT: {N} intervals of {L / N:.2f} m\n")
    sol = opti.solve()
    return s, sol.value(X), sol.value(A)


def racing_line(track, s, X, A, ds=0.5):
    """Racing line sampled every ds metres along itself, in the format read by track.Path."""
    n, xi, V, t = X
    A = (A + np.roll(A, 1, axis=1)) / 2  # accelerations at the nodes, mean of the two intervals around
    ax, ay = np.hstack([A, A[:, :1]])
    at = lambda key: np.interp(s, track["s"], track[key])
    th = at("theta")
    x, y = at("x") - n * np.sin(th), at("y") + n * np.cos(th)
    d = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]  # distance along the line

    line = {
        "x": x,
        "y": y,
        "theta": th + xi,
        "kappa": ay / V**2,  # curvature of the point mass trajectory
        "v": V,
        "t": t,
        "ax": ax,
        "ay": ay,
        "n_left": at("w_left") - n,
        "n_right": -at("w_right") - n,
    }
    sl = np.linspace(0, d[-1], int(np.ceil(d[-1] / ds)) + 1)
    return {"s": sl} | {key: np.interp(sl, d, val) for key, val in line.items()}


if __name__ == "__main__":
    with open("configs/ocp.toml", "rb") as f:
        cfg = tomllib.load(f)
    veh = load_vehicle()
    track = read_csv(cfg["track"])
    print(f"Track {cfg['track']}: {track['s'][-1]:.1f} m")
    print("Building the minimum lap time OCP")

    s, X, A = min_lap_time(track, veh, cfg["ds"], cfg["grip"], cfg["margin"])
    line = racing_line(track, s, X, A)
    os.makedirs("output", exist_ok=True)
    np.savetxt(RACING_LINE, np.column_stack(list(line.values())), delimiter=",", header=",".join(line), comments="")
    print(f"\nRacing line: {line['s'][-1]:.1f} m, lap time {line['t'][-1]:.2f} s -> {RACING_LINE}")

    plots.racing_line(track, line)
    print("Close the figure window to exit")
    plt.show()
