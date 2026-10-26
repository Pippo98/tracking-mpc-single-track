"""Step 2: the MPC drives one lap of the racing line computed in step 1.

The simulated car is the MPC's own single-track model, integrated in time (vehicle.cartesian), starting on the
line at its speed. At every control period: project the car on the line, solve the MPC, apply the first input.
"""

import os
import sys
import tomllib

import casadi as ca
import matplotlib.pyplot as plt
import numpy as np

import plots
from mpc import MPC
from track import RACING_LINE, Path, read_csv
from vehicle import cartesian, load_vehicle, rk4

LOG = "output/mpc_log.csv"


def car_model(veh, dt, steps=10):
    """[X, Y, psi, u, v, r, delta], [ddelta, p] -> state after dt."""
    y, w = ca.SX.sym("y", 7), ca.SX.sym("w", 2)
    f = cartesian(veh)
    return ca.Function("car", [y, w], [rk4(lambda y: f(y, w)[0], y, dt, steps)])


def simulate(path, mpc, car, veh, dt):
    wrap = lambda a: np.mod(a + np.pi, 2 * np.pi) - np.pi
    y = np.array([path("x", 0), path("y", 0), path("theta", 0), path("v", 0), 0, 0, 0])
    s = t = p = 0.0
    log = []
    while s < path.L:
        s, n = path.project(y[0], y[1], s)
        xi = wrap(y[2] - path("theta", s))
        w, info = mpc.solve(path, np.array([n, xi, *y[3:]]), s, p)
        log.append([t, s, n, xi, *y, *w, info["cpu"], info["ok"]])
        if abs(n) > 10 or abs(xi) > 1:
            print(f"\nlost control at s = {s:.0f} m", end="")
            break
        y = np.array(car(y, w)).ravel()
        y[6] = np.clip(y[6], -veh["delta_max"], veh["delta_max"])
        p, t = w[1], t + dt
        print(f"\rs = {s:.0f} / {path.L:.0f} m", end="")
    print()
    keys = "t s n xi X Y psi u v r delta ddelta p cpu ok".split()
    return dict(zip(keys, np.array(log).T))


if __name__ == "__main__":
    with open("configs/mpc.toml", "rb") as f:
        cfg = tomllib.load(f)
    with open("configs/ocp.toml", "rb") as f:
        track = read_csv(tomllib.load(f)["track"])
    veh = load_vehicle()
    if not os.path.exists(RACING_LINE):
        sys.exit(f"{RACING_LINE} not found: run step 1 first (uv run solve_ocp.py)")
    path = Path()
    print(f"Racing line {RACING_LINE}: {path.L:.1f} m, lap time {path.data['t'][-1]:.2f} s")

    print(f"Building the MPC: {cfg['N']} nodes over {cfg['T']} s")
    mpc = MPC(veh, cfg["N"], cfg["T"], cfg["weights"])
    print(f"Driving one lap, solving the MPC every {cfg['dt']} s")
    log = simulate(path, mpc, car_model(veh, cfg["dt"]), veh, cfg["dt"])
    np.savetxt(LOG, np.column_stack(list(log.values())), delimiter=",", header=",".join(log), comments="")

    s, n = log["s"], log["n"]
    print()
    if s[-1] >= path.L:
        print(f"lap time {np.interp(path.L, s, log['t']):.2f} s (racing line {path.data['t'][-1]:.2f} s)")
    print(f"lateral error: max {np.abs(n).max():.3f} m, rms {np.sqrt(np.mean(n**2)):.3f} m")
    print(f"speed error rms {np.sqrt(np.mean((log['u'] - path('v', s)) ** 2)):.3f} m/s")
    cpu = log["cpu"] * 1e3
    print(f"solve time: mean {cpu.mean():.1f} ms, max {cpu.max():.1f} ms, failed {int((log['ok'] == 0).sum())}")
    print(f"\nLog -> {LOG}")

    plots.mpc_lap(track, path, log)
    print("Close the figure window to exit")
    plt.show()
