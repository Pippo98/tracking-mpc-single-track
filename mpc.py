"""Nonlinear MPC tracking the racing line (see notes.md, section 4), solved with CasADi + IPOPT.

State  x = [n, xi, u, v, r, delta]  lateral offset and heading w.r.t. the racing line, then as in vehicle.py
Input  w = [ddelta, p]              steering rate, throttle / brake

The independent variable is the abscissa s along the racing line: node i of the horizon is at z[i],
and the nodes are spaced by T / N seconds at the reference speed.
"""

import time

import casadi as ca
import numpy as np

from vehicle import G, body, rk4, spatial, understeer_gradient


class MPC:
    def __init__(self, veh, N, T, weights):
        self.veh, self.N, self.dt = veh, N, T / N
        W = weights
        f, aux = spatial(veh), body(veh)

        # one step of the model over a distance h, with the path curvature k
        x, w, k, h = ca.SX.sym("x", 6), ca.SX.sym("w", 2), ca.SX.sym("k"), ca.SX.sym("h")
        F = ca.Function("F", [x, w, k, h], [rk4(lambda y: f(y, w, k)[0], x, h)])

        opti = ca.Opti()
        X, U = opti.variable(6, N + 1), opti.variable(2, N)
        S = opti.variable(4, N)  # slacks of the soft constraints
        x0, p0 = opti.parameter(6), opti.parameter()  # current state, last throttle / brake command
        ref = opti.parameter(5, N)  # rows: node spacing, curvature, reference speed, right and left border
        ds = ref[0, :]

        a_max = veh["mu"] * G
        for i in range(N):
            opti.subject_to(X[:, i + 1] == F(X[:, i], U[:, i], ref[1, i], ds[i]))
            alpha_f, alpha_r, ax, ay = ca.vertsplit(aux(X[2:, i], U[:, i])[1])
            opti.subject_to((alpha_f / veh["alpha_max"]) ** 2 <= 1 + S[0, i])  # linear tyre range
            opti.subject_to((alpha_r / veh["alpha_max"]) ** 2 <= 1 + S[1, i])
            opti.subject_to((ax / a_max) ** 2 + (ay / a_max) ** 2 <= 1 + S[2, i])  # friction circle
            opti.subject_to(X[0, i + 1] - ref[4, i] <= S[3, i])  # track borders
            opti.subject_to(ref[3, i] - X[0, i + 1] <= S[3, i])

        opti.subject_to(X[:, 0] == x0)
        opti.subject_to(opti.bounded(-veh["delta_max"], X[5, 1:], veh["delta_max"]))
        opti.subject_to(opti.bounded(-veh["ddelta_max"], U[0, :], veh["ddelta_max"]))
        opti.subject_to(opti.bounded(-1, U[1, :], 1))
        opti.subject_to(X[2, 1:] >= 2.0)  # moving forward: the model is written per metre of path
        opti.subject_to(ca.vec(S) >= 0)

        dp = U[1, :] - ca.horzcat(p0, U[1, :-1])
        cost = (
            W["n"] * X[0, 1:] ** 2
            + W["u"] * (X[2, 1:] - ref[2, :]) ** 2
            + W["xi"] * X[1, 1:] ** 2
            + W["ddelta"] * U[0, :] ** 2
            + W["dp"] * dp**2
        )
        J = ca.dot(ds, cost) + W["n_f"] * X[0, N] ** 2 + W["xi_f"] * X[1, N] ** 2
        J += W["slack"] * (ca.sum1(ca.vec(S)) + ca.sumsqr(S))
        opti.minimize(J)

        opts = {"ipopt.print_level": 0, "ipopt.sb": "yes", "ipopt.max_iter": 100, "ipopt.mu_init": 1e-3,
                "ipopt.bound_push": 1e-6, "ipopt.bound_frac": 1e-6, "print_time": False, "expand": True}
        opti.solver("ipopt", opts)
        # the whole problem as one function: parameters and initial guess in, solution out
        self.solver = opti.to_function("mpc", [x0, p0, ref, X, U, S], [X, U, S])
        self.guess = None

    def references(self, path, s0):
        """Abscissa of the nodes and the reference rows of the problem."""
        z = [s0]
        for _ in range(self.N):
            z.append(z[-1] + self.dt * max(path("v", z[-1]), 5.0))
        z = np.array(z)
        half = self.veh["width"] / 2
        ref = np.vstack([
            np.diff(z),
            path("kappa", (z[:-1] + z[1:]) / 2),
            path("v", z[1:]),
            path("n_right", z[1:]) + half,
            path("n_left", z[1:]) - half,
        ])
        return z, ref

    def initial_guess(self, ref):
        """On the line at the reference speed, steering for steady cornering: delta = k (L + K_us u^2)."""
        k, u = ref[1], ref[2]
        L = self.veh["lf"] + self.veh["lr"]
        X = np.vstack([0 * u, 0 * u, u, 0 * u, u * k, k * (L + understeer_gradient(self.veh) * u**2)])
        return np.hstack([X[:, :1], X]), np.zeros((2, self.N)), np.zeros((4, self.N))

    def solve(self, path, x0, s0, p0):
        """Input to apply now, and the details of the solution."""
        z, ref = self.references(path, s0)
        if self.guess is None:
            Xg, Ug, Sg = self.initial_guess(ref)
        else:  # warm start: previous solution, moved to the new nodes
            zp, Xp, Up, Sp = self.guess
            shift = lambda A, zz, zn: np.array([np.interp(zn, zz, a) for a in A])
            Xg, Ug, Sg = shift(Xp, zp, z), shift(Up, zp[:-1], z[:-1]), shift(Sp, zp[:-1], z[:-1])
        Xg[:, 0] = x0

        t0 = time.perf_counter()
        X, U, S = (np.array(a) for a in self.solver(x0, p0, ref, Xg, Ug, Sg))
        cpu = time.perf_counter() - t0
        self.guess = (z, X, U, S)
        return U[:, 0], dict(z=z, X=X, U=U, cpu=cpu, ok=self.solver.stats()["success"])
