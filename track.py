"""Track and racing line files.

Track (tracks/*.csv), sampled along the centre line, the last row closes the lap:
    s, x, y, theta, kappa, w_left, w_right
    abscissa, position, heading and curvature of the centre line; distance of the borders from it.

Racing line (output/racing_line.csv), written by solve_ocp.py, sampled every 0.5 m along the line itself:
    s, x, y, theta, kappa, v, t, ax, ay, n_left, n_right
    abscissa, position, heading and curvature of the line; speed, time, accelerations;
    lateral position of the left (> 0) and right (< 0) border w.r.t. the line.
"""

import numpy as np

RACING_LINE = "output/racing_line.csv"


def read_csv(file):
    return np.genfromtxt(file, delimiter=",", names=True)


def borders(track):
    """(x, y) of the left and right border."""
    x, y, th = track["x"], track["y"], track["theta"]
    left = x - track["w_left"] * np.sin(th), y + track["w_left"] * np.cos(th)
    right = x + track["w_right"] * np.sin(th), y - track["w_right"] * np.cos(th)
    return left, right


class Path:
    """Racing line: values at any abscissa (wrapped on the lap) and projection of a point on it."""

    def __init__(self, file=RACING_LINE):
        self.data = read_csv(file)
        self.s = self.data["s"]
        self.L = self.s[-1]
        self.ds = self.s[1] - self.s[0]

    def __call__(self, key, s):
        return np.interp(np.mod(s, self.L), self.s, self.data[key])

    def project(self, X, Y, s_guess, window=15.0):
        """Abscissa s and lateral offset n of the point (X, Y), searched within window metres of s_guess.
        s keeps counting after a lap (it is not wrapped)."""
        i0 = int(np.round(np.mod(s_guess, self.L) / self.ds))
        w = int(window / self.ds)
        idx = np.mod(np.arange(i0 - w, i0 + w + 1), len(self.s) - 1)
        x, y, th = self.data["x"][idx], self.data["y"][idx], self.data["theta"][idx]
        j = np.argmin((X - x) ** 2 + (Y - y) ** 2)
        ex, ey = X - x[j], Y - y[j]
        ds = ex * np.cos(th[j]) + ey * np.sin(th[j])
        n = -ex * np.sin(th[j]) + ey * np.cos(th[j])
        s = s_guess + np.mod(self.s[idx[j]] + ds - s_guess + self.L / 2, self.L) - self.L / 2
        return s, n
