"""Convert a track of the TUM racetrack database (https://github.com/TUMFTM/racetrack-database, LGPL-3.0) into
the track format of this repository (see track.py):

    uv run tracks/tumftm.py Silverstone    ->  tracks/silverstone.csv

The database gives the centre line every 5 m with the track widths. The closed centre line is written as a
Fourier series, smoothed with a Gaussian of sigma metres, and sampled every metre; heading and curvature come
from the derivatives of the series.
"""

import sys
import urllib.request

import numpy as np

URL = "https://raw.githubusercontent.com/TUMFTM/racetrack-database/master/tracks/{}.csv"


def convert(name, sigma=2.0, ds=1.0):
    with urllib.request.urlopen(URL.format(name)) as f:
        x, y, w_right, w_left = np.loadtxt(f, delimiter=",", comments="#").T
    n, m = len(x), 8 * len(x)  # points of the database, points of the fine grid
    L = np.sum(np.hypot(np.diff(x, append=x[0]), np.diff(y, append=y[0])))
    k = np.arange((n + 1) // 2)  # harmonics below the Nyquist frequency
    gauss = np.exp(-0.5 * (2 * np.pi * k / L * sigma) ** 2)

    def fine(z, d=0):
        """Smoothed periodic series of z, or its d-th derivative w.r.t. the point index, on the fine grid."""
        C = np.zeros(m // 2 + 1, complex)
        C[: len(k)] = np.fft.rfft(z)[: len(k)] * gauss * (2j * np.pi * k / n) ** d
        return np.fft.irfft(C, m) * m / n

    dx, dy, ddx, ddy = fine(x, 1), fine(y, 1), fine(x, 2), fine(y, 2)
    theta = np.unwrap(np.arctan2(dy, dx))
    turn = 2 * np.pi * np.round((theta[-1] - theta[0]) / (2 * np.pi))
    columns = {
        "x": fine(x),
        "y": fine(y),
        "theta": theta,
        "kappa": (dx * ddy - dy * ddx) / (dx**2 + dy**2) ** 1.5,
        "w_left": fine(w_left),
        "w_right": fine(w_right),
    }
    # close the lap and sample every ds metres of arc length
    speed = np.hypot(dx, dy) * n / m
    s = np.r_[0, np.cumsum((speed + np.roll(speed, -1)) / 2)]
    columns = {key: np.r_[val, val[0] + (turn if key == "theta" else 0)] for key, val in columns.items()}
    s_out = np.linspace(0, s[-1], int(np.round(s[-1] / ds)) + 1)
    return {"s": s_out} | {key: np.interp(s_out, s, val) for key, val in columns.items()}


if __name__ == "__main__":
    name = sys.argv[1]
    print(f"Downloading {name} from the TUM racetrack database")
    track = convert(name)
    file = f"tracks/{name.lower()}.csv"
    fmt = ["%.3f", "%.3f", "%.3f", "%.6f", "%.7f", "%.3f", "%.3f"]
    np.savetxt(file, np.column_stack(list(track.values())), fmt=fmt, delimiter=",", header=",".join(track),
               comments="")
    print(f"{name}: {track['s'][-1]:.1f} m, minimum radius {1 / np.abs(track['kappa']).max():.1f} m -> {file}")
