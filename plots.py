"""Figures of the racing line (step 1) and of the MPC lap (step 2), also saved in output/."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

from track import borders

BLUE, ORANGE = "#2a78d6", "#eb6834"
BORDER = "#898781"
SPEED_CMAP = LinearSegmentedColormap.from_list("speed", ["#86b6ef", "#2a78d6", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb",
    "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": "#52514e",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#e1e0d9",
    "grid.linewidth": 0.6,
    "xtick.color": "#52514e",
    "ytick.color": "#52514e",
    "text.color": "#0b0b0b",
    "lines.linewidth": 1.5,
    "legend.frameon": False,
})


def track_map(ax, track):
    for x, y in borders(track):
        ax.plot(x, y, color=BORDER, lw=0.8)
    ax.plot(track["x"][0], track["y"][0], "o", color="#0b0b0b", ms=6, zorder=3)
    ax.annotate("start", (track["x"][0], track["y"][0]), xytext=(8, 8), textcoords="offset points")
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")


def speed_line(ax, x, y, v):
    """Line coloured by the speed v [m/s], with its colour bar in km/h."""
    points = np.column_stack([x, y])
    segments = np.stack([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=SPEED_CMAP, linewidths=2.5, capstyle="round")
    lc.set_array(v[:-1] * 3.6)
    ax.add_collection(lc)
    plt.colorbar(lc, ax=ax, label="speed [km/h]", shrink=0.8)


def save(fig, file):
    fig.tight_layout()
    fig.savefig(file, dpi=150)
    print(f"Figure -> {file}")


def racing_line(track, line):
    fig, ax = plt.subplots(1, 2, figsize=(14, 7), width_ratios=[1, 1.4])
    fig.suptitle(f"Racing line: lap time {line['t'][-1]:.2f} s")
    track_map(ax[0], track)
    speed_line(ax[0], line["x"], line["y"], line["v"])
    ax[1].plot(line["s"], line["v"] * 3.6, color=BLUE)
    ax[1].set_xlabel("s [m]")
    ax[1].set_ylabel("speed [km/h]")
    save(fig, "output/racing_line.png")


def mpc_lap(track, path, log):
    """Path on the left; speed, throttle / brake and steering along the lap on the right."""
    s = log["s"]
    fig, ax = plt.subplot_mosaic(
        [["path", "speed"], ["path", "pedal"], ["path", "steer"]], figsize=(15, 8), width_ratios=[1, 1.6]
    )
    fig.suptitle("MPC lap")

    track_map(ax["path"], track)
    ax["path"].plot(path.data["x"], path.data["y"], color=ORANGE, lw=4, label="racing line")
    ax["path"].plot(log["X"], log["Y"], color=BLUE, lw=1.2, label="MPC")
    ax["path"].legend(loc="upper left")

    ax["speed"].plot(s, path("v", s) * 3.6, color=ORANGE, label="racing line")
    ax["speed"].plot(s, log["u"] * 3.6, color=BLUE, label="MPC")
    ax["speed"].set_ylabel("speed [km/h]")
    ax["speed"].legend(loc="lower left")
    ax["pedal"].plot(s, log["p"], color=BLUE)
    ax["pedal"].set_ylabel("throttle (> 0) /\nbrake (< 0) [-]")
    ax["pedal"].set_ylim(-1.1, 1.1)
    ax["steer"].plot(s, np.rad2deg(log["delta"]), color=BLUE)
    ax["steer"].set_ylabel("steering angle\n(front wheel) [deg]")
    ax["steer"].set_xlabel("s [m]")
    for key in ("pedal", "steer"):
        ax[key].sharex(ax["speed"])
    for key in ("speed", "pedal"):
        ax[key].tick_params(labelbottom=False)

    save(fig, "output/mpc_lap.png")
