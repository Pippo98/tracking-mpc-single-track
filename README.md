# Tracking MPC

A nonlinear Model Predictive Controller (MPC) that drives a car around Silverstone, following a minimum lap time
racing line. Steering and throttle/brake are computed together, with a dynamic single-track model written
along the path. Everything is written in Python with [CasADi](https://web.casadi.org/) and IPOPT.

The work is split in two steps:

1. **Racing line.** An optimal control problem (OCP) finds the fastest way around the track for a point mass
   with the car's grip and engine limits. The result, the path and the speed to drive, is saved to
   `output/racing_line.csv`.
2. **MPC.** At every control step (0.05 s) the MPC predicts the car over the next 2.5 s and computes the
   steering and throttle/brake that keep it on the racing line at the racing line speed. A simulated car
   drives one lap with it.

## Install

Python 3.11 or newer is required: check with `python3 --version` (the `python3` that comes with macOS is
3.9, too old). Any of the following works:

**uv** ([installation](https://docs.astral.sh/uv/getting-started/installation/)), which also downloads a
suitable Python if needed:

```
uv sync
```

**venv and pip**:

```
python3 -m venv .venv
source .venv/bin/activate              # Windows: .venv\Scripts\activate
pip install casadi numpy matplotlib
```

**conda**:

```
conda create -n tracking-mpc -c conda-forge python=3.12 casadi numpy matplotlib
conda activate tracking-mpc
```

With venv or conda, activate the environment in every new terminal and write `python` instead of `uv run`
in the commands below (e.g. `python solve_ocp.py`).

## Run

```
uv run solve_ocp.py   # step 1: racing line -> output/racing_line.csv   (about 20 s)
uv run run_mpc.py     # step 2: one lap with the MPC                   (about 1 min)
```

Both show a figure at the end, and also save it in `output/`:

| Figure | Content |
|---|---|
| `racing_line.png` | racing line on the track, coloured by speed, and its speed profile |
| `mpc_lap.png` | path driven by the MPC and racing line; speed, throttle/brake and steering along the lap |

Step 2 also prints a summary and saves every control step to `output/mpc_log.csv`. On an Apple M2 Max:

```
lap time 155.07 s (racing line 154.16 s)
lateral error: max 0.072 m, rms 0.008 m
speed error rms 0.575 m/s
solve time: mean 19.9 ms, max 58.0 ms, failed 0
```

The MPC is about 1 s slower than the racing line. The single-track car loses some speed to tyre slip when it
corners under power (the $\Omega v$ term in the $u$ equation, see `notes.md`), and the point mass of step 1
does not have this loss.

## Parameters

All the parameters are in `configs/`:

| File | Content |
|---|---|
| `configs/vehicle.toml` | the car: mass, geometry, tyres, engine, brakes, drag, steering limits |
| `configs/ocp.toml` | step 1: track file, mesh step, share of the grip the racing line may use, margin from the borders |
| `configs/mpc.toml` | step 2: horizon, control period, cost weights |

After changing `vehicle.toml` or `ocp.toml`, run step 1 again before step 2.

## Files

| File | Content |
|---|---|
| `solve_ocp.py` | step 1: minimum lap time OCP of a point mass, conversion to the racing line |
| `run_mpc.py` | step 2: closed-loop simulation of one lap |
| `mpc.py` | the MPC: optimal control problem, references along the line, warm start |
| `vehicle.py` | vehicle parameters, dynamic single-track model with linear tyres |
| `track.py` | track and racing line files, projection of the car on the racing line |
| `plots.py` | figures |
| `tracks/silverstone.csv` | Silverstone centre line and track width, every 1 m |
| `tracks/tumftm.py` | conversion of a track of the TUM racetrack database (see below) |
| `notes.md` | the model and the MPC formulation, with references to the course slides |

## How it works

**Step 1** (`solve_ocp.py`). The car is a point mass with speed $V$, at lateral distance $n$ from the track
centre line and with heading $\xi$ relative to it. Its controls are the accelerations $a_x$, $a_y$, which must
stay inside the friction circle (90% of $\mu g$) and the engine torque and power limits. The equations are
written with the centre line abscissa $s$ as independent variable, so that the track borders are simple
bounds on $n$. The OCP minimises the lap time over one lap, discretised every 2 m (about 2900 intervals),
and IPOPT solves it in about 10 s.

A friction circle is the simplest grip model. [MLTS-point-mass](https://github.com/GiacomoCorradini/MLTS-point-mass),
by Giacomo Corradini, solves the same minimum lap time problem with more complex G-G-V envelopes: any
speed-dependent G-G diagram of the vehicle (e.g. computed from a full vehicle model), with separate grip
scales for braking, traction and cornering.

**Step 2** (`run_mpc.py`, `mpc.py`). The MPC model is the dynamic single-track model with linear tyres,
written along the racing line:
state $[n, \xi, u, v, \Omega, \delta]$ (lateral and heading error from the line, longitudinal and lateral
speed, yaw rate, steering angle), inputs $[\dot\delta, p]$ (steering rate, throttle/brake in $[-1, 1]$).
The cost penalises the lateral, heading and speed errors and the input variations. The tyre slip angles, the
friction circle and the track borders are soft constraints. The simulated car is the same single-track
model, integrated in time. See `notes.md` for the equations.

## Track

`tracks/silverstone.csv` is the Silverstone Grand Prix circuit from the
[TUM racetrack database](https://github.com/TUMFTM/racetrack-database) (LGPL-3.0, see
`tracks/LICENSE-TUMFTM`): centre line from OpenStreetMap, track widths from satellite images.
Only the track geometry comes from the database: the racing line is computed by step 1 (the database's own
race lines are not used), and the lap by the MPC.
`tracks/tumftm.py` downloads a track of the database, smooths its centre line and computes heading and
curvature. Any other track of the database works the same way, e.g.:

```
uv run tracks/tumftm.py Spielberg   # -> tracks/spielberg.csv, then set it as track in configs/ocp.toml
```
