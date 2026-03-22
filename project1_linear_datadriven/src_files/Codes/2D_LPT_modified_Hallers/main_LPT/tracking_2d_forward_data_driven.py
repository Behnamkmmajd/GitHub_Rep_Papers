#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel Trajectory Computation with Gas Fraction and velocity Interpolation at new positions (2D, forward in time, data-driven)
- Reads X.csv, Y.csv and Stacked_*.csv (u,v) from data folder
- Reads gas fraction G from
- Builds interpolants for unsteady field (linear/cubic/Lagrangian)
- Integrates trajectories in parallel and saves results to .mat
"""
# =============================================================================
import os
import sys
import numpy as np
import scipy.io as sio
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
# =============================================================================
# Print header
print("\n2D Forward Data-Driven Trajectory Computation with Gas Fraction Interpolation")
print("=====================================================================\n")

# =============================================================================
# PATHS CONFIGURATION
# =============================================================================
cwd = os.getcwd()
parent_dir = os.path.sep.join(cwd.split(os.path.sep)[:-1])
database_root = os.path.sep.join(cwd.split(os.path.sep)[:-3])

# Database paths
case_name = "dataRe40Bo40"

data_dir_G = os.path.join(database_root, 'database', case_name + '_EG', 'data')

data_dir_UV = os.path.join(database_root, 'database', case_name + '_EUV', 'data')

# Output paths
output_run_name = os.path.join("TrajOnly_" + case_name)
output_root = os.path.join(cwd, 'Results/' , output_run_name)

# Tracking points
grid_name = os.path.join("Grid250x750_" + case_name)
tracking_points_path = os.path.join('cust_shape_dataset/', f"{grid_name}.mat")

# =============================================================================
# USER CONFIGURATION
# =============================================================================
# Time settings
t0_index = 0            # first snapshot (inclusive)
tN_index = 405          # last snapshot (inclusive)
dt_snapshot_dim = 0.05  # snapshot spacing (simulation units)

# Physical scaling
# length_scale = 0.0673435
# vel_scale = 0.256898093

Bubbl_D = 0.4
length_scale = Bubbl_D
vel_scale = (9.81 * length_scale) ** 0.5
time_scale = 0.262141
 
# Interpolation methods
method_V = "linear"
method_G = "linear"

# Parallelism
Ncores = 28
prefer_backend = "processes"

# =============================================================================
# MODULE IMPORTS
# =============================================================================
sys.path.append(os.path.join(parent_dir, "subfunctions", "utils"))
sys.path.append(os.path.join(parent_dir, "subfunctions", "integration"))
from ipynb.fs.defs.Interpolant import interpolant_unsteady  # noqa: F401                                        # type: ignore                   
from ipynb.fs.defs.integration_dFdt import integration_dFdt # noqa: F401                                        # type: ignore                                      
from ipynb.fs.defs.velocity import velocity # noqa: F401                                                        # type: ignore

print("Libraries and modules imported")
print(f"Database (UV): {data_dir_UV}")
print(f"Database (G): {data_dir_G}")
print(f"Tracking points: {tracking_points_path}")
print(f"Output: {output_root}\n")

# =============================================================================
# LOAD GRID AND DATABASE
# =============================================================================
# Load grid coordinates
XX = np.loadtxt(os.path.join(data_dir_G, "x.csv"), delimiter=",").reshape(-1) * length_scale
YY = np.loadtxt(os.path.join(data_dir_G, "y.csv"), delimiter=",").reshape(-1) * length_scale
X, Y = np.meshgrid(XX, YY)
nx, ny = len(XX), len(YY)
print(f"Grid: nx={nx}, ny={ny}")

# Time grid
dt_act = dt_snapshot_dim * time_scale
t0 = t0_index * dt_act
tN = tN_index * dt_act
nt = tN_index - t0_index + 1
time = np.linspace(t0 + dt_act, tN, nt - 1, endpoint=True)

# Initialize velocity and gas fraction arrays
U = np.zeros((ny, nx, nt), dtype=np.float64)
V = np.zeros((ny, nx, nt), dtype=np.float64)
G = np.zeros((ny, nx, nt), dtype=np.float64)

# Load snapshots
for i_snap in range(t0_index, tN_index + 1):
    idx = i_snap - t0_index
    
    # Load velocity
    fname_vel_csv = os.path.join(data_dir_UV, f"Stacked_{i_snap+1}.csv")
    buf = np.loadtxt(fname_vel_csv, delimiter=",").reshape(-1)

    U[:, :, idx] = np.transpose(np.reshape(buf[0:nx*ny] * vel_scale, (nx, ny)))
    V[:, :, idx] = np.transpose(np.reshape(buf[nx*ny:2*nx*ny] * vel_scale, (nx, ny)))
    
    # Load gas fraction
    fname_G_csv = os.path.join(data_dir_G, f"Stacked_{i_snap+1}.csv")
    buf_G = np.loadtxt(fname_G_csv, delimiter=",").reshape(-1)

    G[:, :, idx] = np.transpose(np.reshape(buf_G[0:nx*ny], (nx, ny)))

    # plt.contourf(X, Y, G[:, :, idx], levels=50)

    # plt.axis('equal')

    # plt.pause( 0.001)

    print(f"Loaded Stacked_{i_snap+1}.csv")

# Build interpolants
Interpolant_u = interpolant_unsteady(X, Y, U, method_V)
Interpolant_v = interpolant_unsteady(X, Y, V, method_V)
Interpolant_g = interpolant_unsteady(X, Y, G, method_G)
print(f"Interpolants created (velocity: {method_V}, gas: {method_G})")

# Load tracking points
mat_data = sio.loadmat(tracking_points_path)
m_points = mat_data['T'].T
X0_all = m_points[0:2, :] * length_scale
properties = m_points[2:, :]
Npts = X0_all.shape[1]
print(f"Seed points: N={Npts}")

# =============================================================================
# PARALLEL INTEGRATION
# =============================================================================
periodic = [True, True, False]
bool_unsteady = True
time_data = np.linspace(t0, tN, nt, endpoint=True).reshape(1, -1)   # time data for interpolants                       # type: ignore
assert len(time) == nt - 1, "time array length mismatch"


def split_into_chunks(arr, n_chunks):
    """Split array into n_chunks approximately equal parts."""
    k, m = divmod(len(arr), n_chunks)
    return [arr[i*k + min(i, m):(i+1)*k + min(i+1, m)] for i in range(n_chunks)]


def integrate_chunk(idx_chunk):
    """Integrate trajectories for a chunk of particles."""
    Fmap_local, dFdt_local, props_traj = integration_dFdt(
        time=time,
        x=X0_all[:, idx_chunk],
        properties_init=properties[:, idx_chunk],
        X=X, Y=Y,
        Interpolant_u=Interpolant_u,
        Interpolant_v=Interpolant_v,
        Interpolant_g=Interpolant_g,
        periodic=periodic,
        bool_unsteady=bool_unsteady,
        time_data=time_data,
        verbose=False
    )
    return idx_chunk, Fmap_local, dFdt_local, props_traj


# Split particles into chunks
n_chunks = min(Ncores, Npts)
idx_chunks = [ch for ch in split_into_chunks(np.arange(Npts), n_chunks) if len(ch) > 0]

# Run parallel integration
results = Parallel(n_jobs=n_chunks, prefer=prefer_backend, verbose=50)(
    delayed(integrate_chunk)(idx) for idx in idx_chunks
)

# Assemble results
Nt = len(time)
Fmap = np.zeros((Nt, 2, Npts), dtype=float)
dFdt = np.zeros((Nt-1, 2, Npts), dtype=float)
properties_traj = np.zeros((Nt, 3, Npts), dtype=float)

for idx_chunk, Fmap_local, dFdt_local, props_traj in results:                                                       # type: ignore
    Fmap[:, :, idx_chunk] = Fmap_local
    dFdt[:, :, idx_chunk] = dFdt_local
    properties_traj[:, :, idx_chunk] = props_traj

print("Trajectories integrated (parallel).")

# =============================================================================
# SAVE RESULTS
# =============================================================================
os.makedirs(output_root, exist_ok=True)
out_mat = os.path.join(output_root, f"trajectories_{t0_index}_{tN_index}_{grid_name}_V{method_V}_G{method_G}.mat")

sio.savemat(out_mat, {
    "X0traj": Fmap,
    "X0vel": dFdt,
    "time": np.asarray(time),                                                                                       # type: ignore
    "X0_init": X0_all,
    "properties_initial": properties,
    "properties_trajec": properties_traj,
    "XX": XX,
    "YY": YY
})
print(f"Saved: {out_mat}")
