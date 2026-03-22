# CAE Reconstruction Script — CNN-ROM (Wang et al. 2024)
# ============================================================
# Reconstructs user-selected test snapshots through the trained
# CAE and saves results as NPZ for MATLAB post-processing.
#
# Usage:  python reconstruct.py   (run from simple_cae/)
# ============================================================

# ////////////////////////////////////////////////////////////
# USER CONFIGURATION — modify only this section
# ////////////////////////////////////////////////////////////

# --- Dataset ---
DATASET_FOLDER = "Advection"            # Subfolder inside database/
CHANNELS       = [1]                    # Which channel files to load

# --- Model ---
LATENT_DIM     = 4
RUN_NAME       = "Advection_Ch1_Nt240_B32_LD4_LR0.0001_E5000_Norm"

# --- Snapshots to reconstruct ---
# Each entry is (case_index, timestep_index)
#   case_index:     0-based index into the combined unseen pool (val + test)
#   timestep_index: 0-based timestep within that case
SAMPLES = [
    (0,  0),
    (0,  8),
    (16, 0),
    (32, 8),
    (48, 0),
]

# ////////////////////////////////////////////////////////////
# END OF USER CONFIGURATION
# ////////////////////////////////////////////////////////////

import sys, os
import json
import numpy as np
import torch
import h5py

sys.path.append(os.path.join(os.path.dirname(__file__), 'library_simple_cae'))
from simple_cae.models.cae_model import ConvolutionalAutoencoder

# --- Paths ---
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATABASE_DIR = os.path.join(PROJECT_ROOT, 'database', DATASET_FOLDER)
SPLIT_FILE   = os.path.join(DATABASE_DIR, 'split.json')
RUN_DIR      = os.path.join('runs', RUN_NAME)
MODEL_PATH   = os.path.join(RUN_DIR, 'model_final.pt')

# --- Load split ---
with open(SPLIT_FILE, 'r') as f:
    split = json.load(f)

n_channels         = len(CHANNELS)
timesteps_per_case = split['timesteps_per_case']

# Combine val + test into one unseen pool
val_cases   = split['val_cases']
val_params  = split['val_params']
test_cases  = split['test_cases']
test_params = split['test_params']
n_val = len(val_cases)

pool_cases  = val_cases + test_cases
pool_params = val_params + test_params
pool_source = ['val'] * n_val + ['test'] * len(test_cases)

print(f"[INFO] Unseen pool: {len(pool_cases)} cases "
      f"({n_val} val + {len(test_cases)} test)")

# --- Map (case_index, timestep) to global HDF5 row ---
sample_info = []
global_indices = []
for ci, t_idx in SAMPLES:
    case_id = pool_cases[ci]
    global_row = case_id * timesteps_per_case + t_idx
    sample_info.append({
        'case_index': ci,
        'case_id': case_id,
        'timestep': t_idx,
        'global_row': global_row,
        'params': pool_params[ci],
        'source': pool_source[ci],
    })
    global_indices.append(global_row)

print(f"[INFO] Dataset: {split['dataset']}")
print(f"[INFO] Channels: {CHANNELS}")
print(f"[INFO] Model: {MODEL_PATH}")
print(f"[INFO] Reconstructing {len(SAMPLES)} snapshots:")
for s in sample_info:
    print(f"  [{s['source']}] case {s['case_id']} t={s['timestep']} "
          f"params={s['params']}")

# --- Load data from HDF5 ---
arrays = []
for ch in CHANNELS:
    h5_path = os.path.join(DATABASE_DIR, f'channel_{ch}.h5')
    with h5py.File(h5_path, 'r') as f:
        field = f['field'][global_indices]       # (N, H, W)
    arrays.append(field)

originals = np.stack(arrays, axis=1).astype(np.float32)  # (N, C, H, W)
print(f"[INFO] Loaded originals: {originals.shape}")

# --- Normalize (must match training) ---
def normalize_per_snapshot(data):
    """Per-snapshot, per-channel min-max to [0, 1]."""
    N, C = data.shape[:2]
    data_norm = np.empty_like(data)
    for i in range(N):
        for c in range(C):
            smin = data[i, c].min()
            smax = data[i, c].max()
            denom = smax - smin
            if denom < 1e-12:
                data_norm[i, c] = 0.0
            else:
                data_norm[i, c] = (data[i, c] - smin) / denom
    return data_norm

run_info_path = os.path.join(RUN_DIR, 'run_info.json')
do_norm = False
if os.path.exists(run_info_path):
    with open(run_info_path) as f:
        do_norm = json.load(f).get('normalize', False)

if do_norm:
    originals = normalize_per_snapshot(originals)
    print("[INFO] Applied per-snapshot normalization (matching training)")

# --- Load model ---
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = ConvolutionalAutoencoder(latent_dim=LATENT_DIM,
                                 in_channels=n_channels).to(device)
sd = torch.load(MODEL_PATH, map_location=device, weights_only=True)
model.load_state_dict(sd)
model.eval()
print(f"[INFO] Model loaded on {device}")

# --- Reconstruct ---
x = torch.tensor(originals, dtype=torch.float32).to(device)
with torch.no_grad():
    x_recon = model(x)
reconstructions = x_recon.cpu().numpy()

# --- Per-sample metrics ---
for i, s in enumerate(sample_info):
    mse = np.mean((originals[i] - reconstructions[i]) ** 2)
    mae = np.mean(np.abs(originals[i] - reconstructions[i]))
    print(f"  Sample {i}: MSE={mse:.4e}, MAE={mae:.4e}")

# --- Save NPZ ---
out_path = os.path.join(RUN_DIR, 'reconstruction_test.npz')
np.savez(out_path,
         originals=originals,                                    # (N, C, H, W)
         reconstructions=reconstructions,                        # (N, C, H, W)
         case_indices=np.array([s['case_id'] for s in sample_info]),
         timestep_indices=np.array([s['timestep'] for s in sample_info]),
         params=np.array([s['params'] for s in sample_info]),
         sources=np.array([s['source'] for s in sample_info]),
         channels=np.array(CHANNELS))

print(f"\n[INFO] Saved: {out_path}")
