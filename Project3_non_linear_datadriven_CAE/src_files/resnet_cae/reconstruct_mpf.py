# CAE Reconstruction Script — MPF Dataset (KTH ResNet)
# ============================================================
# Reconstructs MPF snapshots through the old KTH ResNet CAE
# using the EXACT normalization from the original training.
#
# The MPF encoder requires per-snapshot G normalization and
# per-case U/V fluctuation normalization (not raw HDF5).
# This script replicates that pipeline, reconstructs, and
# saves results as NPZ for MATLAB post-processing.
#
# Usage:  python reconstruct_mpf.py   (run from resnet_cae/)
# ============================================================

# ////////////////////////////////////////////////////////////
# USER CONFIGURATION — modify only this section
# ////////////////////////////////////////////////////////////

# --- Encoder (must match the trained model) ---
DIM            = 16
DIM_MULTS      = (1, 2, 4, 8, 16, 32, 64, 128)   # 8 levels
Z_CHANNELS     = 4
BLOCK_TYPE     = 1
ENCODER_CHANNELS = 3                                # G, U, V

# --- Run to load ---
RUN_NAME = "2DMPF_GUV_Nt2424_B16_NL8_LR1e-05_L21e-06_E300"

# --- Test snapshots to reconstruct ---
# Each entry is (case_source, case_index, timestep_index)
#   case_source:  "train", "val", or "test"
#   case_index:   0-based index into split.json {source}_cases
#   timestep_index: 0-based timestep within that case (0..100)
TEST_SAMPLES = [
    ("train",  0,   0),    # 1st train case, t=0
    ("train",  0,  50),    # 1st train case, t=50
    ("train",  0, 100),    # 1st train case, t=100
    ("val",    0,   0),    # 1st val case, t=0
    ("val",    0,  50),    # 1st val case, t=50
    ("val",    0, 100),    # 1st val case, t=100
    ("test",   0,   0),    # 1st test case, t=0
    ("test",   0,  50),    # 1st test case, t=50
    ("test",   0, 100),    # 1st test case, t=100
    ("test",   1,  50),    # 2nd test case, t=50
]

# ////////////////////////////////////////////////////////////
# END OF USER CONFIGURATION
# ////////////////////////////////////////////////////////////

import sys, os
import json
import numpy as np
import torch
import h5py

sys.path.append(os.path.join(os.path.dirname(__file__), 'library_resnet_cae'))
from resnet_cae.models.baseline_model import ConvAutoencoderBaseline

# --- Paths ---
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATABASE_DIR = os.path.join(PROJECT_ROOT, 'database', 'MPF')
SPLIT_FILE   = os.path.join(DATABASE_DIR, 'split.json')
RUN_DIR      = os.path.join('runs', RUN_NAME)
MODEL_PATH   = os.path.join(RUN_DIR, 'model_final.pt')

# --- Load split ---
with open(SPLIT_FILE, 'r') as f:
    split = json.load(f)

timesteps_per_case = split['timesteps_per_case']

# ============================================================
# MPF Normalization (per-case)
# ============================================================

def normalize_case_mpf(case_idx, database_dir, tpc):
    """
    Load one MPF case from HDF5, apply the original normalization.
    Returns: np.ndarray (tpc, 3, 256, 512) float32 in [0, 1].
    """
    start = case_idx * tpc

    with h5py.File(os.path.join(database_dir, 'channel_1.h5'), 'r') as f:
        c_all = f['field'][start:start + tpc, :, :]
    with h5py.File(os.path.join(database_dir, 'channel_2.h5'), 'r') as f:
        u_all = f['field'][start:start + tpc, :, :]
    with h5py.File(os.path.join(database_dir, 'channel_3.h5'), 'r') as f:
        v_all = f['field'][start:start + tpc, :, :]

    X = np.zeros((tpc, 3, 256, 512), dtype=np.float32)

    # ---- G channel: per-snapshot to [0, 1] ----
    for t in range(tpc):
        G = c_all[t].T.astype(np.float32)
        g_min = np.min(G)
        if g_min < 0.0:
            GNor = np.abs(g_min) + G
            GNor = GNor / np.max(GNor)
        elif g_min > 0.0:
            GNor = G - np.abs(g_min)
            GNor = GNor / np.max(GNor)
        else:
            GNor = G / np.max(G)
        X[t, 0] = GNor

    # ---- U, V channels: per-case fluctuation to [0, 1] ----
    ns = 512 * 256
    U_flat = u_all.reshape(tpc, ns).T
    V_flat = v_all.reshape(tpc, ns).T

    UP = U_flat - np.mean(U_flat, axis=1, keepdims=True)
    VP = V_flat - np.mean(V_flat, axis=1, keepdims=True)

    UNor = (UP - UP.min()) / (UP.max() - UP.min() + 1e-8)
    VNor = (VP - VP.min()) / (VP.max() - VP.min() + 1e-8)

    for t in range(tpc):
        X[t, 1] = UNor[:, t].reshape(512, 256).T
        X[t, 2] = VNor[:, t].reshape(512, 256).T

    return X


# ============================================================
# Resolve samples
# ============================================================
source_map = {
    'train': (split['train_cases'], split.get('train_params', [])),
    'val':   (split['val_cases'],   split.get('val_params', [])),
    'test':  (split['test_cases'],  split.get('test_params', [])),
}

sample_info = []
for source, ci, ti in TEST_SAMPLES:
    cases, params = source_map[source]
    case_id = cases[ci]
    p = params[ci] if params else []
    sample_info.append({
        'source': source,
        'source_index': ci,
        'case_id': case_id,
        'timestep': ti,
        'params': p,
    })

print(f"[INFO] Dataset: MPF")
print(f"[INFO] Model: {MODEL_PATH}")
print(f"[INFO] Reconstructing {len(sample_info)} snapshots:")
for s in sample_info:
    print(f"  {s['source']}[{s['source_index']}] case={s['case_id']} "
          f"t={s['timestep']} params={s['params']}")

# ============================================================
# Load model
# ============================================================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = ConvAutoencoderBaseline(
    dim=DIM, dim_mults=DIM_MULTS, channels=ENCODER_CHANNELS,
    z_channels=Z_CHANNELS, block_type=BLOCK_TYPE,
).float().to(device)

sd = torch.load(MODEL_PATH, map_location=device, weights_only=True)
if list(sd.keys())[0].startswith("module."):
    sd = {k.replace("module.", "", 1): v for k, v in sd.items()}
model.load_state_dict(sd, strict=True)
model.eval()
print(f"[INFO] Model loaded on {device}")

# ============================================================
# Normalize & Reconstruct (case by case, to reuse normalization)
# ============================================================
# Group samples by case_id so we normalize each case only once
from collections import defaultdict

case_samples = defaultdict(list)   # case_id -> [(out_idx, timestep), ...]
for i, s in enumerate(sample_info):
    case_samples[s['case_id']].append((i, s['timestep']))

n_samples = len(sample_info)
originals       = np.zeros((n_samples, 3, 256, 512), dtype=np.float32)
reconstructions = np.zeros((n_samples, 3, 256, 512), dtype=np.float32)

for case_id, entries in case_samples.items():
    print(f"  Processing case {case_id} ({len(entries)} snapshots)...")
    X_case = normalize_case_mpf(case_id, DATABASE_DIR, timesteps_per_case)

    timesteps = [t for _, t in entries]
    X_batch = X_case[timesteps]

    x = torch.tensor(X_batch, dtype=torch.float32).to(device)
    with torch.no_grad():
        x_recon = model(x).cpu().numpy()

    for j, (out_idx, t) in enumerate(entries):
        originals[out_idx]       = X_batch[j]
        reconstructions[out_idx] = x_recon[j]

# ============================================================
# Metrics
# ============================================================
channel_names = ['G (gas)', 'U (vel-x)', 'V (vel-y)']
print("\n  Per-sample reconstruction quality:")
for i, s in enumerate(sample_info):
    mse = np.mean((originals[i] - reconstructions[i]) ** 2)
    mae = np.mean(np.abs(originals[i] - reconstructions[i]))
    print(f"  [{i}] {s['source']}[{s['source_index']}] t={s['timestep']:3d}: "
          f"MSE={mse:.4e}, MAE={mae:.4e}")

print("\n  Per-channel metrics (all samples):")
for ch in range(3):
    mse = np.mean((originals[:, ch] - reconstructions[:, ch]) ** 2)
    mae = np.mean(np.abs(originals[:, ch] - reconstructions[:, ch]))
    print(f"  Ch{ch} ({channel_names[ch]}): MSE={mse:.4e}, MAE={mae:.4e}")

# ============================================================
# Save NPZ
# ============================================================
out_path = os.path.join(RUN_DIR, 'reconstruction_test.npz')
np.savez(out_path,
         originals=originals,                                       # (N, 3, 256, 512)
         reconstructions=reconstructions,                           # (N, 3, 256, 512)
         channel_names=np.array(channel_names),
         case_ids=np.array([s['case_id'] for s in sample_info]),
         timesteps=np.array([s['timestep'] for s in sample_info]),
         sources=np.array([s['source'] for s in sample_info]),
         params=np.array([s['params'] for s in sample_info]),
         )

print(f"\n[INFO] Saved: {out_path}")
print(f"  Shape: originals={originals.shape}, reconstructions={reconstructions.shape}")
