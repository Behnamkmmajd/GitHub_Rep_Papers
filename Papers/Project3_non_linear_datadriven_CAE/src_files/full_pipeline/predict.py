# Full Pipeline Prediction — (mu, t) → MLP → z → Decoder → u(x,y)
# ============================================================
# Predicts unseen test snapshots using the trained MLP + decoder.
# Saves results as NPZ for MATLAB plotting.
#
# Usage:  python predict.py   (run from full_pipeline/)
# ////////////////////////////////////////////////////////////
# USER CONFIGURATION — modify only this section
# ////////////////////////////////////////////////////////////

# --- Dataset ---
DATASET_FOLDER = "Advection"
CHANNELS       = [1]

# --- Encoder / decoder type ---
ENCODER_TYPE   = "KTH"                  # "CNN" or "KTH"

# CNN settings (used when ENCODER_TYPE == "CNN")
LATENT_DIM     = 4
CNN_RUN_NAME   = "Advection_Ch1_Nt240_B32_LD4_LR0.0001_E5000_Norm"

# KTH settings (used when ENCODER_TYPE == "KTH")
DIM            = 16
DIM_MULTS      = (1, 2, 4, 8, 16, 32)
Z_CHANNELS     = 4
BLOCK_TYPE     = 1
KTH_RUN_NAME   = "Advection_Ch1_Nt240_B8_NL6_LR1e-05_L21e-06_E300_Norm"

# --- MLP run ---
MLP_RUN_NAME   = "Advection_Ch1_KTH_MLP_HD64_LD1024_LR0.001_E5000"
HIDDEN_DIM     = 64

# --- Snapshot selection ---
# Set to None to predict ALL test snapshots, or pick specific
# (case_index, time_index) pairs from test_cases / timestep_values.
# Example: [(0, 0), (0, 8), (5, 0), (5, 8)]
#   case_index = position in test_cases list (0 = first test case)
#   time_index = position in timestep_values  (0 = first timestep)
SAMPLE_PAIRS   = [(3, 12), (5, 14), (10, 5), (32, 5)]   # None = all
# SAMPLE_PAIRS = None

# --- Prediction batch size ---
PRED_BATCH     = 64

# ////////////////////////////////////////////////////////////
# END OF USER CONFIGURATION
# ////////////////////////////////////////////////////////////

import sys, os, json
import numpy as np
import torch
import h5py

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# ============================================================
# Import models from source projects
# ============================================================
sys.path.append(os.path.join(PROJECT_ROOT, 'latent_mlp', 'library_latent_mlp'))
from latent_mlp.models.latent_mlp import LatentMLP

if ENCODER_TYPE == "CNN":
    sys.path.append(os.path.join(PROJECT_ROOT, 'simple_cae',
                                 'library_simple_cae'))
    from simple_cae.models.cae_model import ConvolutionalAutoencoder
    DECODER_MODEL_PATH = os.path.join(PROJECT_ROOT, 'simple_cae',
                                      'runs', CNN_RUN_NAME, 'model_final.pt')
elif ENCODER_TYPE == "KTH":
    sys.path.append(os.path.join(PROJECT_ROOT, 'resnet_cae', 'library_resnet_cae'))
    from resnet_cae.models.baseline_model import ConvAutoencoderBaseline
    DECODER_MODEL_PATH = os.path.join(PROJECT_ROOT, 'resnet_cae',
                                      'runs', KTH_RUN_NAME, 'model_final.pt')
else:
    raise ValueError(f"Unknown ENCODER_TYPE: {ENCODER_TYPE!r}")

# ============================================================
# Paths
# ============================================================
DATABASE_DIR   = os.path.join(PROJECT_ROOT, 'database', DATASET_FOLDER)
SPLIT_FILE     = os.path.join(DATABASE_DIR, 'split.json')
MLP_RUN_DIR    = os.path.join(PROJECT_ROOT, 'latent_mlp', 'runs', MLP_RUN_NAME)
MLP_MODEL_PATH = os.path.join(MLP_RUN_DIR, 'model_final.pt')
NORM_PATH      = os.path.join(MLP_RUN_DIR, 'normalization_params.npz')

# Output directory
OUTPUT_DIR = os.path.join('runs', f'{DATASET_FOLDER}_{ENCODER_TYPE}_prediction')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Device
# ============================================================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[INFO] Device: {device}")

# ============================================================
# Load split.json
# ============================================================
with open(SPLIT_FILE, 'r') as f:
    split = json.load(f)

timesteps_per_case = split['timesteps_per_case']
timestep_values    = np.array(split['timestep_values'])
param_names        = split['param_names']
n_params           = len(param_names)
n_channels         = len(CHANNELS)

# Combine val + test into one unseen pool
_val_cases  = split['val_cases']
_val_params = split['val_params']
_test_cases = split['test_cases']
_test_params = split['test_params']

test_cases  = _val_cases + _test_cases
test_params = _val_params + _test_params
test_indices = split['val_indices'] + split['test_indices']

n_val_pool  = len(_val_cases)
n_test_pool = len(_test_cases)
print(f"[INFO] Unseen pool: {len(test_cases)} cases "
      f"({n_val_pool} val + {n_test_pool} test), "
      f"snapshots (all): {len(test_indices)}")

# --- Apply snapshot selection ---
if SAMPLE_PAIRS is not None:
    sel_indices = []
    sel_positions = []    # positions within test_indices
    for ci, ti in SAMPLE_PAIRS:
        case = test_cases[ci]
        global_idx = case * timesteps_per_case + ti
        pos = test_indices.index(global_idx)
        sel_indices.append(global_idx)
        sel_positions.append(pos)
    test_indices = sel_indices
    # Rebuild test_cases/test_params to only include selected cases
    used_cases = sorted(set(test_cases[ci] for ci, ti in SAMPLE_PAIRS))
    test_params = [test_params[test_cases.index(c)] for c in used_cases]
    test_cases  = used_cases
    print(f"[INFO] Selected {len(test_indices)} snapshots: {SAMPLE_PAIRS}")
else:
    print(f"[INFO] Predicting all {len(test_indices)} test snapshots")

# ============================================================
# Load ground-truth fields
# ============================================================
arrays = []
for ch in CHANNELS:
    with h5py.File(os.path.join(DATABASE_DIR, f'channel_{ch}.h5'), 'r') as f:
        arrays.append(f['field'][test_indices])
U_true = np.stack(arrays, axis=1).astype(np.float32)   # (N_test, C, H, W)
print(f"[INFO] U_true: {U_true.shape}")

# --- Normalize ground truth if encoder was trained with normalization ---
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

mlp_ri = json.load(open(os.path.join(MLP_RUN_DIR, 'run_info.json')))
field_normalize = mlp_ri.get('field_normalize', False)
if field_normalize:
    U_true = normalize_per_snapshot(U_true)
    print("[INFO] Normalized ground truth (per-snapshot min-max [0,1])")

# ============================================================
# Build (mu..., t) input array for test set
# ============================================================
N_test = len(test_indices)
X_test = np.zeros((N_test, n_params + 1), dtype=np.float32)
for i, idx in enumerate(test_indices):
    case_id   = idx // timesteps_per_case
    t_idx     = idx % timesteps_per_case
    local_pos = test_cases.index(case_id)
    X_test[i, :n_params] = test_params[local_pos]
    X_test[i, n_params]  = timestep_values[t_idx]

# Normalize inputs using training statistics
norm = np.load(NORM_PATH)
X_test_norm = (X_test - norm['input_mean']) / norm['input_std']

print(f"[INFO] X_test: {X_test.shape}  (columns: {param_names + ['t']})")

# ============================================================
# Load MLP
# ============================================================
latent_dim_actual = mlp_ri['latent_dim']
input_dim         = mlp_ri['input_dim']

mlp = LatentMLP(input_dim=input_dim, hidden_dim=HIDDEN_DIM,
                latent_dim=latent_dim_actual).to(device)
mlp.load_state_dict(torch.load(MLP_MODEL_PATH, map_location=device, weights_only=True))
mlp.eval()
print(f"[INFO] MLP loaded: input_dim={input_dim}, latent_dim={latent_dim_actual}")

# ============================================================
# Load Decoder
# ============================================================
if ENCODER_TYPE == "CNN":
    decoder = ConvolutionalAutoencoder(latent_dim=LATENT_DIM,
                                       in_channels=n_channels).to(device)
    sd = torch.load(DECODER_MODEL_PATH, map_location=device, weights_only=True)
    decoder.load_state_dict(sd)
else:
    decoder = ConvAutoencoderBaseline(
        dim=DIM, dim_mults=DIM_MULTS, channels=n_channels,
        z_channels=Z_CHANNELS, block_type=BLOCK_TYPE,
    ).float().to(device)
    sd = torch.load(DECODER_MODEL_PATH, map_location=device, weights_only=True)
    if list(sd.keys())[0].startswith("module."):
        sd = {k.replace("module.", "", 1): v for k, v in sd.items()}
    decoder.load_state_dict(sd, strict=True)

decoder.eval()
print(f"[INFO] Decoder loaded: {ENCODER_TYPE}")

# ============================================================
# Predict: (mu, t) → MLP → z → Decoder → u_pred
# ============================================================
print("\n[INFO] Predicting test snapshots...")

# Infer spatial latent shape for KTH (needed to un-flatten z)
if ENCODER_TYPE == "KTH":
    with torch.no_grad():
        dummy = torch.zeros(1, n_channels, *split['spatial_shape'],
                            dtype=torch.float32, device=device)
        z_shape = decoder.encode(dummy).shape[1:]   # (z_ch, H', W')
    print(f"[INFO] KTH latent shape: {z_shape}")

X_test_t = torch.tensor(X_test_norm, dtype=torch.float32)
U_pred_list = []

with torch.no_grad():
    for i in range(0, N_test, PRED_BATCH):
        x_batch = X_test_t[i:i+PRED_BATCH].to(device)
        z_flat  = mlp(x_batch)                       # (B, latent_dim_actual)

        if ENCODER_TYPE == "KTH":
            z = z_flat.reshape(-1, *z_shape)          # un-flatten to spatial
        else:
            z = z_flat

        u_batch = decoder.decode(z)                   # (B, C, H, W)
        U_pred_list.append(u_batch.cpu().numpy())

U_pred = np.concatenate(U_pred_list, axis=0)          # (N_test, C, H, W)
print(f"[INFO] U_pred: {U_pred.shape}")

# ============================================================
# Save results
# ============================================================
np.savez(os.path.join(OUTPUT_DIR, 'predictions.npz'),
         U_true=U_true,                    # (N_test, C, H, W)
         U_pred=U_pred,                    # (N_test, C, H, W)
         X_test=X_test,                    # (N_test, n_params+1)  [mu..., t]
         test_cases=np.array(test_cases),
         timesteps_per_case=timesteps_per_case)

# Per-snapshot error metrics
err = U_true - U_pred
mae  = np.mean(np.abs(err), axis=(1, 2, 3))
mse  = np.mean(err**2,      axis=(1, 2, 3))
linf = np.max(np.abs(err),  axis=(1, 2, 3))

print(f"\n  MAE  — mean: {mae.mean():.6f}, max: {mae.max():.6f}")
print(f"  MSE  — mean: {mse.mean():.6f}, max: {mse.max():.6f}")
print(f"  Linf — mean: {linf.mean():.6f}, max: {linf.max():.6f}")

np.savez(os.path.join(OUTPUT_DIR, 'error_metrics.npz'),
         mae=mae, mse=mse, linf=linf,
         X_test=X_test)

print(f"\n[INFO] All results saved to: {OUTPUT_DIR}")
