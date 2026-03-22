# Full Pipeline Prediction — MPF Dataset
# ============================================================
# (DR, VR, Re, Bo, t) → MLP → z → Decoder → u_pred (normalized)
#
# The MPF encoder was trained with custom per-case normalization,
# so predictions are in normalized [0,1] space.  We normalize the
# ground truth identically for a fair comparison.
#
# Usage:  python predict_mpf.py   (run from full_pipeline/)
# ============================================================

# ////////////////////////////////////////////////////////////
# USER CONFIGURATION
# ////////////////////////////////////////////////////////////

DATASET_FOLDER = "MPF"

# --- KTH ResNet CAE (decoder) ---
DIM            = 16
DIM_MULTS      = (1, 2, 4, 8, 16, 32, 64, 128)
Z_CHANNELS     = 4
BLOCK_TYPE     = 1
ENCODER_CHANNELS = 3
KTH_RUN_NAME   = "2DMPF_GUV_Nt2424_B16_NL8_LR1e-05_L21e-06_E300"

# --- MLP ---
MLP_RUN_NAME   = "MPF_KTH_MLP_HD128_LD128_LR0.001_E5000"
HIDDEN_DIM     = 128

# --- Snapshot selection ---
# (case_index_in_test_list, time_index)
# test_cases = [3831, 3864]  →  case_index 0 or 1
# timestep 0..100
SAMPLE_PAIRS   = [
    (0, 0), (0, 25), (0, 50), (0, 75), (0, 100),
    (1, 0), (1, 25), (1, 50), (1, 75), (1, 100),
]
# SAMPLE_PAIRS = None   # uncomment to predict ALL 202 test snapshots

PRED_BATCH     = 4       # small — 450M param decoder

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
# Imports
# ============================================================
sys.path.append(os.path.join(PROJECT_ROOT, 'latent_mlp', 'library_latent_mlp'))
from latent_mlp.models.latent_mlp import LatentMLP

sys.path.append(os.path.join(PROJECT_ROOT, 'resnet_cae', 'library_resnet_cae'))
from resnet_cae.models.baseline_model import ConvAutoencoderBaseline  # type: ignore

# ============================================================
# Paths
# ============================================================
DATABASE_DIR       = os.path.join(PROJECT_ROOT, 'database', DATASET_FOLDER)
SPLIT_FILE         = os.path.join(DATABASE_DIR, 'split.json')
DECODER_MODEL_PATH = os.path.join(PROJECT_ROOT, 'resnet_cae', 'runs',
                                  KTH_RUN_NAME, 'model_final.pt')
MLP_RUN_DIR        = os.path.join(PROJECT_ROOT, 'latent_mlp', 'runs', MLP_RUN_NAME)
MLP_MODEL_PATH     = os.path.join(MLP_RUN_DIR, 'model_final.pt')
NORM_PATH          = os.path.join(MLP_RUN_DIR, 'normalization_params.npz')

OUTPUT_DIR = os.path.join('runs', f'{DATASET_FOLDER}_KTH_prediction')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Device
# ============================================================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[INFO] Device: {device}")

# ============================================================
# Load Split
# ============================================================
with open(SPLIT_FILE, 'r') as f:
    split = json.load(f)

tpc             = split['timesteps_per_case']
timestep_values = np.array(split['timestep_values'], dtype=np.float32)
param_names     = split['param_names']
n_params        = len(param_names)
test_cases      = split['test_cases']
test_params     = split['test_params']

print(f"[INFO] Test cases: {test_cases}")
print(f"[INFO] Test params: {test_params}")

# ============================================================
# Snapshot selection
# ============================================================
if SAMPLE_PAIRS is not None:
    case_time_pairs = [(test_cases[ci], ti) for ci, ti in SAMPLE_PAIRS]
    unique_cases = sorted(set(c for c, _ in case_time_pairs))
    print(f"[INFO] Selected {len(SAMPLE_PAIRS)} snapshots from {len(unique_cases)} cases")
else:
    case_time_pairs = [(c, t) for c in test_cases for t in range(tpc)]
    unique_cases = test_cases[:]
    print(f"[INFO] All {len(case_time_pairs)} test snapshots from {len(unique_cases)} cases")

N_test = len(case_time_pairs)

# ============================================================
# MPF Normalization (same as training)
# ============================================================

def normalize_case_mpf(case_idx):
    """Load one case, apply the original normalization → (tpc, 3, 256, 512)."""
    start = case_idx * tpc

    with h5py.File(os.path.join(DATABASE_DIR, 'channel_1.h5'), 'r') as f:
        c_all = f['field'][start:start + tpc, :, :]
    with h5py.File(os.path.join(DATABASE_DIR, 'channel_2.h5'), 'r') as f:
        u_all = f['field'][start:start + tpc, :, :]
    with h5py.File(os.path.join(DATABASE_DIR, 'channel_3.h5'), 'r') as f:
        v_all = f['field'][start:start + tpc, :, :]

    X = np.zeros((tpc, 3, 256, 512), dtype=np.float32)

    # G channel: per-snapshot to [0, 1]
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

    # U, V channels: per-case fluctuation to [0, 1]
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
# Build normalized ground truth
# ============================================================
print("\n[INFO] Loading & normalizing ground truth...")

# Cache normalized cases
case_cache = {}
for case_id in unique_cases:
    print(f"  Normalizing case {case_id}...")
    case_cache[case_id] = normalize_case_mpf(case_id)

# Assemble selected snapshots
U_true = np.zeros((N_test, 3, 256, 512), dtype=np.float32)
for i, (case_id, t_idx) in enumerate(case_time_pairs):
    U_true[i] = case_cache[case_id][t_idx]

print(f"[INFO] U_true: {U_true.shape}")

# ============================================================
# Build MLP input: [DR, VR, Re, Bo, t]
# ============================================================
X_test = np.zeros((N_test, n_params + 1), dtype=np.float32)
for i, (case_id, t_idx) in enumerate(case_time_pairs):
    ci = test_cases.index(case_id)
    X_test[i, :n_params] = test_params[ci]
    X_test[i, n_params]  = timestep_values[t_idx]

# Normalize with training z-score stats
norm = np.load(NORM_PATH)
X_test_norm = (X_test - norm['input_mean']) / norm['input_std']
print(f"[INFO] X_test: {X_test.shape}")

# ============================================================
# Load MLP
# ============================================================
mlp_info = json.load(open(os.path.join(MLP_RUN_DIR, 'run_info.json')))
latent_dim_actual = mlp_info['latent_dim']
input_dim         = mlp_info['input_dim']

mlp = LatentMLP(input_dim=input_dim, hidden_dim=HIDDEN_DIM,
                latent_dim=latent_dim_actual).to(device)
mlp.load_state_dict(torch.load(MLP_MODEL_PATH, map_location=device,
                                weights_only=True))
mlp.eval()
print(f"[INFO] MLP loaded: {input_dim} → {latent_dim_actual}")

# ============================================================
# Load Decoder
# ============================================================
decoder = ConvAutoencoderBaseline(
    dim=DIM, dim_mults=DIM_MULTS, channels=ENCODER_CHANNELS,
    z_channels=Z_CHANNELS, block_type=BLOCK_TYPE,
).float().to(device)

sd = torch.load(DECODER_MODEL_PATH, map_location=device, weights_only=True)
if list(sd.keys())[0].startswith("module."):
    sd = {k.replace("module.", "", 1): v for k, v in sd.items()}
decoder.load_state_dict(sd, strict=True)
decoder.eval()
print(f"[INFO] Decoder loaded: {KTH_RUN_NAME}")

# Infer latent spatial shape
with torch.no_grad():
    dummy = torch.zeros(1, ENCODER_CHANNELS, 256, 512,
                        dtype=torch.float32, device=device)
    z_shape = decoder.encode(dummy).shape[1:]   # (4, 4, 8)
print(f"[INFO] Latent shape: {z_shape}  (flat={np.prod(z_shape)})")

# ============================================================
# Predict: (mu, t) → MLP → z → Decoder → u_pred
# ============================================================
print("\n[INFO] Predicting...")

X_test_t = torch.tensor(X_test_norm, dtype=torch.float32)
U_pred_list = []

with torch.no_grad():
    for i in range(0, N_test, PRED_BATCH):
        x_batch = X_test_t[i:i+PRED_BATCH].to(device)
        z_flat  = mlp(x_batch)
        z       = z_flat.reshape(-1, *z_shape)
        u_batch = decoder.decode(z)
        U_pred_list.append(u_batch.cpu().numpy())

U_pred = np.concatenate(U_pred_list, axis=0)
print(f"[INFO] U_pred: {U_pred.shape}")

# Free decoder memory
del decoder
if device.type == 'cuda':
    torch.cuda.empty_cache()

# ============================================================
# Error Metrics (in normalized space)
# ============================================================
err = U_true - U_pred
mae_per   = np.mean(np.abs(err), axis=(1, 2, 3))
mse_per   = np.mean(err**2,      axis=(1, 2, 3))
linf_per  = np.max(np.abs(err),  axis=(1, 2, 3))

# Per-channel errors (averaged across samples)
channel_names = ['G (gas)', 'U (vel-x)', 'V (vel-y)']
print(f"\n{'='*60}")
print("Per-sample metrics (normalized space)")
print(f"{'='*60}")
for i, (case_id, t_idx) in enumerate(case_time_pairs):
    print(f"  [{i:2d}] case={case_id} t={t_idx:3d}  "
          f"MSE={mse_per[i]:.6f}  MAE={mae_per[i]:.6f}  Linf={linf_per[i]:.6f}")

print(f"\n{'='*60}")
print("Per-channel metrics (all samples)")
print(f"{'='*60}")
for ch in range(3):
    ch_err = U_true[:, ch] - U_pred[:, ch]
    ch_mse = np.mean(ch_err**2)
    ch_mae = np.mean(np.abs(ch_err))
    print(f"  {channel_names[ch]:12s}: MSE={ch_mse:.6f}, MAE={ch_mae:.6f}")

print(f"\n  Overall — MAE: {mae_per.mean():.6f}, MSE: {mse_per.mean():.6f}, "
      f"Linf: {linf_per.mean():.6f}")

# ============================================================
# Save results
# ============================================================
# Build metadata arrays
case_ids_arr  = np.array([c for c, _ in case_time_pairs], dtype=np.int64)
timesteps_arr = np.array([t for _, t in case_time_pairs], dtype=np.int64)
params_arr    = np.array([test_params[test_cases.index(c)]
                          for c, _ in case_time_pairs], dtype=np.float32)

np.savez(os.path.join(OUTPUT_DIR, 'pipeline_predictions.npz'),
         originals=U_true,
         predictions=U_pred,
         X_test=X_test,
         case_ids=case_ids_arr,
         timesteps=timesteps_arr,
         params=params_arr,
         channel_names=np.array(channel_names),
         mae=mae_per, mse=mse_per, linf=linf_per)

print(f"\n[INFO] Saved: {OUTPUT_DIR}/pipeline_predictions.npz")
print(f"[INFO] Done!")
