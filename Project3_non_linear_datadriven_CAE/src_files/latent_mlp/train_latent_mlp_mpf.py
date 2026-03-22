# Latent MLP Training Script — MPF Dataset (Stage 2)
# ============================================================
# Dedicated script for the MPF multiphase-flow dataset.
#
# The MPF encoder was trained with a NON-TRIVIAL normalization
# pipeline (per-snapshot G normalization, per-case U/V fluctuation
# normalization) that differs from the standard raw-HDF5 path used
# by the Advection scripts.  This script replicates that exact
# normalization so latent codes match the trained encoder.
#
# Encoder: KTH ResNet CAE (3-channel: G, U, V)
# Run:     2DMPF_GUV_Nt2424_B16_NL8_LR1e-05_L21e-06_E300
#
# Usage:  python train_latent_mlp_mpf.py   (run from latent_mlp/)
# ============================================================

# ////////////////////////////////////////////////////////////
# USER CONFIGURATION — modify only this section
# ////////////////////////////////////////////////////////////

# --- Dataset ---
DATASET_FOLDER = "MPF"

# --- Encoder (KTH ResNet CAE, 3-channel) ---
DIM            = 16
DIM_MULTS      = (1, 2, 4, 8, 16, 32, 64, 128)   # 8 levels
Z_CHANNELS     = 4
BLOCK_TYPE     = 1
ENCODER_CHANNELS = 3                                # G, U, V
KTH_RUN_NAME   = "2DMPF_GUV_Nt2424_B16_NL8_LR1e-05_L21e-06_E300"

# --- MLP ---
HIDDEN_DIM     = 128
LEARNING_RATE  = 1e-3
BATCH_SIZE     = 64
MAX_EPOCHS     = 5000
PATIENCE       = 500

# --- Reproducibility ---
SEED           = 42

# ////////////////////////////////////////////////////////////
# END OF USER CONFIGURATION
# ////////////////////////////////////////////////////////////

# ============================================================
# Libraries
# ============================================================
import sys, os
import json
import time
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import h5py

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# Load Modules
# ============================================================
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

sys.path.append(os.path.join(SCRIPT_DIR, 'library_latent_mlp'))
from latent_mlp.models.latent_mlp import LatentMLP
from latent_mlp.utils import count_parameters

sys.path.append(os.path.join(PROJECT_ROOT, 'resnet_cae', 'library_resnet_cae'))
from resnet_cae.models.baseline_model import ConvAutoencoderBaseline   # type: ignore

print('[INFO] Modules loaded\n-----------------')

# ============================================================
# Paths
# ============================================================
DATABASE_DIR       = os.path.join(PROJECT_ROOT, 'database', DATASET_FOLDER)
SPLIT_FILE         = os.path.join(DATABASE_DIR, 'split.json')
ENCODER_MODEL_PATH = os.path.join(PROJECT_ROOT, 'resnet_cae', 'runs',
                                  KTH_RUN_NAME, 'model_final.pt')

# ============================================================
# Random Seeds
# ============================================================
np.random.seed(SEED)
torch.manual_seed(SEED)

# ============================================================
# Device Setup
# ============================================================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[INFO] Device: {device}")
if device.type == 'cuda':
    print(f"  GPU: {torch.cuda.get_device_name(0)}")

# ============================================================
# Load Split
# ============================================================
print("\n" + "=" * 60)
print("Loading Split & Metadata")
print("=" * 60)

with open(SPLIT_FILE, 'r') as f:
    split = json.load(f)

timesteps_per_case = split['timesteps_per_case']
timestep_values    = np.array(split['timestep_values'], dtype=np.float32)
param_names        = split['param_names']
n_params           = len(param_names)

train_cases  = split['train_cases']
val_cases    = split['val_cases']
train_params = split['train_params']
val_params   = split['val_params']

n_train = len(train_cases) * timesteps_per_case
n_val   = len(val_cases)   * timesteps_per_case

print(f"  Dataset:          {split['dataset']}")
print(f"  Encoder run:      {KTH_RUN_NAME}")
print(f"  Param names:      {param_names}  ({n_params}D)")
print(f"  Timesteps/case:   {timesteps_per_case}")
print(f"  Train:            {len(train_cases)} cases = {n_train} snapshots")
print(f"  Val:              {len(val_cases)} cases = {n_val} snapshots")

# ============================================================
# MPF Normalization
# ============================================================
# Replicates the exact normalization from the original training
# script (AE_2D_edited_multiple_channel.py) so that the latent
# codes match what the encoder was trained on.
# ============================================================

def normalize_case_mpf(case_idx, database_dir, tpc):
    """
    Load one MPF case from HDF5 and apply the original normalization.

    Returns: np.ndarray of shape (tpc, 3, 256, 512) float32 in [0, 1]
    """
    start = case_idx * tpc

    with h5py.File(os.path.join(database_dir, 'channel_1.h5'), 'r') as f:
        c_all = f['field'][start:start + tpc, :, :]       # (T, 512, 256)
    with h5py.File(os.path.join(database_dir, 'channel_2.h5'), 'r') as f:
        u_all = f['field'][start:start + tpc, :, :]
    with h5py.File(os.path.join(database_dir, 'channel_3.h5'), 'r') as f:
        v_all = f['field'][start:start + tpc, :, :]

    X = np.zeros((tpc, 3, 256, 512), dtype=np.float32)

    # ---- G channel: per-snapshot to [0, 1] ----
    for t in range(tpc):
        G = c_all[t].T.astype(np.float32)                 # (512,256) → (256,512)
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
    U_flat = u_all.reshape(tpc, ns).T      # (ns, T)
    V_flat = v_all.reshape(tpc, ns).T

    UP = U_flat - np.mean(U_flat, axis=1, keepdims=True)
    VP = V_flat - np.mean(V_flat, axis=1, keepdims=True)

    UNor = (UP - UP.min()) / (UP.max() - UP.min() + 1e-8)
    VNor = (VP - VP.min()) / (VP.max() - VP.min() + 1e-8)

    for t in range(tpc):
        X[t, 1] = UNor[:, t].reshape(512, 256).T          # (ns,)→(512,256)→(256,512)
        X[t, 2] = VNor[:, t].reshape(512, 256).T

    return X

# ============================================================
# Load Encoder
# ============================================================
print("\n" + "=" * 60)
print("Loading KTH Encoder")
print("=" * 60)

encoder = ConvAutoencoderBaseline(
    dim=DIM, dim_mults=DIM_MULTS, channels=ENCODER_CHANNELS,
    z_channels=Z_CHANNELS, block_type=BLOCK_TYPE,
).float().to(device)

sd = torch.load(ENCODER_MODEL_PATH, map_location=device, weights_only=True)
if list(sd.keys())[0].startswith("module."):
    sd = {k.replace("module.", "", 1): v for k, v in sd.items()}
encoder.load_state_dict(sd, strict=True)
encoder.eval()
print(f"  Loaded: {ENCODER_MODEL_PATH}")

# ============================================================
# Compute Latent Codes (case by case)
# ============================================================
print("\n" + "=" * 60)
print("Computing Latent Codes (per-case normalization → encode)")
print("=" * 60)

ENCODE_BATCH = 4   # small batches — 450M param model


def encode_cases(case_list, label):
    """
    For each case: load + normalize + encode → flat latent vectors.
    Returns: Z of shape (n_cases * tpc, latent_dim_flat)
    """
    all_latents = []
    for ci, case_idx in enumerate(case_list):
        X = normalize_case_mpf(case_idx, DATABASE_DIR, timesteps_per_case)
        case_latents = []
        with torch.no_grad():
            for i in range(0, timesteps_per_case, ENCODE_BATCH):
                batch = torch.from_numpy(X[i:i + ENCODE_BATCH]).to(device)
                z = encoder.encode(batch)
                z_flat = z.reshape(z.size(0), -1)
                case_latents.append(z_flat.cpu().numpy())
        all_latents.append(np.concatenate(case_latents, axis=0))
        if (ci + 1) % 5 == 0 or ci == len(case_list) - 1:
            print(f"  [{label}] {ci + 1}/{len(case_list)} cases encoded")
    return np.concatenate(all_latents, axis=0)


Z_train = encode_cases(train_cases, "Train")
Z_val   = encode_cases(val_cases,   "Val")
latent_dim_actual = Z_train.shape[1]

print(f"  Z_train: {Z_train.shape}")
print(f"  Z_val:   {Z_val.shape}")
print(f"  Latent dim: {latent_dim_actual}")

# Free encoder memory
del encoder
if device.type == 'cuda':
    torch.cuda.empty_cache()

# ============================================================
# Build MLP Input Arrays: [DR, VR, Re, Bo, t]
# ============================================================
print("\n" + "=" * 60)
print("Building MLP Inputs")
print("=" * 60)


def build_inputs(case_list, params_list):
    """Build (n_cases * tpc, n_params + 1) array with [params..., t]."""
    n = len(case_list) * timesteps_per_case
    X = np.zeros((n, n_params + 1), dtype=np.float32)
    row = 0
    for ci, case_idx in enumerate(case_list):
        p = params_list[ci]
        for t_idx in range(timesteps_per_case):
            X[row, :n_params] = p
            X[row, n_params]  = timestep_values[t_idx]
            row += 1
    return X


X_train = build_inputs(train_cases, train_params)
X_val   = build_inputs(val_cases,   val_params)
print(f"  X_train: {X_train.shape}  (columns: {param_names + ['t']})")
print(f"  X_val:   {X_val.shape}")

# ============================================================
# Normalize MLP Inputs (z-score)
# ============================================================
input_mean = X_train.mean(axis=0)
input_std  = X_train.std(axis=0)
input_std[input_std < 1e-12] = 1.0

X_train_norm = (X_train - input_mean) / input_std
X_val_norm   = (X_val   - input_mean) / input_std

# ============================================================
# Results Directory
# ============================================================
run_name = (f"{split['dataset']}_KTH_MLP_HD{HIDDEN_DIM}_LD{latent_dim_actual}"
            f"_LR{LEARNING_RATE}_E{MAX_EPOCHS}")
results_directory = os.path.join('runs', run_name)
os.makedirs(results_directory, exist_ok=True)
os.makedirs(os.path.join(results_directory, 'metrics'), exist_ok=True)
print(f"\n[INFO] Run: {run_name}")

# Save normalization
np.savez(os.path.join(results_directory, 'normalization_params.npz'),
         input_mean=input_mean, input_std=input_std)

with open(os.path.join(results_directory, 'normalization_info.txt'), 'w') as f:
    col_names = param_names + ['t']
    f.write("Input normalization: z-score per column\n")
    for i, name in enumerate(col_names):
        f.write(f"  {name}: mean={input_mean[i]:.6f}, std={input_std[i]:.6f}\n")
    f.write("\nField normalization (applied before encoding):\n")
    f.write("  G (ch1): per-snapshot min-max to [0, 1]\n")
    f.write("  U (ch2): per-case temporal-mean subtracted, per-case min-max to [0, 1]\n")
    f.write("  V (ch3): per-case temporal-mean subtracted, per-case min-max to [0, 1]\n")
    f.write("  Spatial: HDF5 (512,256) transposed to (256,512) for model input\n")

# ============================================================
# PyTorch DataLoaders
# ============================================================
X_train_t = torch.tensor(X_train_norm, dtype=torch.float32)
Y_train_t = torch.tensor(Z_train,      dtype=torch.float32)
X_val_t   = torch.tensor(X_val_norm,   dtype=torch.float32)
Y_val_t   = torch.tensor(Z_val,        dtype=torch.float32)

train_loader = DataLoader(TensorDataset(X_train_t, Y_train_t),
                          batch_size=BATCH_SIZE, shuffle=True)
val_loader   = DataLoader(TensorDataset(X_val_t, Y_val_t),
                          batch_size=BATCH_SIZE, shuffle=False)

# ============================================================
# MLP Model
# ============================================================
print("\n" + "=" * 60)
print("Creating MLP Model")
print("=" * 60)

input_dim = n_params + 1   # (DR, VR, Re, Bo, t)
mlp = LatentMLP(input_dim=input_dim, hidden_dim=HIDDEN_DIM,
                latent_dim=latent_dim_actual).to(device)

total_params_mlp, trainable_params_mlp = count_parameters(mlp)
print(f"  Architecture: {input_dim} -> {HIDDEN_DIM} -> {HIDDEN_DIM} -> {latent_dim_actual}")
print(f"  Parameters: {trainable_params_mlp:,}")

# ============================================================
# run_info.json
# ============================================================
run_info = {
    'dataset': split['dataset'],
    'encoder_type': 'KTH',
    'encoder_run': KTH_RUN_NAME,
    'encoder_model_path': ENCODER_MODEL_PATH,
    'encoder_channels': ENCODER_CHANNELS,
    'encoder_dim_mults': list(DIM_MULTS),
    'latent_dim': latent_dim_actual,
    'input_dim': input_dim,
    'hidden_dim': HIDDEN_DIM,
    'learning_rate': LEARNING_RATE,
    'batch_size': BATCH_SIZE,
    'max_epochs': MAX_EPOCHS,
    'patience': PATIENCE,
    'seed': SEED,
    'total_parameters': total_params_mlp,
    'trainable_parameters': trainable_params_mlp,
    'train_cases': train_cases,
    'val_cases': val_cases,
    'n_train_snapshots': n_train,
    'n_val_snapshots': n_val,
    'split_file': SPLIT_FILE,
    'field_normalization': 'MPF custom (per-snapshot G, per-case UV fluctuation)',
}
with open(os.path.join(results_directory, 'run_info.json'), 'w') as f:
    json.dump(run_info, f, indent=2)
print(f"[INFO] Saved run_info.json")

# ============================================================
# GPU Memory Summary
# ============================================================
if device.type == 'cuda':
    with torch.no_grad():
        _ = mlp(torch.randn(1, input_dim).to(device))
    with open(os.path.join(results_directory, 'memory_summary.txt'), 'w') as f:
        f.write(torch.cuda.memory_summary())
    print(f"[INFO] Saved memory_summary.txt")

# ============================================================
# Training
# ============================================================
criterion = nn.MSELoss()
optimizer = optim.Adam(mlp.parameters(), lr=LEARNING_RATE)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                  factor=0.5, patience=100)

print("\n" + "=" * 60)
print("Training MLP")
print("=" * 60)

train_losses = []
val_losses   = []
best_val_loss = float('inf')
best_epoch    = 0
epochs_no_improve = 0
start_time = time.time()

for epoch in range(MAX_EPOCHS):
    # --- Train ---
    mlp.train()
    epoch_train_loss = 0.0
    for X_batch, Y_batch in train_loader:
        X_batch = X_batch.to(device)
        Y_batch = Y_batch.to(device)
        pred = mlp(X_batch)
        loss = criterion(pred, Y_batch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch_train_loss += loss.item() * X_batch.size(0)
    epoch_train_loss /= len(train_loader.dataset)
    train_losses.append(epoch_train_loss)

    # --- Validate ---
    mlp.eval()
    epoch_val_loss = 0.0
    with torch.no_grad():
        for X_batch, Y_batch in val_loader:
            X_batch = X_batch.to(device)
            Y_batch = Y_batch.to(device)
            pred = mlp(X_batch)
            loss = criterion(pred, Y_batch)
            epoch_val_loss += loss.item() * X_batch.size(0)
    epoch_val_loss /= len(val_loader.dataset)
    val_losses.append(epoch_val_loss)

    scheduler.step(epoch_val_loss)

    # --- Early stopping ---
    if epoch_val_loss < best_val_loss:
        best_val_loss = epoch_val_loss
        best_epoch = epoch
        epochs_no_improve = 0
        torch.save({
            'epoch': epoch + 1,
            'model': mlp.state_dict(),
            'opt': optimizer.state_dict(),
        }, os.path.join(results_directory, 'model-best.pt'))
    else:
        epochs_no_improve += 1

    if epoch % 100 == 0 or epoch == MAX_EPOCHS - 1:
        lr_now = optimizer.param_groups[0]['lr']
        print(f"  Epoch {epoch:5d} | Train: {epoch_train_loss:.6f} | "
              f"Val: {epoch_val_loss:.6f} | LR: {lr_now:.6f}")

    if epochs_no_improve >= PATIENCE:
        print(f"\n  Early stopping at epoch {epoch} (best: epoch {best_epoch})")
        break

total_time = time.time() - start_time
print(f"\n  Training complete! Time: {total_time:.1f}s")
print(f"  Best val loss: {best_val_loss:.6f} at epoch {best_epoch}")

# ============================================================
# Save Final Model (best weights)
# ============================================================
best_ckpt = torch.load(os.path.join(results_directory, 'model-best.pt'),
                        map_location=device)
mlp.load_state_dict(best_ckpt['model'])

torch.save(mlp.state_dict(), os.path.join(results_directory, 'model_final.pt'))
print(f"[INFO] Saved model_final.pt (from best epoch {best_epoch})")

# ============================================================
# Loss History
# ============================================================
loss_history = {'Train': train_losses, 'Val': val_losses}
with open(os.path.join(results_directory, 'loss_history.json'), 'w') as f:
    json.dump(loss_history, f)

fig, ax = plt.subplots(figsize=(3, 3), dpi=200)
ax.plot(train_losses, label='Train')
ax.plot(val_losses, label='Val')
ax.axvline(best_epoch, color='red', linestyle='--', alpha=0.5,
           label=f'Best (ep {best_epoch})')
ax.set_yscale('log')
ax.set_xlabel('Epoch')
ax.set_ylabel('Loss')
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(results_directory, 'loss_history.png'))
plt.close(fig)
print(f"[INFO] Saved loss_history.json and loss_history.png")

# ============================================================
# Evaluation: Latent Code Comparison (Val Set)
# ============================================================
print("\n" + "=" * 60)
print("Evaluation — Latent Comparison (Val Set)")
print("=" * 60)

mlp.eval()
with torch.no_grad():
    Z_pred_val = mlp(X_val_t.to(device)).cpu().numpy()

n_plot = min(latent_dim_actual, 4)
fig, axes = plt.subplots(1, n_plot, figsize=(4 * n_plot, 4))
if n_plot == 1:
    axes = [axes]
for i, ax in enumerate(axes):
    ax.scatter(Z_val[:, i], Z_pred_val[:, i], alpha=0.3, s=10)
    z_min, z_max = Z_val[:, i].min(), Z_val[:, i].max()
    ax.plot([z_min, z_max], [z_min, z_max], 'r--', label='Perfect')
    ax.set_xlabel(f'True z[{i}]')
    ax.set_ylabel(f'Predicted z[{i}]')
    ax.set_title(f'Latent dim {i}')
    ax.legend()
    ax.grid(True, alpha=0.3)
plt.suptitle('KTH MPF Latent Code Prediction (Val Set)', fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(results_directory, 'latent_comparison.png'), dpi=150)
plt.close(fig)
print(f"[INFO] Saved latent_comparison.png")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
print(f"[INFO] All results saved to: {results_directory}")
