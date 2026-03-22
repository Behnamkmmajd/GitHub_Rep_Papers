# Latent MLP Training Script — Stage 2
# ============================================================
# Trains an MLP to map physical parameters (mu, t) → latent
# codes z produced by a trained Stage-1 encoder (CNN or KTH).
#
# Reads HDF5 field data + split.json (same as the CAE training
# scripts), loads the trained encoder, encodes every training
# and validation snapshot to get target latent codes, then
# trains the MLP to predict those codes from (mu, t).
#
# Usage:  python train_latent_mlp.py   (run from latent_mlp/)
# ////////////////////////////////////////////////////////////
# USER CONFIGURATION — modify only this section
# ////////////////////////////////////////////////////////////

# --- Dataset ---
DATASET_FOLDER = "Advection"            # Subfolder inside database/
CHANNELS       = [1]                    # Which channel files to load

# --- Encoder selection ---
ENCODER_TYPE   = "KTH"                  # "CNN" or "KTH"

# CNN encoder settings (used when ENCODER_TYPE == "CNN")
LATENT_DIM     = 4
CNN_RUN_NAME   = "Advection_Ch1_Nt240_B32_LD4_LR0.0001_E5000_Norm"

# KTH encoder settings (used when ENCODER_TYPE == "KTH")
DIM            = 16
DIM_MULTS      = (1, 2, 4, 8, 16, 32)          # Must match the trained model
Z_CHANNELS     = 4
BLOCK_TYPE     = 1
KTH_RUN_NAME   = "Advection_Ch1_Nt240_B8_NL6_LR1e-05_L21e-06_E300_Norm"

# --- MLP ---
HIDDEN_DIM     = 64
LEARNING_RATE  = 1e-3
BATCH_SIZE     = 32
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
# Load MLP Modules
# ============================================================
sys.path.append(os.path.join(os.path.dirname(__file__), 'library_latent_mlp'))

from latent_mlp.models.latent_mlp import LatentMLP
from latent_mlp.utils import count_parameters

print('[INFO] Modules loaded\n-----------------')

# ============================================================
# Load Encoder from Source Project
# ============================================================
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

if ENCODER_TYPE == "CNN":
    sys.path.append(os.path.join(PROJECT_ROOT, 'simple_cae',
                                 'library_simple_cae'))
    from simple_cae.models.cae_model import ConvolutionalAutoencoder       # type: ignore
    ENCODER_MODEL_PATH = os.path.join(PROJECT_ROOT, 'simple_cae',
                                      'runs', CNN_RUN_NAME, 'model_final.pt')
elif ENCODER_TYPE == "KTH":
    sys.path.append(os.path.join(PROJECT_ROOT, 'resnet_cae', 'library_resnet_cae'))
    from resnet_cae.models.baseline_model import ConvAutoencoderBaseline   # type: ignore
    ENCODER_MODEL_PATH = os.path.join(PROJECT_ROOT, 'resnet_cae',
                                      'runs', KTH_RUN_NAME, 'model_final.pt')
else:
    raise ValueError(f"Unknown ENCODER_TYPE: {ENCODER_TYPE!r}. Use 'CNN' or 'KTH'.")

# ============================================================
# Paths
# ============================================================
DATABASE_DIR = os.path.join(PROJECT_ROOT, 'database', DATASET_FOLDER)
SPLIT_FILE   = os.path.join(DATABASE_DIR, 'split.json')

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
# Load Split & Data
# ============================================================
print("\n" + "=" * 60)
print("Loading Data")
print("=" * 60)

with open(SPLIT_FILE, 'r') as f:
    split = json.load(f)

n_channels         = len(CHANNELS)
timesteps_per_case = split['timesteps_per_case']
timestep_values    = np.array(split['timestep_values'])
param_names        = split['param_names']
n_params           = len(param_names)
spatial_shape      = split['spatial_shape']

train_indices = split['train_indices']
val_indices   = split['val_indices']
train_cases   = split['train_cases']
val_cases     = split['val_cases']
train_params  = split['train_params']    # list of [mu, ...] per case
val_params    = split['val_params']

print(f"  Dataset:          {split['dataset']}")
print(f"  Channels:         {CHANNELS}")
print(f"  Encoder:          {ENCODER_TYPE}")
print(f"  Param names:      {param_names}")
print(f"  Timesteps/case:   {timesteps_per_case}")
print(f"  Timestep values:  {timestep_values}")
print(f"  Train snapshots:  {len(train_indices)}")
print(f"  Val snapshots:    {len(val_indices)}")

# --- Load HDF5 fields ---
def load_fields(indices):
    """Load multi-channel field data from HDF5 for given row indices."""
    arrays = []
    for ch in CHANNELS:
        h5_path = os.path.join(DATABASE_DIR, f'channel_{ch}.h5')
        with h5py.File(h5_path, 'r') as f:
            arrays.append(f['field'][indices])           # (N, H, W)
    return np.stack(arrays, axis=1).astype(np.float32)   # (N, C, H, W)

U_train = load_fields(train_indices)
U_val   = load_fields(val_indices)
print(f"  U_train: {U_train.shape}")
print(f"  U_val:   {U_val.shape}")

# --- Check if encoder was trained with normalization ---
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

if ENCODER_TYPE == "CNN":
    encoder_run_dir = os.path.join(PROJECT_ROOT, 'simple_cae', 'runs', CNN_RUN_NAME)
else:
    encoder_run_dir = os.path.join(PROJECT_ROOT, 'resnet_cae', 'runs', KTH_RUN_NAME)

encoder_ri_path = os.path.join(encoder_run_dir, 'run_info.json')
field_normalize = False
if os.path.exists(encoder_ri_path):
    with open(encoder_ri_path) as f:
        field_normalize = json.load(f).get('normalize', False)

if field_normalize:
    print("[INFO] Encoder trained with normalization — applying per-snapshot min-max [0,1]")
    U_train = normalize_per_snapshot(U_train)
    U_val   = normalize_per_snapshot(U_val)

# --- Build (params..., t) input arrays ---
def build_inputs(indices, case_list, params_list):
    """Build MLP input array [mu..., t] for each snapshot index."""
    n = len(indices)
    X = np.zeros((n, n_params + 1), dtype=np.float32)
    for i, idx in enumerate(indices):
        case_id   = idx // timesteps_per_case
        t_idx     = idx % timesteps_per_case
        local_pos = case_list.index(case_id)
        X[i, :n_params] = params_list[local_pos]
        X[i, n_params]  = timestep_values[t_idx]
    return X

X_train = build_inputs(train_indices, train_cases, train_params)
X_val   = build_inputs(val_indices,   val_cases,   val_params)
print(f"  X_train: {X_train.shape}  (columns: {param_names + ['t']})")

# ============================================================
# Load Encoder & Compute Latent Codes
# ============================================================
print("\n" + "=" * 60)
print(f"Loading {ENCODER_TYPE} Encoder & Computing Latent Codes")
print("=" * 60)

if ENCODER_TYPE == "CNN":
    encoder = ConvolutionalAutoencoder(latent_dim=LATENT_DIM,
                                      in_channels=n_channels).to(device)
    sd = torch.load(ENCODER_MODEL_PATH, map_location=device, weights_only=True)
    encoder.load_state_dict(sd)
else:
    encoder = ConvAutoencoderBaseline(
        dim=DIM, dim_mults=DIM_MULTS, channels=n_channels,
        z_channels=Z_CHANNELS, block_type=BLOCK_TYPE,
    ).float().to(device)
    sd = torch.load(ENCODER_MODEL_PATH, map_location=device, weights_only=True)
    if list(sd.keys())[0].startswith("module."):
        sd = {k.replace("module.", "", 1): v for k, v in sd.items()}
    encoder.load_state_dict(sd, strict=True)

encoder.eval()
print(f"  Loaded: {ENCODER_MODEL_PATH}")


def encode_all(U, batch_size=64):
    """Encode field data through the frozen encoder -> flat latent vectors."""
    latent_parts = []
    with torch.no_grad():
        for i in range(0, len(U), batch_size):
            batch = torch.tensor(U[i:i+batch_size], dtype=torch.float32).to(device)
            z = encoder.encode(batch)
            z_flat = z.reshape(z.size(0), -1)   # flatten spatial latent (KTH) or no-op (CNN)
            latent_parts.append(z_flat.cpu().numpy())
    return np.concatenate(latent_parts, axis=0)


Z_train = encode_all(U_train)
Z_val   = encode_all(U_val)
latent_dim_actual = Z_train.shape[1]

print(f"  Z_train: {Z_train.shape}")
print(f"  Z_val:   {Z_val.shape}")
print(f"  Latent dim (actual): {latent_dim_actual}")

# Free encoder memory
del encoder, U_train, U_val
if device.type == 'cuda':
    torch.cuda.empty_cache()

# ============================================================
# Normalize MLP Inputs (z-score)
# ============================================================
input_mean = X_train.mean(axis=0)
input_std  = X_train.std(axis=0)
input_std[input_std < 1e-12] = 1.0     # avoid division by zero for constant columns

X_train_norm = (X_train - input_mean) / input_std
X_val_norm   = (X_val   - input_mean) / input_std

# ============================================================
# Results Directory
# ============================================================
ch_str = "".join(str(c) for c in CHANNELS)
run_name = (f"{split['dataset']}_Ch{ch_str}_{ENCODER_TYPE}"
            f"_MLP_HD{HIDDEN_DIM}_LD{latent_dim_actual}"
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

input_dim = n_params + 1   # (mu..., t)
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
    'channels': CHANNELS,
    'encoder_type': ENCODER_TYPE,
    'encoder_model_path': ENCODER_MODEL_PATH,
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
    'n_train_snapshots': len(train_indices),
    'n_val_snapshots': len(val_indices),
    'split_file': SPLIT_FILE,
    'field_normalize': field_normalize,
    'field_normalization_type': 'per_snapshot_minmax_01' if field_normalize else 'none',
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
ax.axvline(best_epoch, color='red', linestyle='--', alpha=0.5, label=f'Best (ep {best_epoch})')
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
plt.suptitle(f'{ENCODER_TYPE} Latent Code Prediction (Val Set)', fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(results_directory, 'latent_comparison.png'), dpi=150)
plt.close(fig)
print(f"[INFO] Saved latent_comparison.png")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
print(f"[INFO] All results saved to: {results_directory}")
