# CAE Training Script — CNN-ROM (Wang et al. 2024)
# ============================================================
# Trains a Convolutional Autoencoder on HDF5 field data using
# a standardized split.json for train/val/test partitioning.
# Works with any dataset (Advection, MPF, etc.) and any number
# of channels.
#
# Usage:  python train_cae.py   (run from simple_cae/)
# ============================================================

# ////////////////////////////////////////////////////////////
# USER CONFIGURATION — modify only this section
# ////////////////////////////////////////////////////////////

# --- Dataset ---
DATASET_FOLDER = "Advection"            # Subfolder inside database/
CHANNELS       = [1]                    # Which channel files to load
                                        #   [1]        → channel_1.h5 (single)
                                        #   [1, 2, 3]  → channels 1–3 (multi)

# --- Architecture ---
LATENT_DIM     = 4                      # Latent vector size

# --- Normalization ---
NORMALIZE      = True                   # Per-snapshot min-max to [0, 1]

# --- Training ---
BATCH_SIZE     = 32
LEARNING_RATE  = 1e-4
NUM_EPOCHS     = 5000
PATIENCE       = 500                    # Early stopping patience (epochs)

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
import matplotlib.pyplot as plt
import pandas as pd
import h5py

# ============================================================
# Load CNN-ROM Modules
# ============================================================
sys.path.append(os.path.join(os.path.dirname(__file__), 'library_simple_cae'))

from simple_cae.models.cae_model import ConvolutionalAutoencoder
from simple_cae.metrics import MetricType, compute_metrics_single_array
from simple_cae.utils import count_parameters

print('[INFO] Modules loaded\n-----------------')

# ============================================================
# Paths
# ============================================================
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATABASE_DIR = os.path.join(PROJECT_ROOT, 'database', DATASET_FOLDER)
SPLIT_FILE   = os.path.join(DATABASE_DIR, 'split.json')

# ============================================================
# Read split.json
# ============================================================
with open(SPLIT_FILE, 'r') as f:
    split = json.load(f)

n_channels         = len(CHANNELS)
timesteps_per_case = split['timesteps_per_case']
spatial_shape      = tuple(split['spatial_shape'])       # (H, W)
train_indices      = split['train_indices']
val_indices        = split['val_indices']
# test_indices are NOT loaded — reserved for final evaluation only

print(f"[INFO] Dataset:    {split['dataset']}")
print(f"[INFO] Channels:   {CHANNELS} ({n_channels} ch)")
print(f"[INFO] Spatial:    {spatial_shape}")
print(f"[INFO] Train:      {len(train_indices)} snapshots  "
      f"({len(split['train_cases'])} cases)")
print(f"[INFO] Val:        {len(val_indices)} snapshots  "
      f"({len(split['val_cases'])} cases)")
print(f"[INFO] Test:       {len(split['test_indices'])} snapshots  "
      f"({len(split['test_cases'])} cases)  [NOT loaded]")

# ============================================================
# Results Directory
# ============================================================
ch_tag = "".join(str(c) for c in CHANNELS)
norm_tag = "_Norm" if NORMALIZE else ""
run_name = (f"{split['dataset']}_Ch{ch_tag}_"
            f"Nt{len(train_indices)}_B{BATCH_SIZE}_"
            f"LD{LATENT_DIM}_LR{LEARNING_RATE}_E{NUM_EPOCHS}{norm_tag}")
results_directory = os.path.join('runs', run_name)
os.makedirs(results_directory, exist_ok=True)
os.makedirs(os.path.join(results_directory, 'metrics'), exist_ok=True)
os.makedirs(os.path.join(results_directory, 'intensity_histograms'), exist_ok=True)
print(f"[INFO] Results:    {results_directory}")

# ============================================================
# Random Seeds & Device
# ============================================================
np.random.seed(SEED)
torch.manual_seed(SEED)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[INFO] Device:     {device}")
if device.type == 'cuda':
    print(f"  GPU: {torch.cuda.get_device_name(0)}")
    print(f"  Memory: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")

# ============================================================
# Load Data from HDF5
# ============================================================
print("\n" + "=" * 60)
print("Loading Data")
print("=" * 60)


def load_snapshots_from_h5(database_dir, channels, indices):
    """
    Load selected snapshot rows from one or more channel HDF5 files.

    Returns
    -------
    np.ndarray of shape (N, C, H, W), float32
    """
    arrays = []
    for ch in channels:
        h5_path = os.path.join(database_dir, f'channel_{ch}.h5')
        with h5py.File(h5_path, 'r') as f:
            # Load only the requested rows (efficient HDF5 fancy indexing)
            field = f['field'][sorted(indices)]             # (N, H, W)
        arrays.append(field)

    # Stack channels: list of (N, H, W) -> (N, C, H, W)
    data = np.stack(arrays, axis=1)
    return data


# ============================================================
# Normalization
# ============================================================

def normalize_per_snapshot(data):
    """
    Per-snapshot, per-channel min-max normalization to [0, 1].

    Parameters
    ----------
    data : np.ndarray, shape (N, C, H, W)

    Returns
    -------
    data_norm : np.ndarray, same shape, values in [0, 1]
    """
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


train_data = load_snapshots_from_h5(DATABASE_DIR, CHANNELS, train_indices)
val_data   = load_snapshots_from_h5(DATABASE_DIR, CHANNELS, val_indices)

if NORMALIZE:
    print("[INFO] Applying per-snapshot min-max normalization to [0, 1]")
    train_data = normalize_per_snapshot(train_data)
    val_data   = normalize_per_snapshot(val_data)
    print(f"  Train range: [{train_data.min():.4f}, {train_data.max():.4f}]")
    print(f"  Val range:   [{val_data.min():.4f}, {val_data.max():.4f}]")

print(f"  Train: {train_data.shape}  ({train_data.nbytes / 1e6:.1f} MB)")
print(f"  Val:   {val_data.shape}  ({val_data.nbytes / 1e6:.1f} MB)")

# Convert to PyTorch tensors and DataLoaders
train_tensor = torch.tensor(train_data, dtype=torch.float32)
val_tensor   = torch.tensor(val_data, dtype=torch.float32)

train_loader = DataLoader(TensorDataset(train_tensor),
                          batch_size=BATCH_SIZE, shuffle=True)
val_loader   = DataLoader(TensorDataset(val_tensor),
                          batch_size=BATCH_SIZE, shuffle=False)

# ============================================================
# Model
# ============================================================
print("\n" + "=" * 60)
print("Creating Model")
print("=" * 60)

model = ConvolutionalAutoencoder(latent_dim=LATENT_DIM,
                                 in_channels=n_channels).to(device)

total_params, trainable_params = count_parameters(model)
H, W = spatial_shape
print(f"  Channels:       {n_channels}")
print(f"  Latent dim:     {LATENT_DIM}")
print(f"  Compression:    {n_channels * H * W:,} -> {LATENT_DIM}  "
      f"({n_channels * H * W // LATENT_DIM:,}:1)")
print(f"  Parameters:     {trainable_params:,}")

# ============================================================
# run_info.json
# ============================================================
run_info = {
    'dataset': split['dataset'],
    'channels': CHANNELS,
    'n_channels': n_channels,
    'spatial_shape': list(spatial_shape),
    'latent_dim': LATENT_DIM,
    'compression_ratio': float(n_channels * H * W / LATENT_DIM),
    'batch_size': BATCH_SIZE,
    'learning_rate': LEARNING_RATE,
    'num_epochs': NUM_EPOCHS,
    'patience': PATIENCE,
    'seed': SEED,
    'total_parameters': total_params,
    'trainable_parameters': trainable_params,
    'train_cases': split['train_cases'],
    'val_cases': split['val_cases'],
    'test_cases': split['test_cases'],
    'n_train_snapshots': len(train_indices),
    'n_val_snapshots': len(val_indices),
    'n_test_snapshots': len(split['test_indices']),
    'split_file': SPLIT_FILE,
    'normalize': NORMALIZE,
    'normalization_type': 'per_snapshot_minmax_01' if NORMALIZE else 'none',
}
with open(os.path.join(results_directory, 'run_info.json'), 'w') as f:
    json.dump(run_info, f, indent=2)
print(f"[INFO] Saved run_info.json")

# ============================================================
# Loss & Optimizer
# ============================================================
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)

print(f"\n  Loss:       MSE")
print(f"  Optimizer:  Adam (lr={LEARNING_RATE})")
print(f"  Epochs:     {NUM_EPOCHS}")
print(f"  Patience:   {PATIENCE}")

# ============================================================
# GPU Memory Summary
# ============================================================
if device.type == 'cuda':
    with torch.no_grad():
        _ = model(train_tensor[:1].to(device))
    with open(os.path.join(results_directory, 'memory_summary.txt'), 'w') as f:
        f.write(torch.cuda.memory_summary())
    print(f"[INFO] Saved memory_summary.txt")

# ============================================================
# Training Loop
# ============================================================
print("\n" + "=" * 60)
print("Training")
print("=" * 60)

train_losses = []
val_losses   = []
best_val_loss = float('inf')
patience_counter = 0
start_time = time.time()

metric_types = [MetricType.MAE, MetricType.MSE, MetricType.LINF]
eval_every   = max(1, NUM_EPOCHS // 10)

for epoch in range(NUM_EPOCHS):

    # --- Train ---
    model.train()
    epoch_train_loss = 0.0
    for (x,) in train_loader:
        x = x.to(device)
        x_recon = model(x)
        loss = criterion(x_recon, x)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch_train_loss += loss.item()
    avg_train_loss = epoch_train_loss / len(train_loader)
    train_losses.append(avg_train_loss)

    # --- Validate ---
    model.eval()
    epoch_val_loss = 0.0
    with torch.no_grad():
        for (x,) in val_loader:
            x = x.to(device)
            x_recon = model(x)
            loss = criterion(x_recon, x)
            epoch_val_loss += loss.item()
    avg_val_loss = epoch_val_loss / len(val_loader)
    val_losses.append(avg_val_loss)

    # --- Early stopping ---
    if avg_val_loss < best_val_loss:
        best_val_loss = avg_val_loss
        patience_counter = 0
        torch.save({'epoch': epoch + 1,
                     'model': model.state_dict(),
                     'opt': optimizer.state_dict()},
                    os.path.join(results_directory, 'model-best.pt'))
    else:
        patience_counter += 1

    # --- Print progress ---
    if (epoch + 1) % 100 == 0 or epoch == 0:
        elapsed = time.time() - start_time
        print(f"  Epoch {epoch+1:5d}/{NUM_EPOCHS} | "
              f"Train: {avg_train_loss:.6f} | "
              f"Val: {avg_val_loss:.6f} | "
              f"Time: {elapsed:.1f}s")

    # --- Periodic metrics ---
    if (epoch + 1) % eval_every == 0 or epoch == 0 or epoch == NUM_EPOCHS - 1:
        model.eval()

        # Validation metrics
        val_metrics = []
        with torch.no_grad():
            for (x,) in val_loader:
                x = x.to(device)
                x_recon = model(x)
                for j in range(x.shape[0]):
                    for ch in range(n_channels):
                        orig = x[j, ch].cpu().numpy()
                        pred = x_recon[j, ch].cpu().numpy()
                        row = compute_metrics_single_array(orig, pred, metric_types)
                        row['channel'] = ch
                        val_metrics.append(row)
        pd.DataFrame(val_metrics).to_csv(
            os.path.join(results_directory, 'metrics', f'val_metrics_{epoch+1}.csv'),
            index=False)

        # Train metrics (subset)
        train_metrics = []
        with torch.no_grad():
            for i, (x,) in enumerate(train_loader):
                if i >= 3:
                    break
                x = x.to(device)
                x_recon = model(x)
                for j in range(x.shape[0]):
                    for ch in range(n_channels):
                        orig = x[j, ch].cpu().numpy()
                        pred = x_recon[j, ch].cpu().numpy()
                        row = compute_metrics_single_array(orig, pred, metric_types)
                        row['channel'] = ch
                        train_metrics.append(row)
        pd.DataFrame(train_metrics).to_csv(
            os.path.join(results_directory, 'metrics', f'train_metrics_{epoch+1}.csv'),
            index=False)

        # Intensity histograms
        preds_list, data_list = [], []
        with torch.no_grad():
            for i, (x,) in enumerate(val_loader):
                if i >= 3:
                    break
                x = x.to(device)
                x_recon = model(x)
                data_list.append(x.cpu().numpy())
                preds_list.append(x_recon.cpu().numpy())
        all_data  = np.concatenate(data_list)
        all_preds = np.concatenate(preds_list)

        fig, axes = plt.subplots(1, n_channels, figsize=(3 * n_channels, 3), dpi=150)
        if n_channels == 1:
            axes = [axes]
        for ch, ax in enumerate(axes):
            ax.hist(all_data[:, ch].flatten(), bins=100, alpha=0.5,
                    label='Data', color='blue')
            ax.hist(all_preds[:, ch].flatten(), bins=100, alpha=0.5,
                    label='Pred', color='red')
            ax.set_yscale('log')
            ax.set_xlabel('Intensity')
            ax.set_ylabel('Frequency')
            ax.set_title(f'Ch {CHANNELS[ch]}')
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(results_directory, 'intensity_histograms',
                                 f'hist_val_{epoch+1}.png'))
        plt.close(fig)

    # --- Stop ---
    if patience_counter >= PATIENCE:
        print(f"\n  Early stopping at epoch {epoch+1} "
              f"(no improvement for {PATIENCE} epochs)")
        break

total_time = time.time() - start_time
print(f"\n  Training complete! Time: {total_time:.1f}s")
print(f"  Best val loss: {best_val_loss:.6f}")

# ============================================================
# Save model_final.pt (best weights)
# ============================================================
best_ckpt = torch.load(os.path.join(results_directory, 'model-best.pt'),
                        map_location=device)
model.load_state_dict(best_ckpt['model'])

model_weights_path = os.path.join(results_directory, 'model_final.pt')
torch.save(model.state_dict(), model_weights_path)
print(f"[INFO] Saved model_final.pt (best epoch: {best_ckpt['epoch']})")

# ============================================================
# Loss History
# ============================================================
loss_history = {'Train': train_losses, 'Val': val_losses}
with open(os.path.join(results_directory, 'loss_history.json'), 'w') as f:
    json.dump(loss_history, f)

fig, ax = plt.subplots(figsize=(3, 3), dpi=200)
ax.plot(train_losses, label='Train')
ax.plot(val_losses, label='Val')
ax.axvline(best_ckpt['epoch'], color='red', ls='--', alpha=0.5,
           label=f"Best (ep {best_ckpt['epoch']})")
ax.set_yscale('log')
ax.set_xlabel('Epoch')
ax.set_ylabel('Loss')
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(results_directory, 'loss_history.png'))
plt.close(fig)
print(f"[INFO] Saved loss_history.json + loss_history.png")

# ============================================================
# Normalization Info
# ============================================================
with open(os.path.join(results_directory, 'normalization_info.txt'), 'w') as f:
    f.write(f"Dataset: {split['dataset']}\n")
    f.write(f"Channels: {CHANNELS}\n")
    if NORMALIZE:
        f.write("Normalization: per-snapshot min-max to [0, 1]\n")
        f.write("  x_norm = (x - x_min) / (x_max - x_min)  per snapshot, per channel\n")
    else:
        f.write("No normalization applied (raw field values)\n")
    f.write(f"Spatial grid: {H} x {W}\n")
    f.write(f"Training: {len(train_indices)} snapshots "
            f"({len(split['train_cases'])} cases)\n")
    f.write(f"Validation: {len(val_indices)} snapshots "
            f"({len(split['val_cases'])} cases)\n")
print(f"[INFO] Saved normalization_info.txt")

# ============================================================
# Reconstruction Visualization
# ============================================================
print("\n" + "=" * 60)
print("Reconstruction Samples")
print("=" * 60)

model.eval()
# Pick a few val samples
n_show = min(4, val_tensor.shape[0])
show_idx = np.linspace(0, val_tensor.shape[0] - 1, n_show, dtype=int)

fig, axes = plt.subplots(2 * n_channels, n_show,
                         figsize=(3 * n_show, 3 * 2 * n_channels))
if n_show == 1:
    axes = axes[:, np.newaxis]

with torch.no_grad():
    for col, idx in enumerate(show_idx):
        x = val_tensor[idx:idx+1].to(device)
        x_recon = model(x)
        for ch in range(n_channels):
            orig_np  = x[0, ch].cpu().numpy()
            recon_np = x_recon[0, ch].cpu().numpy()
            vmin, vmax = orig_np.min(), orig_np.max()

            row_orig  = ch * 2
            row_recon = ch * 2 + 1

            axes[row_orig, col].imshow(orig_np.T, origin='lower',
                                       cmap='plasma', vmin=vmin, vmax=vmax)
            axes[row_orig, col].axis('off')
            if col == 0:
                axes[row_orig, col].set_ylabel(f'Ch{CHANNELS[ch]} Orig', fontsize=9)

            axes[row_recon, col].imshow(recon_np.T, origin='lower',
                                        cmap='plasma', vmin=vmin, vmax=vmax)
            axes[row_recon, col].axis('off')
            if col == 0:
                axes[row_recon, col].set_ylabel(f'Ch{CHANNELS[ch]} Recon', fontsize=9)

            if ch == 0:
                mse = np.mean((orig_np - recon_np) ** 2)
                axes[row_orig, col].set_title(f'Val[{idx}]\nMSE={mse:.2e}',
                                              fontsize=8)

plt.suptitle('CAE Reconstruction (Validation Set)', fontsize=12)
plt.tight_layout()
fig.savefig(os.path.join(results_directory, 'reconstruction_test.png'), dpi=150)
plt.close(fig)
print(f"[INFO] Saved reconstruction_test.png")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
print(f"[INFO] All results saved to: {results_directory}")
