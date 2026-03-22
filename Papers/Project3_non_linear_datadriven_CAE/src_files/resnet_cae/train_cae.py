# CAE Training Script — KTH ResNet Autoencoder
# ============================================================
# Trains a ResNet-based Convolutional Autoencoder on HDF5 field
# data using a standardized split.json for train/val/test
# partitioning.  Works with any dataset and any number of channels.
#
# Usage:  python train_cae.py   (run from resnet_cae/)
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
DIM            = 16                     # Base feature dimension
DIM_MULTS      = (1, 2, 4, 8, 16, 32)             # 6 down/up levels
Z_CHANNELS     = 4                      # Latent feature channels
BLOCK_TYPE     = 1                      # 0=Basic, 1=ResNet

# --- Normalization ---
NORMALIZE      = True                   # Per-snapshot min-max to [0, 1]

# --- Training ---
BATCH_SIZE     = 8
LEARNING_RATE  = 1e-5
L2_REG         = 1e-6
NUM_EPOCHS     = 300
SAVE_EVERY     = 2                      # Evaluate metrics every N epochs

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
import numpy as np
import torch
import h5py
from torch.utils.data import Dataset
from torch import nn

# ============================================================
# Load AE Modules
# ============================================================
sys.path.append(os.path.join(os.path.dirname(__file__), 'library_resnet_cae'))

from resnet_cae.models.baseline_model import ConvAutoencoderBaseline                   # type: ignore
from resnet_cae.trainer_ae import MyAETrainer                                          # type: ignore
from resnet_cae.metrics import MetricType                                              # type: ignore

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
spatial_shape       = tuple(split['spatial_shape'])       # (H, W)
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
# Dataset Class
# ============================================================
class FieldDataset(Dataset):
    def __init__(self, data: torch.Tensor):
        """data: (N, C, H, W) float32"""
        self.data_norm = data

    def __len__(self):
        return self.data_norm.shape[0]

    def __getitem__(self, idx):
        return self.data_norm[idx]

    def unnormalise_array(self, array):
        """Identity (raw field values, no normalization applied)."""
        return array

# ============================================================
# Results Directory
# ============================================================
np.random.seed(SEED)
torch.manual_seed(SEED)

ch_tag = "".join(str(c) for c in CHANNELS)
norm_tag = "_Norm" if NORMALIZE else ""
run_name = (f"{split['dataset']}_Ch{ch_tag}_"
            f"Nt{len(train_indices)}_B{BATCH_SIZE}_"
            f"NL{len(DIM_MULTS)}_LR{LEARNING_RATE}_L2{L2_REG}_E{NUM_EPOCHS}{norm_tag}")
results_directory = os.path.join('runs', run_name)
os.makedirs(results_directory, exist_ok=True)
print(f"[INFO] Results:    {results_directory}")

# ============================================================
# Load Data from HDF5
# ============================================================
print("\n" + "=" * 60)
print("Loading Data")
print("=" * 60)


def load_snapshots_from_h5(database_dir, channels, indices):
    """Load selected rows from channel HDF5 files → (N, C, H, W)."""
    arrays = []
    for ch in channels:
        h5_path = os.path.join(database_dir, f'channel_{ch}.h5')
        with h5py.File(h5_path, 'r') as f:
            field = f['field'][sorted(indices)]
        arrays.append(field)
    return np.stack(arrays, axis=1)


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

train_tensor = torch.tensor(train_data, dtype=torch.float32)
val_tensor   = torch.tensor(val_data, dtype=torch.float32)

train_ds = FieldDataset(train_tensor)
val_ds   = FieldDataset(val_tensor)

print(f"[INFO] Datasets ready: train={len(train_ds)}, val={len(val_ds)}")

# ============================================================
# Model Definition
# ============================================================
print("\n" + "=" * 60)
print("Creating Model")
print("=" * 60)

model = ConvAutoencoderBaseline(
    dim=DIM,
    dim_mults=DIM_MULTS,
    channels=n_channels,
    z_channels=Z_CHANNELS,
    block_type=BLOCK_TYPE
).float()

n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
H, W = spatial_shape
print(f"  Channels:    {n_channels}")
print(f"  Dim:         {DIM}")
print(f"  Levels:      {len(DIM_MULTS)}")
print(f"  z_channels:  {Z_CHANNELS}")
print(f"  Parameters:  {n_params:,}")

# ============================================================
# Hardware Setup
# ============================================================
if torch.cuda.is_available():
    torch.cuda.set_device(0)
    _ = torch.randn(1, device="cuda")
    model = model.to("cuda")

    if torch.cuda.device_count() > 1:
        print(f"[INFO] Using {torch.cuda.device_count()} GPUs with DataParallel.")
        model = torch.nn.DataParallel(model)
    else:
        print(f"[INFO] Using 1 GPU: {torch.cuda.get_device_name(0)}")
    CPU_ONLY = False
else:
    print("[INFO] No GPU detected. Using CPU only.")
    CPU_ONLY = True

# ============================================================
# Training
# ============================================================
print("\n" + "=" * 60)
print("Training")
print("=" * 60)

trainer = MyAETrainer(
    model=model,
    dataset_train=train_ds,
    dataset_val=val_ds,
    train_batch_size=BATCH_SIZE,
    train_lr=LEARNING_RATE,
    l2_reg=L2_REG,
    train_num_epochs=NUM_EPOCHS,
    save_and_sample_every=SAVE_EVERY,
    results_folder=results_directory,
    loss=nn.MSELoss(),
    metric_types=[MetricType.MAE, MetricType.MSE, MetricType.LINF],
    cpu_only=CPU_ONLY,
    low_data_mode=True,
    num_dl_workers=0
)

print("[INFO] Starting training...")
trainer.train()
print("[INFO] Training finished. Final mean val metrics:", trainer.mean_val_metrics)

# ============================================================
# Save model_final.pt & Normalization Info
# ============================================================
model_weights_path = os.path.join(results_directory, "model_final.pt")
torch.save(model.state_dict(), model_weights_path)
print(f"[INFO] Saved model_final.pt")

with open(os.path.join(results_directory, "normalization_info.txt"), "w") as f:
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

# ============================================================
# run_info supplement (trainer already writes basic run_info.json)
# ============================================================
run_info_path = os.path.join(results_directory, 'run_info.json')
with open(run_info_path, 'r') as f:
    run_info = json.load(f)

run_info.update({
    'dataset': split['dataset'],
    'channels': CHANNELS,
    'n_channels': n_channels,
    'spatial_shape': list(spatial_shape),
    'dim': DIM,
    'dim_mults': list(DIM_MULTS),
    'z_channels': Z_CHANNELS,
    'block_type': BLOCK_TYPE,
    'batch_size': BATCH_SIZE,
    'learning_rate': LEARNING_RATE,
    'l2_reg': L2_REG,
    'num_epochs': NUM_EPOCHS,
    'seed': SEED,
    'trainable_parameters': n_params,
    'train_cases': split['train_cases'],
    'val_cases': split['val_cases'],
    'test_cases': split['test_cases'],
    'n_train_snapshots': len(train_indices),
    'n_val_snapshots': len(val_indices),
    'n_test_snapshots': len(split['test_indices']),
    'split_file': SPLIT_FILE,
    'normalize': NORMALIZE,
    'normalization_type': 'per_snapshot_minmax_01' if NORMALIZE else 'none',
})
with open(run_info_path, 'w') as f:
    json.dump(run_info, f, indent=2)

print(f"[INFO] Updated run_info.json")
print(f"[INFO] Saved normalization_info.txt")
print(f"\n[INFO] All results saved to: {results_directory}")
