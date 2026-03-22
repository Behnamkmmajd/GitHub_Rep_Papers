"""
Phase 2 — Batch processing: Compute circularity, rise velocity, and fragment
count for all 5,000 cases from channel_1.h5 (phase field).

Usage:
    conda run -p ~/miniconda3 --no-capture-output python batch_process_all_cases.py

Output:
    results/steady_state_metrics.npz
        - circularity      (5000, 101)  float32
        - y_com            (5000, 101)  float32  (unwrapped)
        - rise_velocity    (5000, 100)  float32  (finite differences, 100 intervals)
        - n_fragments      (5000, 101)  int16
        - case_params      (5000, 4)    int32    [DR, VR, Re, Bo]
        - case_names       (5000,)      str
        - dt               scalar       float64

Runtime estimate: ~1–2 hours on HDD, mostly I/O-bound.
"""

import os
import sys
import time
import numpy as np
import h5py
from scipy import ndimage
from skimage import measure
from tqdm import tqdm

# ─── Configuration ──────────────────────────────────────────────────────────
DATABASE_DIR = os.path.join(os.path.dirname(__file__), '..', 'database', 'MPF')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
DT = 0.01
H, W = 512, 256
SNAPSHOTS_PER_CASE = 101
C_THRESHOLD = 0.0      # phase field in [-0.5, 0.5]; c < 0 → bubble
MIN_BUBBLE_AREA = 50   # pixels

# Batch size: how many cases to read at once (balance RAM vs I/O efficiency)
# 50 cases × 101 snapshots × 512 × 256 × 4 bytes ≈ 2.6 GB
BATCH_SIZE = 50

os.makedirs(RESULTS_DIR, exist_ok=True)


# ─── Core functions (same as prototype, inlined for standalone use) ─────────

def circular_mean_y(mask, H):
    ys = np.where(mask)[0].astype(np.float64)
    if len(ys) == 0:
        return np.nan
    theta = 2.0 * np.pi * ys / H
    y_center = H / (2.0 * np.pi) * np.arctan2(np.mean(np.sin(theta)),
                                                 np.mean(np.cos(theta)))
    if y_center < 0:
        y_center += H
    return y_center


def roll_to_center(field_2d, y_center, H):
    shift = int(round(H // 2 - y_center))
    return np.roll(field_2d, shift, axis=0)


def largest_component_mask(binary, min_area=MIN_BUBBLE_AREA):
    labeled, n = ndimage.label(binary)
    if n == 0:
        return np.zeros_like(binary, dtype=bool), 0
    sizes = ndimage.sum(binary, labeled, range(1, n + 1))
    largest_label = np.argmax(sizes) + 1
    n_fragments = int(np.sum(np.array(sizes) >= min_area))
    return labeled == largest_label, n_fragments


def compute_circularity(component_mask):
    labeled_single = ndimage.label(component_mask)[0]
    props = measure.regionprops(labeled_single)
    if len(props) == 0:
        return np.nan
    region = max(props, key=lambda r: r.area)
    area = region.area
    if area < MIN_BUBBLE_AREA:
        return np.nan
    perimeter = region.perimeter
    if perimeter == 0:
        return np.nan
    return 4.0 * np.pi * area / (perimeter ** 2)


def process_snapshot(phase_field):
    bubble_mask = phase_field < C_THRESHOLD
    y_com = circular_mean_y(bubble_mask, H)
    rolled_mask = roll_to_center(bubble_mask, y_com, H)
    comp_mask, n_frags = largest_component_mask(rolled_mask)
    circ = compute_circularity(comp_mask)
    return circ, y_com, n_frags


def unwrap_y_com(y_com_series, H):
    unwrapped = np.copy(y_com_series)
    for i in range(1, len(unwrapped)):
        diff = unwrapped[i] - unwrapped[i - 1]
        if diff > H / 2:
            unwrapped[i:] -= H
        elif diff < -H / 2:
            unwrapped[i:] += H
    return unwrapped


def process_case(snapshots):
    """Process all snapshots for one case.

    Args:
        snapshots: (101, 512, 256) float32 array

    Returns:
        circularity:    (101,) float32
        y_com:          (101,) float32  (unwrapped)
        rise_velocity:  (100,) float32  (forward differences)
        n_fragments:    (101,) int16
    """
    n = snapshots.shape[0]
    circularity = np.zeros(n, dtype=np.float32)
    y_com_raw = np.zeros(n, dtype=np.float64)
    n_fragments = np.zeros(n, dtype=np.int16)

    for j in range(n):
        circ, yc, nf = process_snapshot(snapshots[j])
        circularity[j] = circ
        y_com_raw[j] = yc
        n_fragments[j] = nf

    # Unwrap periodic y_com
    y_com = unwrap_y_com(y_com_raw, H).astype(np.float32)

    # Rise velocity via finite differences
    rise_velocity = np.diff(y_com) / DT  # (100,)

    return circularity, y_com, rise_velocity, n_fragments


# ─── Main ───────────────────────────────────────────────────────────────────

def main():
    h5_path = os.path.join(DATABASE_DIR, 'channel_1.h5')
    out_path = os.path.join(RESULTS_DIR, 'steady_state_metrics.npz')

    # Check for existing partial results to allow resumption
    start_case = 0
    with h5py.File(h5_path, 'r') as f:
        n_cases = int(f.attrs['n_cases'])
        case_params = f['case_params'][:]       # (5000, 4) int32
        case_names_raw = f['case_names'][:]     # (5000,)
        case_names = np.array([n.decode() if isinstance(n, bytes) else n
                               for n in case_names_raw])

    print(f"Dataset: {n_cases} cases, {SNAPSHOTS_PER_CASE} snapshots each")
    print(f"Output: {out_path}")
    print(f"Batch size: {BATCH_SIZE} cases (~{BATCH_SIZE * SNAPSHOTS_PER_CASE * H * W * 4 / 1e9:.1f} GB per batch)")

    # Allocate output arrays
    all_circularity = np.zeros((n_cases, SNAPSHOTS_PER_CASE), dtype=np.float32)
    all_y_com = np.zeros((n_cases, SNAPSHOTS_PER_CASE), dtype=np.float32)
    all_rise_vel = np.zeros((n_cases, SNAPSHOTS_PER_CASE - 1), dtype=np.float32)
    all_n_frags = np.zeros((n_cases, SNAPSHOTS_PER_CASE), dtype=np.int16)

    # Check for partial checkpoint
    ckpt_path = os.path.join(RESULTS_DIR, 'checkpoint.npz')
    if os.path.exists(ckpt_path):
        print(f"Found checkpoint: {ckpt_path}")
        ckpt = np.load(ckpt_path)
        start_case = int(ckpt['last_completed_case']) + 1
        if start_case < n_cases:
            all_circularity[:start_case] = ckpt['circularity'][:start_case]
            all_y_com[:start_case] = ckpt['y_com'][:start_case]
            all_rise_vel[:start_case] = ckpt['rise_velocity'][:start_case]
            all_n_frags[:start_case] = ckpt['n_fragments'][:start_case]
            print(f"Resuming from case {start_case}")
        else:
            print("Checkpoint shows all cases done. Saving final output.")

    t_start = time.time()
    n_batches = (n_cases - start_case + BATCH_SIZE - 1) // BATCH_SIZE

    with h5py.File(h5_path, 'r') as f:
        field_ds = f['field']
        pbar = tqdm(total=n_cases - start_case, desc="Processing cases",
                    unit="case")

        for batch_idx in range(n_batches):
            batch_start = start_case + batch_idx * BATCH_SIZE
            batch_end = min(batch_start + BATCH_SIZE, n_cases)
            batch_n = batch_end - batch_start

            # Read all snapshots for this batch of cases at once
            row_start = batch_start * SNAPSHOTS_PER_CASE
            row_end = batch_end * SNAPSHOTS_PER_CASE
            batch_data = field_ds[row_start:row_end]  # (batch_n*101, 512, 256)

            for local_i in range(batch_n):
                case_i = batch_start + local_i
                snap_start = local_i * SNAPSHOTS_PER_CASE
                snap_end = snap_start + SNAPSHOTS_PER_CASE
                case_snaps = batch_data[snap_start:snap_end]

                circ, ycom, vrise, nfrags = process_case(case_snaps)
                all_circularity[case_i] = circ
                all_y_com[case_i] = ycom
                all_rise_vel[case_i] = vrise
                all_n_frags[case_i] = nfrags

                pbar.update(1)

            # Save checkpoint every batch
            np.savez_compressed(
                ckpt_path,
                circularity=all_circularity,
                y_com=all_y_com,
                rise_velocity=all_rise_vel,
                n_fragments=all_n_frags,
                last_completed_case=batch_end - 1
            )

        pbar.close()

    elapsed = time.time() - t_start
    print(f"\nDone in {elapsed / 60:.1f} min ({elapsed / 3600:.2f} hr)")
    print(f"Rate: {(n_cases - start_case) / elapsed:.1f} cases/sec, "
          f"{(n_cases - start_case) * SNAPSHOTS_PER_CASE / elapsed:.0f} "
          f"snapshots/sec")

    # Save final output
    np.savez_compressed(
        out_path,
        circularity=all_circularity,
        y_com=all_y_com,
        rise_velocity=all_rise_vel,
        n_fragments=all_n_frags,
        case_params=case_params,
        case_names=case_names,
        dt=DT
    )
    print(f"Saved: {out_path}")

    # Clean up checkpoint
    if os.path.exists(ckpt_path):
        os.remove(ckpt_path)
        print("Removed checkpoint file.")

    # Quick stats
    print("\n─── Quick Stats ───────────────────────────────────────")
    valid = ~np.isnan(all_circularity[:, -1])
    print(f"  Cases with valid circularity at final snapshot: "
          f"{valid.sum()} / {n_cases}")
    breakup_cases = np.any(all_n_frags > 1, axis=1).sum()
    print(f"  Cases with breakup (any snapshot): {breakup_cases} / {n_cases}")


if __name__ == '__main__':
    main()
