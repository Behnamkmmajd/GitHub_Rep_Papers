"""
Phase 1 — Prototype: Validate circularity & rise velocity on a single case.

Usage:
    conda run -p ~/miniconda3 --no-capture-output python prototype_single_case.py

Outputs:
    figures/prototype_case{CASE_ID}_fields.png
    figures/prototype_case{CASE_ID}_timeseries.png
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage import measure

# ─── Configuration ──────────────────────────────────────────────────────────
DATABASE_DIR = os.path.join(os.path.dirname(__file__), '..', 'database', 'MPF')
CASE_ID = 1775       # which case to prototype on (calm case: Re=17, Bo=11)
DT = 0.01            # non-dimensional time step between snapshots
H, W = 512, 256      # grid dimensions (height, width)
SNAPSHOTS_PER_CASE = 101
C_THRESHOLD = 0.0    # phase-field threshold: c < thr → bubble  (field in [-0.5, 0.5])
MIN_BUBBLE_AREA = 50 # minimum area (pixels) to count as a real fragment

os.makedirs(os.path.join(os.path.dirname(__file__), 'figures'), exist_ok=True)


# ─── Helper functions ───────────────────────────────────────────────────────

def circular_mean_y(mask, H):
    """Compute the center-of-mass y-coordinate on a periodic [0, H) domain.

    Uses the circular-mean trick: map y → angle, average sin/cos, map back.
    """
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
    """Roll a 2D array vertically so that y_center sits at row H//2."""
    shift = int(round(H // 2 - y_center))
    return np.roll(field_2d, shift, axis=0), shift


def largest_component_mask(binary, min_area=MIN_BUBBLE_AREA):
    """Return mask of the largest connected component (area >= min_area)."""
    labeled, n = ndimage.label(binary)
    if n == 0:
        return np.zeros_like(binary, dtype=bool), 0
    sizes = ndimage.sum(binary, labeled, range(1, n + 1))
    largest_label = np.argmax(sizes) + 1
    # Count fragments above threshold
    n_fragments = int(np.sum(np.array(sizes) >= min_area))
    return labeled == largest_label, n_fragments


def compute_circularity(component_mask):
    """Circularity = 4π·A / P². Returns NaN if no valid region.
    
    Uses regionprops for robust area and perimeter computation.
    Perimeter uses the corrected method (accounts for diagonal adjacency).
    """
    labeled_single = ndimage.label(component_mask)[0]
    props = measure.regionprops(labeled_single)
    if len(props) == 0:
        return np.nan
    # Take the largest region (should be only one after largest_component_mask)
    region = max(props, key=lambda r: r.area)
    area = region.area
    if area < MIN_BUBBLE_AREA:
        return np.nan
    perimeter = region.perimeter
    if perimeter == 0:
        return np.nan
    return 4.0 * np.pi * area / (perimeter ** 2)


def process_snapshot(phase_field, H):
    """Process a single snapshot → (circularity, y_com, n_fragments)."""
    bubble_mask = phase_field < C_THRESHOLD

    # Circular center of mass
    y_com = circular_mean_y(bubble_mask, H)

    # Roll field so bubble is centered (avoid split across periodic boundary)
    rolled_mask, _ = roll_to_center(bubble_mask, y_com, H)

    # Largest connected component
    comp_mask, n_frags = largest_component_mask(rolled_mask)

    # Circularity of the main bubble
    circ = compute_circularity(comp_mask)

    return circ, y_com, n_frags


def unwrap_y_com(y_com_series, H):
    """Unwrap the y center-of-mass series to remove periodic jumps."""
    unwrapped = np.copy(y_com_series)
    for i in range(1, len(unwrapped)):
        diff = unwrapped[i] - unwrapped[i - 1]
        if diff > H / 2:
            unwrapped[i:] -= H
        elif diff < -H / 2:
            unwrapped[i:] += H
    return unwrapped


# ─── Main ───────────────────────────────────────────────────────────────────

def main():
    h5_path = os.path.join(DATABASE_DIR, 'channel_1.h5')
    print(f"Opening {h5_path} ...")

    with h5py.File(h5_path, 'r') as f:
        # Load case params
        params = f['case_params'][CASE_ID]
        case_name = f['case_names'][CASE_ID]
        if isinstance(case_name, bytes):
            case_name = case_name.decode()
        print(f"Case {CASE_ID}: {case_name}")
        print(f"  DR={params[0]}, VR={params[1]}, Re={params[2]}, Bo={params[3]}")

        # Load all 101 snapshots for this case
        start_idx = CASE_ID * SNAPSHOTS_PER_CASE
        end_idx = start_idx + SNAPSHOTS_PER_CASE
        print(f"  Loading field[{start_idx}:{end_idx}] ...")
        snapshots = f['field'][start_idx:end_idx]  # (101, 512, 256)

    print(f"  Loaded shape: {snapshots.shape}, dtype: {snapshots.dtype}")
    print(f"  Value range: [{snapshots.min():.4f}, {snapshots.max():.4f}]")

    # ── Process each snapshot ──
    circularity = np.zeros(SNAPSHOTS_PER_CASE)
    y_com = np.zeros(SNAPSHOTS_PER_CASE)
    n_fragments = np.zeros(SNAPSHOTS_PER_CASE, dtype=int)

    for j in range(SNAPSHOTS_PER_CASE):
        circ, yc, nf = process_snapshot(snapshots[j], H)
        circularity[j] = circ
        y_com[j] = yc
        n_fragments[j] = nf
        if j % 20 == 0:
            print(f"  snapshot {j:3d}: circ={circ:.4f}, y_com={yc:.1f}, "
                  f"n_frag={nf}")

    # Unwrap y_com and compute rise velocity
    y_com_unwrapped = unwrap_y_com(y_com, H)
    rise_velocity = np.gradient(y_com_unwrapped, DT)

    time = np.arange(SNAPSHOTS_PER_CASE) * DT

    # ── Figure 1: Phase field snapshots ──
    fig, axes = plt.subplots(1, 5, figsize=(20, 8))
    snapshot_indices = [0, 25, 50, 75, 100]
    for ax, si in zip(axes, snapshot_indices):
        im = ax.imshow(snapshots[si], origin='lower', cmap='RdBu_r',
                       vmin=0, vmax=1, aspect='equal')
        ax.set_title(f't = {si * DT:.2f}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        # Overlay c=0.5 contour
        contours = measure.find_contours(snapshots[si].astype(np.float64),
                                              C_THRESHOLD)
        for c in contours:
            ax.plot(c[:, 1], c[:, 0], 'k-', linewidth=0.8)
    fig.colorbar(im, ax=axes, shrink=0.6, label='Phase field c')
    fig.suptitle(f'Case {CASE_ID}: {case_name}\n'
                 f'DR={params[0]}, VR={params[1]}, Re={params[2]}, Bo={params[3]}',
                 fontsize=13)
    fig.tight_layout()
    fig_path = os.path.join(os.path.dirname(__file__), 'figures',
                            f'prototype_case{CASE_ID}_fields.png')
    fig.savefig(fig_path, dpi=150)
    print(f"\nSaved: {fig_path}")
    plt.close(fig)

    # ── Figure 2: Time series ──
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(time, circularity, 'b-o', markersize=2)
    axes[0].set_ylabel('Circularity')
    axes[0].set_ylim(0, 1.1)
    axes[0].axhline(1.0, color='gray', ls='--', alpha=0.5, label='perfect circle')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(time, rise_velocity, 'r-o', markersize=2)
    axes[1].set_ylabel('Rise velocity (dy_com/dt)')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(time, n_fragments, 'g-o', markersize=2)
    axes[2].set_ylabel('# fragments')
    axes[2].set_xlabel('Time (non-dim)')
    axes[2].set_ylim(0, max(n_fragments.max() + 1, 3))
    axes[2].grid(True, alpha=0.3)

    fig.suptitle(f'Case {CASE_ID}: {case_name}\n'
                 f'DR={params[0]}, VR={params[1]}, Re={params[2]}, Bo={params[3]}',
                 fontsize=13)
    fig.tight_layout()
    fig_path = os.path.join(os.path.dirname(__file__), 'figures',
                            f'prototype_case{CASE_ID}_timeseries.png')
    fig.savefig(fig_path, dpi=150)
    print(f"Saved: {fig_path}")
    plt.close(fig)

    # ── Summary ──
    print("\n─── Summary ───────────────────────────────────────────")
    print(f"  Circularity range: [{np.nanmin(circularity):.4f}, "
          f"{np.nanmax(circularity):.4f}]")
    print(f"  Rise velocity range: [{np.nanmin(rise_velocity):.2f}, "
          f"{np.nanmax(rise_velocity):.2f}]")
    print(f"  Max fragments: {n_fragments.max()}")
    breakup_snaps = np.where(n_fragments > 1)[0]
    if len(breakup_snaps) > 0:
        print(f"  Breakup detected at snapshots: {breakup_snaps}")
    else:
        print(f"  No breakup detected")

    # Quick steady-state check: rolling rel. std over last 10 snapshots
    W = 10
    if SNAPSHOTS_PER_CASE > W:
        tail = circularity[-W:]
        rel_std = np.nanstd(tail) / max(np.nanmean(tail), 1e-12)
        print(f"  Circularity rel. std (last {W} snapshots): {rel_std:.4f} "
              f"({'STEADY' if rel_std < 0.01 else 'NOT STEADY'} at 1% thr)")


if __name__ == '__main__':
    main()
