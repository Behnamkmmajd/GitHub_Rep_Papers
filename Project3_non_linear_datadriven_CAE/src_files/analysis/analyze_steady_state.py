"""
Phase 3+4 — Steady-state detection + visualization.

Reads the output of batch_process_all_cases.py and produces:
  1. Per-case steady-state classification (steady / transient / breakup)
  2. Histograms of N_steady (number of steady-state snapshots)
  3. Parameter-regime breakdowns (Re vs Bo, DR vs VR heatmaps)
  4. Threshold sensitivity sweep
  5. Example gallery of representative cases

Usage:
    conda run -p ~/miniconda3 --no-capture-output python analyze_steady_state.py

Inputs:
    results/steady_state_metrics.npz  (from batch_process_all_cases.py)

Outputs:
    figures/histogram_n_steady.png
    figures/histogram_steady_onset.png
    figures/classification_pie.png
    figures/regime_Re_vs_Bo.png
    figures/regime_DR_vs_VR.png
    figures/threshold_sensitivity.png
    figures/gallery_examples.png
    results/steady_state_labels.npz
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

# ─── Configuration ──────────────────────────────────────────────────────────
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), 'figures')
DATABASE_DIR = os.path.join(os.path.dirname(__file__), '..', 'database', 'MPF')

SNAPSHOTS_PER_CASE = 101
DT = 0.01

# Steady-state detection parameters (primary)
WINDOW = 10          # rolling window size (snapshots)
THRESHOLD = 0.01     # relative std threshold (1%)

os.makedirs(FIGURES_DIR, exist_ok=True)


# ─── Steady-state detection ────────────────────────────────────────────────

def detect_steady_onset(signal, window, threshold):
    """Find the earliest snapshot from which the signal is 'steady' until the end.

    'Steady' = rolling relative std-dev < threshold for all remaining windows.

    Returns:
        onset: int or None (snapshot index where steady state begins)
        n_steady: int (number of steady-state snapshots, 0 if never steady)
    """
    n = len(signal)
    if n < window:
        return None, 0

    # Precompute rolling rel_std for each window ending at j
    is_steady = np.zeros(n, dtype=bool)
    for j in range(window - 1, n):
        seg = signal[j - window + 1:j + 1]
        seg_clean = seg[~np.isnan(seg)]
        if len(seg_clean) < window // 2:
            is_steady[j] = False
            continue
        mean_val = np.mean(seg_clean)
        if abs(mean_val) < 1e-12:
            is_steady[j] = False
            continue
        rel_std = np.std(seg_clean) / abs(mean_val)
        is_steady[j] = rel_std < threshold

    # Find onset: earliest j where all j..end are steady
    # Walk backward from end to find the last non-steady snapshot
    if not is_steady[-1]:
        return None, 0

    onset = n - 1
    for j in range(n - 1, window - 2, -1):
        if is_steady[j]:
            onset = j
        else:
            break

    # Convert onset from window-end index to actual snapshot index
    # onset is the window-end index; the snapshot range is [onset-window+1, end]
    # But for counting steady snapshots, we use onset as the start of steady region
    actual_onset = onset - window + 1  # first snapshot in the steady window
    n_steady = n - actual_onset

    return actual_onset, n_steady


def classify_cases(circularity, n_fragments, window, threshold):
    """Classify all cases as steady / transient / breakup.

    Returns:
        labels: (n_cases,) array of strings
        onset:  (n_cases,) int (-1 if not steady)
        n_steady: (n_cases,) int
    """
    n_cases = circularity.shape[0]
    labels = np.empty(n_cases, dtype='U12')
    onset_arr = np.full(n_cases, -1, dtype=int)
    n_steady_arr = np.zeros(n_cases, dtype=int)

    for i in range(n_cases):
        has_breakup = np.any(n_fragments[i] > 1)
        ss_onset, n_ss = detect_steady_onset(circularity[i], window, threshold)

        if ss_onset is not None and n_ss > 0:
            labels[i] = 'steady'
            onset_arr[i] = ss_onset
            n_steady_arr[i] = n_ss
        elif has_breakup:
            labels[i] = 'breakup'
        else:
            labels[i] = 'transient'

    return labels, onset_arr, n_steady_arr


# ─── Visualization ──────────────────────────────────────────────────────────

def plot_histogram_n_steady(n_steady, labels, threshold, window):
    fig, ax = plt.subplots(figsize=(10, 5))
    steady_mask = labels == 'steady'
    vals = n_steady[steady_mask]
    if len(vals) > 0:
        bins = np.arange(0, SNAPSHOTS_PER_CASE + 2) - 0.5
        ax.hist(vals, bins=bins, color='steelblue', edgecolor='white', alpha=0.8)
    ax.set_xlabel('Number of steady-state snapshots', fontsize=12)
    ax.set_ylabel('Number of cases', fontsize=12)
    ax.set_title(f'Distribution of steady-state snapshots\n'
                 f'(window={window}, threshold={threshold*100:.1f}%,'
                 f' {steady_mask.sum()} steady out of {len(labels)} cases)',
                 fontsize=13)
    ax.axvline(x=10, color='red', ls='--', alpha=0.6, label='10 snapshots')
    ax.axvline(x=20, color='orange', ls='--', alpha=0.6, label='20 snapshots')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'histogram_n_steady.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def plot_histogram_steady_onset(onset_arr, labels):
    fig, ax = plt.subplots(figsize=(10, 5))
    steady_mask = labels == 'steady'
    vals = onset_arr[steady_mask]
    vals = vals[vals >= 0]
    if len(vals) > 0:
        bins = np.arange(0, SNAPSHOTS_PER_CASE + 2) - 0.5
        ax.hist(vals, bins=bins, color='darkorange', edgecolor='white', alpha=0.8)
    ax.set_xlabel('Steady-state onset snapshot', fontsize=12)
    ax.set_ylabel('Number of cases', fontsize=12)
    ax.set_title('When does steady state begin?', fontsize=13)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'histogram_steady_onset.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def plot_classification_pie(labels):
    counts = {}
    for lab in ['steady', 'transient', 'breakup']:
        counts[lab] = np.sum(labels == lab)

    fig, ax = plt.subplots(figsize=(7, 7))
    colors = {'steady': 'steelblue', 'transient': 'darkorange', 'breakup': 'firebrick'}
    wedge_labels = [f"{k}\n({v}, {v/len(labels)*100:.1f}%)"
                    for k, v in counts.items()]
    ax.pie(list(counts.values()), labels=wedge_labels,
           colors=[colors[k] for k in counts],
           autopct='', startangle=90, textprops={'fontsize': 13})
    ax.set_title(f'Case classification (n={len(labels)})', fontsize=14)
    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'classification_pie.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def plot_regime_map(case_params, n_steady_arr, labels):
    """2D scatter: Re vs Bo and DR vs VR, colored by n_steady."""
    DR = case_params[:, 0]
    VR = case_params[:, 1]
    Re = case_params[:, 2]
    Bo = case_params[:, 3]

    # ── Re vs Bo ──
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Color by n_steady (0 for non-steady cases)
    sc = axes[0].scatter(Re, Bo, c=n_steady_arr, cmap='viridis', s=5, alpha=0.6,
                         norm=Normalize(vmin=0, vmax=SNAPSHOTS_PER_CASE))
    axes[0].set_xlabel('Re', fontsize=12)
    axes[0].set_ylabel('Bo', fontsize=12)
    axes[0].set_title('Re vs Bo — colored by N_steady', fontsize=13)
    fig.colorbar(sc, ax=axes[0], label='N steady-state snapshots')
    axes[0].grid(True, alpha=0.3)

    # ── DR vs VR ──
    sc2 = axes[1].scatter(DR, VR, c=n_steady_arr, cmap='viridis', s=5, alpha=0.6,
                          norm=Normalize(vmin=0, vmax=SNAPSHOTS_PER_CASE))
    axes[1].set_xlabel('DR (Density Ratio)', fontsize=12)
    axes[1].set_ylabel('VR (Viscosity Ratio)', fontsize=12)
    axes[1].set_title('DR vs VR — colored by N_steady', fontsize=13)
    fig.colorbar(sc2, ax=axes[1], label='N steady-state snapshots')
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'regime_scatter.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")

    # ── Categorical version: steady / transient / breakup ──
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    color_map = {'steady': 'steelblue', 'transient': 'darkorange', 'breakup': 'firebrick'}
    for lab, color in color_map.items():
        mask = labels == lab
        axes[0].scatter(Re[mask], Bo[mask], c=color, s=5, alpha=0.4, label=lab)
        axes[1].scatter(DR[mask], VR[mask], c=color, s=5, alpha=0.4, label=lab)

    axes[0].set_xlabel('Re', fontsize=12)
    axes[0].set_ylabel('Bo', fontsize=12)
    axes[0].set_title('Re vs Bo — by classification', fontsize=13)
    axes[0].legend(markerscale=3, fontsize=11)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel('DR (Density Ratio)', fontsize=12)
    axes[1].set_ylabel('VR (Viscosity Ratio)', fontsize=12)
    axes[1].set_title('DR vs VR — by classification', fontsize=13)
    axes[1].legend(markerscale=3, fontsize=11)
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'regime_classification.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def plot_threshold_sensitivity(circularity, n_fragments):
    """Sweep threshold and window to see robustness."""
    thresholds = [0.005, 0.01, 0.02, 0.03, 0.05]
    windows = [5, 10, 15, 20]

    results = np.zeros((len(windows), len(thresholds)), dtype=int)

    for wi, w in enumerate(windows):
        for ti, thr in enumerate(thresholds):
            labs, _, _ = classify_cases(circularity, n_fragments, w, thr)
            results[wi, ti] = np.sum(labs == 'steady')

    fig, ax = plt.subplots(figsize=(9, 6))
    for wi, w in enumerate(windows):
        ax.plot([t * 100 for t in thresholds], results[wi],
                '-o', label=f'window = {w}', markersize=6)

    ax.set_xlabel('Threshold (% relative std)', fontsize=12)
    ax.set_ylabel('Number of "steady" cases (out of 5000)', fontsize=12)
    ax.set_title('Threshold sensitivity: how many cases are classified as steady?',
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'threshold_sensitivity.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def plot_gallery(circularity, rise_velocity, n_fragments, case_params,
                 case_names, labels, onset_arr, n_steady_arr):
    """Pick 3 steady + 3 transient + 3 breakup cases and plot their time series."""
    import h5py as h5

    gallery_cases = {}
    for lab in ['steady', 'transient', 'breakup']:
        mask = labels == lab
        idxs = np.where(mask)[0]
        if len(idxs) == 0:
            continue
        # Pick 3 spread across the population
        if lab == 'steady':
            # Pick cases with high, medium, low n_steady
            sorted_idx = idxs[np.argsort(n_steady_arr[idxs])]
            picks = [sorted_idx[0], sorted_idx[len(sorted_idx)//2], sorted_idx[-1]]
        else:
            # Random spread
            rng = np.random.default_rng(42)
            picks = rng.choice(idxs, size=min(3, len(idxs)), replace=False)
        gallery_cases[lab] = list(picks)

    n_rows = sum(len(v) for v in gallery_cases.values())
    if n_rows == 0:
        print("No cases to show in gallery.")
        return

    time_circ = np.arange(SNAPSHOTS_PER_CASE) * DT
    time_vel = np.arange(SNAPSHOTS_PER_CASE - 1) * DT + DT / 2

    # Load phase field snapshots for gallery cases
    h5_path = os.path.join(DATABASE_DIR, 'channel_1.h5')

    fig, axes = plt.subplots(n_rows, 3, figsize=(18, 4.5 * n_rows))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    row = 0
    with h5.File(h5_path, 'r') as f:
        for lab, idxs in gallery_cases.items():
            for case_i in idxs:
                # Load initial and final snapshot
                snap_idx_init = case_i * SNAPSHOTS_PER_CASE
                snap_idx_end = snap_idx_init + SNAPSHOTS_PER_CASE - 1
                snap_init = f['field'][snap_idx_init]
                snap_end = f['field'][snap_idx_end]

                params = case_params[case_i]
                name = case_names[case_i]
                row_label = (f"[{lab.upper()}] {name}\n"
                             f"DR={params[0]}, VR={params[1]}, "
                             f"Re={params[2]}, Bo={params[3]}")

                # Col 0: Phase field at t=0 and t=end side by side
                ax = axes[row, 0]
                ax.imshow(snap_init, origin='lower', cmap='RdBu_r',
                          vmin=-0.5, vmax=0.5, aspect='auto')
                ax.set_title(f't=0.00', fontsize=10)
                ax.set_ylabel(row_label, fontsize=8, rotation=0, labelpad=120,
                              ha='right', va='center')

                # Col 1: Phase field at final snapshot
                ax = axes[row, 1]
                ax.imshow(snap_end, origin='lower', cmap='RdBu_r',
                          vmin=-0.5, vmax=0.5, aspect='auto')
                t_end = (SNAPSHOTS_PER_CASE - 1) * DT
                ax.set_title(f't={t_end:.2f}', fontsize=10)

                # Col 2: Circularity + velocity time series
                ax = axes[row, 2]
                ax.plot(time_circ, circularity[case_i], 'b-', linewidth=1.2,
                        label='Circularity')
                ax.set_ylabel('Circularity', color='b', fontsize=10)
                ax.set_ylim(0, 1.1)
                ax.tick_params(axis='y', labelcolor='b')
                ax.set_xlabel('Time', fontsize=10)
                ax.grid(True, alpha=0.3)

                if onset_arr[case_i] >= 0:
                    t_onset = onset_arr[case_i] * DT
                    ax.axvline(t_onset, color='green', ls='--', alpha=0.7,
                               label=f'onset={onset_arr[case_i]}')

                # Twin axis for rise velocity
                ax2 = ax.twinx()
                ax2.plot(time_vel, rise_velocity[case_i], 'r-', linewidth=0.8,
                         alpha=0.6, label='Rise vel.')
                ax2.set_ylabel('Rise velocity', color='r', fontsize=10)
                ax2.tick_params(axis='y', labelcolor='r')

                # Combined legend
                lines1, labels1 = ax.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                ax.legend(lines1 + lines2, labels1 + labels2,
                          fontsize=8, loc='upper left')

                row += 1

    fig.suptitle('Gallery: Representative Cases', fontsize=15, y=1.01)
    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'gallery_examples.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {path}")


def plot_marginals(case_params, n_steady_arr, labels):
    """Marginal distributions: n_steady conditioned on binned Re, Bo, DR, VR."""
    param_names = ['DR', 'VR', 'Re', 'Bo']
    n_bins = 5

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.ravel()

    for pi, (pname, ax) in enumerate(zip(param_names, axes)):
        vals = case_params[:, pi]
        bin_edges = np.quantile(vals, np.linspace(0, 1, n_bins + 1))
        bin_edges = np.unique(bin_edges)  # handle degenerate bins

        for bi in range(len(bin_edges) - 1):
            lo, hi = bin_edges[bi], bin_edges[bi + 1]
            if bi == len(bin_edges) - 2:
                mask = (vals >= lo) & (vals <= hi)
            else:
                mask = (vals >= lo) & (vals < hi)
            steady_mask = mask & (labels == 'steady')
            n_total = mask.sum()
            n_steady = steady_mask.sum()
            frac = n_steady / max(n_total, 1) * 100
            ax.bar(bi, frac, color='steelblue', edgecolor='white', alpha=0.8)
            ax.text(bi, frac + 1, f'{n_steady}/{n_total}', ha='center',
                    fontsize=8)

        bin_labels = [f'{bin_edges[i]:.0f}–{bin_edges[i+1]:.0f}'
                      for i in range(len(bin_edges) - 1)]
        ax.set_xticks(range(len(bin_labels)))
        ax.set_xticklabels(bin_labels, fontsize=9)
        ax.set_xlabel(pname, fontsize=12)
        ax.set_ylabel('% steady cases', fontsize=11)
        ax.set_ylim(0, 110)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_title(f'Steady fraction by {pname} bin', fontsize=12)

    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, 'marginal_steady_fraction.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


# ─── Main ───────────────────────────────────────────────────────────────────

def main():
    metrics_path = os.path.join(RESULTS_DIR, 'steady_state_metrics.npz')
    if not os.path.exists(metrics_path):
        print(f"ERROR: {metrics_path} not found.")
        print("Run batch_process_all_cases.py first.")
        return

    print(f"Loading: {metrics_path}")
    data = np.load(metrics_path, allow_pickle=True)
    circularity = data['circularity']     # (5000, 101)
    y_com = data['y_com']                 # (5000, 101)
    rise_velocity = data['rise_velocity'] # (5000, 100)
    n_fragments = data['n_fragments']     # (5000, 101)
    case_params = data['case_params']     # (5000, 4)
    case_names = data['case_names']       # (5000,)

    n_cases = circularity.shape[0]
    print(f"Loaded {n_cases} cases\n")

    # ── Phase 3: Classify ──
    print(f"Steady-state detection: window={WINDOW}, threshold={THRESHOLD*100}%")
    labels, onset_arr, n_steady_arr = classify_cases(
        circularity, n_fragments, WINDOW, THRESHOLD)

    for lab in ['steady', 'transient', 'breakup']:
        count = np.sum(labels == lab)
        print(f"  {lab:>10s}: {count:5d}  ({count/n_cases*100:5.1f}%)")

    steady_mask = labels == 'steady'
    if steady_mask.any():
        ss_n = n_steady_arr[steady_mask]
        print(f"\nSteady-state snapshot statistics:")
        print(f"  Mean N_steady:   {ss_n.mean():.1f}")
        print(f"  Median N_steady: {np.median(ss_n):.1f}")
        print(f"  Min N_steady:    {ss_n.min()}")
        print(f"  Max N_steady:    {ss_n.max()}")
        print(f"  Total steady-state snapshots across all cases: {ss_n.sum()}")

    # Save labels
    labels_path = os.path.join(RESULTS_DIR, 'steady_state_labels.npz')
    np.savez_compressed(
        labels_path,
        labels=labels,
        onset=onset_arr,
        n_steady=n_steady_arr,
        window=WINDOW,
        threshold=THRESHOLD,
        case_params=case_params,
        case_names=case_names
    )
    print(f"\nSaved labels: {labels_path}")

    # ── Phase 4: Plots ──
    print("\nGenerating plots ...")
    plot_histogram_n_steady(n_steady_arr, labels, THRESHOLD, WINDOW)
    plot_histogram_steady_onset(onset_arr, labels)
    plot_classification_pie(labels)
    plot_regime_map(case_params, n_steady_arr, labels)
    plot_marginals(case_params, n_steady_arr, labels)

    print("\nRunning threshold sensitivity sweep (may take a moment) ...")
    plot_threshold_sensitivity(circularity, n_fragments)

    print("\nGenerating example gallery ...")
    plot_gallery(circularity, rise_velocity, n_fragments, case_params,
                 case_names, labels, onset_arr, n_steady_arr)

    print("\n✓ All done. Check figures/ and results/ directories.")


if __name__ == '__main__':
    main()
