"""PNG renderers for the L3 CT slice and label overlay."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # non-interactive backend; safe for headless servers

import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt      # noqa: E402
import numpy as np                   # noqa: E402

# Discrete colormap for label IDs 0..4
_CMAP = mcolors.ListedColormap(["black", "cyan", "orange", "yellow", "red"])
_NORM = mcolors.BoundaryNorm([0, 0.5, 1.5, 2.5, 3.5, 4.5], 5)


def _window(arr, wc=40, ww=400):
    lo = wc - ww / 2
    hi = wc + ww / 2
    return np.clip((arr.astype(float) - lo) / (hi - lo), 0, 1)


def save_png_ct(ct_slice, path: str):
    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    ax.imshow(_window(ct_slice), cmap="gray", origin="upper")
    ax.set_title("L3 CT Slice (soft-tissue window)", fontsize=10)
    ax.axis("off")
    fig.tight_layout(pad=0)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_png_label(label_map, ct_slice, path: str):
    bg = _window(ct_slice, wc=-20, ww=400)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), dpi=150)
    axes[0].imshow(bg, cmap="gray", origin="upper")
    axes[0].set_title("CT slice", fontsize=9)
    axes[0].axis("off")

    axes[1].imshow(bg, cmap="gray", origin="upper")
    rgba = _CMAP(_NORM(label_map))
    rgba[label_map == 0, 3] = 0
    axes[1].imshow(rgba, origin="upper", alpha=0.65)
    axes[1].set_title("Label overlay", fontsize=9)
    axes[1].axis("off")

    from matplotlib.patches import Patch
    legend = [
        Patch(facecolor="cyan",   label="1 - subcutaneous fat"),
        Patch(facecolor="yellow", label="3 - IMAT (< -30 HU)"),
        Patch(facecolor="red",    label="4 - muscle (>= -30 HU)"),
    ]
    axes[1].legend(handles=legend, loc="lower right", fontsize=7, framealpha=0.7)
    fig.tight_layout(pad=0.5)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
