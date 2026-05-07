#!/usr/bin/env python3
"""
Selected publication-quality MD analysis plots from cpptraj outputs across replicates.

This script intentionally generates ONLY the selected outputs requested:

figures_md/hbond_avg_top_contacts_frac_replicate_analysis_1.png
figures_md/hbond_avg_top_contacts_frac_replicate_analysis_2.png
figures_md/hbond_avg_top_contacts_frac_replicate_analysis_3.png
figures_md/hbonds_distribution_boxplot.png
figures_md/hbonds_mean.png
figures_md/hbonds_mean_std.png
figures_md/hbonds_replicates.png
figures_md/hbonds_total_replicate_analysis_1.png
figures_md/hbonds_total_replicate_analysis_2.png
figures_md/hbonds_total_replicate_analysis_3.png
figures_md/md_replicates_summary.csv
figures_md/rg_protein_mean_std.png
figures_md/rg_protein_replicates.png
figures_md/rmsd_ligand_mean_std.png
figures_md/rmsd_ligand_replicates.png
figures_md/rmsf_ca_mean_std.png
figures_md/rmsf_ca_replicates.png

Additional outputs requested:
figures_md/hbond_avg_top_contacts_frac_subplot_3x1.png
figures_md/md_metrics_mean_std_subplot_2x2.png
figures_md/hbond_lig_prot_avg_detailed.xlsx

Expected directory structure:
    replicate_analysis_1/
        rmsd_ligand.dat
        rmsf_ca.dat
        rg_protein.dat
        hbond_lig_prot.dat
        hbond_lig_prot_avg.dat

    replicate_analysis_2/
        ...

    replicate_analysis_3/
        ...

Also accepts input folders named:
    analise_replica_1/
    analysis_replica_1/

Even if the input folder is named analise_replica_1, the output files will use
replicate_analysis_1 in the filename, matching the requested output names.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# If each replicate trajectory should be mapped from 0 to TOTAL_TIME_NS.
# Example: if each production trajectory has 100 ns, keep 100.0.
# If your trajectory has 1 ns, change to 1.0.
TOTAL_TIME_NS = 100.0

# If you know the time spacing between saved frames, set it here.
# Example: TIME_PER_FRAME_NS = 0.1
# If not None, this has priority over TOTAL_TIME_NS.
TIME_PER_FRAME_NS = None

OUTPUT_DIR = "figures_md"

# Only these replicate indices will be used.
REPLICATE_INDICES = (1, 2, 3)

# Accepted input directory patterns.
REPLICATE_DIR_CANDIDATES = (
    "replicate_analysis_{i}",
    "analise_replica_{i}",
    "analysis_replica_{i}",
)

# Input filenames expected inside each replicate directory.
RMSD_FILE = "rmsd_ligand.dat"
RMSF_FILE = "rmsf_ca.dat"
RG_FILE = "rg_protein.dat"
HBOND_TIMESERIES_FILE = "hbond_lig_prot.dat"
HBOND_AVG_FILE = "hbond_lig_prot_avg.dat"

# Number of H-bond contacts shown in individual and subplot barplots.
HBOND_TOP_N = 25

# If the 3x1 subplot becomes too large, reduce this to 10 or 15.
HBOND_TOP_N_SUBPLOT = HBOND_TOP_N

# Output names for the two additional subplot figures and the detailed Excel file.
HBOND_AVG_SUBPLOT_NAME = "hbond_avg_top_contacts_frac_subplot_3x1"
MD_METRICS_SUBPLOT_NAME = "md_metrics_mean_std_subplot_2x2"
HBOND_AVG_DETAILED_XLSX = "hbond_lig_prot_avg_detailed.xlsx"


# =============================================================================
# MATPLOTLIB CONFIGURATION
# =============================================================================

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 1.2,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
})

# Keep warnings visible by default. Version-specific Matplotlib compatibility is
# handled locally where needed.


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class Replicate:
    index: int
    path: Path

    @property
    def suffix(self) -> str:
        return f"replicate_analysis_{self.index}"

    @property
    def label(self) -> str:
        return f"Replicate {self.index}"


# =============================================================================
# GENERAL UTILITIES
# =============================================================================

def nature_axes(ax: plt.Axes) -> None:
    """Apply a clean publication-style axis format."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4, width=1)


def ensure_output_dir(path: str | Path = OUTPUT_DIR) -> Path:
    outdir = Path(path)
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def save_png(fig: plt.Figure, outdir: Path, basename: str) -> Path:
    filepath = outdir / f"{basename}.png"
    fig.savefig(filepath)
    plt.close(fig)
    print(f"[OK] Saved: {filepath}")
    return filepath


def find_replicates(base: str | Path = ".") -> list[Replicate]:
    """Find replicate folders and map them to canonical output suffixes."""
    base = Path(base)
    reps: list[Replicate] = []

    for i in REPLICATE_INDICES:
        found_path = None
        for pattern in REPLICATE_DIR_CANDIDATES:
            candidate = base / pattern.format(i=i)
            if candidate.is_dir():
                found_path = candidate
                break

        if found_path is None:
            print(f"[WARN] Replicate {i} directory not found.")
        else:
            reps.append(Replicate(index=i, path=found_path))

    return reps


def build_time_axis(n_points: int, frames: np.ndarray | None = None) -> np.ndarray:
    """Build time axis in ns."""
    if n_points <= 0:
        return np.array([])

    if TIME_PER_FRAME_NS is not None:
        if frames is not None and len(frames) == n_points:
            frames = np.asarray(frames, dtype=float)
            return (frames - frames[0]) * float(TIME_PER_FRAME_NS)
        return np.arange(n_points, dtype=float) * float(TIME_PER_FRAME_NS)

    if TOTAL_TIME_NS is not None:
        return np.linspace(0.0, float(TOTAL_TIME_NS), n_points)

    raise ValueError("Define TIME_PER_FRAME_NS or TOTAL_TIME_NS.")


def load_numeric_dat(filepath: Path) -> pd.DataFrame:
    """Load whitespace-separated numeric .dat file while ignoring comments."""
    rows: list[list[float]] = []

    with open(filepath, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def truncate_series_to_min_length(series: list[np.ndarray]) -> np.ndarray:
    """Trim all time series to the shortest length and stack into a matrix."""
    min_len = min(len(x) for x in series)
    return np.array([x[:min_len] for x in series], dtype=float)


# =============================================================================
# CPPTRAJ FILE LOADERS
# =============================================================================

def load_rmsd_file(filepath: Path) -> pd.DataFrame | None:
    df = load_numeric_dat(filepath)
    if df.empty or df.shape[1] < 2:
        return None

    out = df.iloc[:, :2].copy()
    out.columns = ["Frame", "RMSD"]
    out = out.apply(pd.to_numeric, errors="coerce").dropna()
    return out


def load_rmsf_file(filepath: Path) -> pd.DataFrame | None:
    df = load_numeric_dat(filepath)
    if df.empty or df.shape[1] < 2:
        return None

    out = df.iloc[:, :2].copy()
    out.columns = ["Residue", "RMSF"]
    out = out.apply(pd.to_numeric, errors="coerce").dropna()
    out["Residue"] = out["Residue"].astype(int)
    return out


def load_rg_file(filepath: Path) -> pd.DataFrame | None:
    df = load_numeric_dat(filepath)
    if df.empty or df.shape[1] < 2:
        return None

    out = df.iloc[:, :2].copy()
    out.columns = ["Frame", "Rg"]
    out = out.apply(pd.to_numeric, errors="coerce").dropna()
    return out


def load_hbond_timeseries(filepath: Path) -> pd.DataFrame | None:
    """
    Load hbond_lig_prot.dat.

    Expected main format:
        #Frame      HBALL[UU]
             1           22
             2           22

    If a file has 4 or more numeric columns, the fourth column is treated as
    total H-bonds for compatibility with direction-based outputs.
    """
    df = load_numeric_dat(filepath)
    if df.empty:
        return None

    if df.shape[1] >= 4:
        out = pd.DataFrame({
            "Frame": df.iloc[:, 0],
            "HBonds": df.iloc[:, 3],
        })
    elif df.shape[1] >= 2:
        out = pd.DataFrame({
            "Frame": df.iloc[:, 0],
            "HBonds": df.iloc[:, 1],
        })
    else:
        out = pd.DataFrame({
            "Frame": np.arange(1, len(df) + 1),
            "HBonds": df.iloc[:, 0],
        })

    out = out.apply(pd.to_numeric, errors="coerce").dropna(subset=["Frame", "HBonds"])
    return out


# =============================================================================
# HBOND AVG PARSING FOR ARTICLE-READY TABLE
# =============================================================================

def split_cpptraj_residue_atom(token: str) -> dict[str, str | int | None]:
    """
    Parse cpptraj-like atom tokens.

    Example:
        ASP_176@OD1_2680

    Returns:
        residue_full = ASP_176
        residue_name = ASP
        residue_id = 176
        atom_name = OD1
        atom_index = 2680
    """
    token = str(token).strip()

    if "@" in token:
        residue_part, atom_part = token.split("@", 1)
    else:
        residue_part, atom_part = token, ""

    residue_name = residue_part
    residue_id: int | None = None
    if "_" in residue_part:
        left, right = residue_part.rsplit("_", 1)
        if right.lstrip("-").isdigit():
            residue_name = left
            residue_id = int(right)

    atom_name = atom_part
    atom_index: int | None = None
    if "_" in atom_part:
        left, right = atom_part.rsplit("_", 1)
        if right.lstrip("-").isdigit():
            atom_name = left
            atom_index = int(right)

    return {
        "raw": token,
        "residue_full": residue_part,
        "residue_name": residue_name,
        "residue_id": residue_id,
        "atom_name": atom_name,
        "atom_index": atom_index,
    }


def load_hbond_avg(filepath: Path, replicate: Replicate | None = None) -> pd.DataFrame:
    """
    Load cpptraj hbond_lig_prot_avg.dat and split acceptor/donor fields.

    Expected columns:
        #Acceptor DonorH Donor Frames Frac AvgDist AvgAng

    Example:
        ASP_176@OD1_2680 0CU_630@H4O_9667 0CU_630@O4_9666 1000 1.0000 2.6293 163.8671
    """
    rows: list[dict[str, object]] = []

    with open(filepath, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 7:
                continue

            acceptor_token, donor_h_token, donor_token = parts[0], parts[1], parts[2]

            try:
                frames = int(float(parts[3]))
                frac = float(parts[4])
                avg_dist = float(parts[5])
                avg_ang = float(parts[6])
            except ValueError:
                continue

            acc = split_cpptraj_residue_atom(acceptor_token)
            don_h = split_cpptraj_residue_atom(donor_h_token)
            don = split_cpptraj_residue_atom(donor_token)

            replicate_index = replicate.index if replicate is not None else None
            replicate_label = replicate.label if replicate is not None else "Current directory"
            replicate_suffix = replicate.suffix if replicate is not None else "current_directory"

            interaction = f"{acc['residue_full']}@{acc['atom_name']} ← {don['residue_full']}@{don['atom_name']}"

            rows.append({
                "Replicate_index": replicate_index,
                "Replicate": replicate_label,
                "Replicate_suffix": replicate_suffix,
                "Interaction": interaction,

                "Acceptor_raw": acc["raw"],
                "Acceptor_residue": acc["residue_full"],
                "Acceptor_residue_name": acc["residue_name"],
                "Acceptor_residue_id": acc["residue_id"],
                "Acceptor_atom": acc["atom_name"],
                "Acceptor_atom_index": acc["atom_index"],

                "DonorH_raw": don_h["raw"],
                "DonorH_residue": don_h["residue_full"],
                "DonorH_residue_name": don_h["residue_name"],
                "DonorH_residue_id": don_h["residue_id"],
                "DonorH_atom": don_h["atom_name"],
                "DonorH_atom_index": don_h["atom_index"],

                "Donor_raw": don["raw"],
                "Donor_residue": don["residue_full"],
                "Donor_residue_name": don["residue_name"],
                "Donor_residue_id": don["residue_id"],
                "Donor_atom": don["atom_name"],
                "Donor_atom_index": don["atom_index"],

                "Frames": frames,
                "Frac": frac,
                "AvgDist_A": avg_dist,
                "AvgAng_deg": avg_ang,
            })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.sort_values(
        ["Replicate_index", "Frac", "Frames"],
        ascending=[True, False, False],
        na_position="last",
    ).reset_index(drop=True)

    # Rank contacts within each replicate.
    if "Replicate" in df.columns:
        df.insert(0, "Rank", df.groupby("Replicate").cumcount() + 1)

    return df


# =============================================================================
# RMSD PLOTS
# =============================================================================

def collect_rmsd_data(replicates: list[Replicate]) -> dict[str, object] | None:
    series: list[np.ndarray] = []
    labels: list[str] = []

    for rep in replicates:
        filepath = rep.path / RMSD_FILE
        if not filepath.exists():
            print(f"[WARN] Missing {filepath}")
            continue

        df = load_rmsd_file(filepath)
        if df is None or df.empty:
            print(f"[WARN] Could not parse {filepath}")
            continue

        series.append(df["RMSD"].to_numpy(dtype=float))
        labels.append(rep.label)

    if not series:
        return None

    matrix = truncate_series_to_min_length(series)
    time_ns = build_time_axis(matrix.shape[1])
    return {
        "time_ns": time_ns,
        "matrix": matrix,
        "labels": labels,
        "mean": matrix.mean(axis=0),
        "std": matrix.std(axis=0, ddof=1) if matrix.shape[0] > 1 else np.zeros(matrix.shape[1]),
    }


def plot_rmsd_replicates(data: dict[str, object], outdir: Path) -> None:
    time_ns = data["time_ns"]
    matrix = data["matrix"]
    labels = data["labels"]

    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    for i, label in enumerate(labels):
        ax.plot(time_ns, matrix[i], lw=1.0, alpha=0.85, label=label)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Ligand RMSD (Å)")
    ax.set_title("Ligand RMSD per replicate")
    ax.legend(frameon=False)
    nature_axes(ax)
    save_png(fig, outdir, "rmsd_ligand_replicates")


def plot_rmsd_mean_std(data: dict[str, object], outdir: Path) -> None:
    time_ns = data["time_ns"]
    mean = data["mean"]
    std = data["std"]

    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.plot(time_ns, mean, lw=1.8, label="Mean")
    ax.fill_between(time_ns, mean - std, mean + std, alpha=0.25, label="± SD")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Ligand RMSD (Å)")
    ax.set_title("Ligand RMSD")
    ax.legend(frameon=False)
    nature_axes(ax)
    save_png(fig, outdir, "rmsd_ligand_mean_std")


# =============================================================================
# RG PLOTS
# =============================================================================

def collect_rg_data(replicates: list[Replicate]) -> dict[str, object] | None:
    series: list[np.ndarray] = []
    labels: list[str] = []

    for rep in replicates:
        filepath = rep.path / RG_FILE
        if not filepath.exists():
            print(f"[WARN] Missing {filepath}")
            continue

        df = load_rg_file(filepath)
        if df is None or df.empty:
            print(f"[WARN] Could not parse {filepath}")
            continue

        series.append(df["Rg"].to_numpy(dtype=float))
        labels.append(rep.label)

    if not series:
        return None

    matrix = truncate_series_to_min_length(series)
    time_ns = build_time_axis(matrix.shape[1])
    return {
        "time_ns": time_ns,
        "matrix": matrix,
        "labels": labels,
        "mean": matrix.mean(axis=0),
        "std": matrix.std(axis=0, ddof=1) if matrix.shape[0] > 1 else np.zeros(matrix.shape[1]),
    }


def plot_rg_replicates(data: dict[str, object], outdir: Path) -> None:
    time_ns = data["time_ns"]
    matrix = data["matrix"]
    labels = data["labels"]

    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    for i, label in enumerate(labels):
        ax.plot(time_ns, matrix[i], lw=1.0, alpha=0.85, label=label)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Radius of gyration (Å)")
    ax.set_title("Protein Rg per replicate")
    ax.legend(frameon=False)
    nature_axes(ax)
    save_png(fig, outdir, "rg_protein_replicates")


def plot_rg_mean_std(data: dict[str, object], outdir: Path) -> None:
    time_ns = data["time_ns"]
    mean = data["mean"]
    std = data["std"]

    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.plot(time_ns, mean, lw=1.8, label="Mean")
    ax.fill_between(time_ns, mean - std, mean + std, alpha=0.25, label="± SD")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Radius of gyration (Å)")
    ax.set_title("Protein radius of gyration")
    ax.legend(frameon=False)
    nature_axes(ax)
    save_png(fig, outdir, "rg_protein_mean_std")


# =============================================================================
# RMSF PLOTS
# =============================================================================

def collect_rmsf_data(replicates: list[Replicate]) -> dict[str, object] | None:
    dfs: list[pd.DataFrame] = []
    labels: list[str] = []

    for rep in replicates:
        filepath = rep.path / RMSF_FILE
        if not filepath.exists():
            print(f"[WARN] Missing {filepath}")
            continue

        df = load_rmsf_file(filepath)
        if df is None or df.empty:
            print(f"[WARN] Could not parse {filepath}")
            continue

        dfs.append(df)
        labels.append(rep.label)

    if not dfs:
        return None

    merged = dfs[0].rename(columns={"RMSF": labels[0]})
    for df, label in zip(dfs[1:], labels[1:]):
        merged = merged.merge(
            df.rename(columns={"RMSF": label}),
            on="Residue",
            how="inner",
        )

    value_cols = [c for c in merged.columns if c != "Residue"]
    matrix = merged[value_cols].to_numpy(dtype=float).T

    return {
        "residues": merged["Residue"].to_numpy(dtype=int),
        "matrix": matrix,
        "labels": value_cols,
        "mean": matrix.mean(axis=0),
        "std": matrix.std(axis=0, ddof=1) if matrix.shape[0] > 1 else np.zeros(matrix.shape[1]),
    }


def plot_rmsf_replicates(data: dict[str, object], outdir: Path) -> None:
    residues = data["residues"]
    matrix = data["matrix"]
    labels = data["labels"]

    fig, ax = plt.subplots(figsize=(5.8, 3.4))
    for i, label in enumerate(labels):
        ax.plot(residues, matrix[i], lw=1.0, alpha=0.85, label=label)

    ax.set_xlabel("Residue index")
    ax.set_ylabel("Cα RMSF (Å)")
    ax.set_title("Protein Cα RMSF per replicate")
    ax.legend(frameon=False, ncol=2)
    nature_axes(ax)
    save_png(fig, outdir, "rmsf_ca_replicates")


def plot_rmsf_mean_std(data: dict[str, object], outdir: Path) -> None:
    residues = data["residues"]
    mean = data["mean"]
    std = data["std"]

    fig, ax = plt.subplots(figsize=(5.8, 3.4))
    ax.plot(residues, mean, lw=1.8, label="Mean")
    ax.fill_between(residues, mean - std, mean + std, alpha=0.25, label="± SD")

    ax.set_xlabel("Residue index")
    ax.set_ylabel("Cα RMSF (Å)")
    ax.set_title("Protein Cα RMSF")
    ax.legend(frameon=False)
    nature_axes(ax)
    save_png(fig, outdir, "rmsf_ca_mean_std")


# =============================================================================
# HBOND TIME-SERIES PLOTS
# =============================================================================

def collect_hbond_timeseries_data(replicates: list[Replicate]) -> dict[str, object] | None:
    series: list[np.ndarray] = []
    labels: list[str] = []
    per_replicate: dict[int, pd.DataFrame] = {}

    for rep in replicates:
        filepath = rep.path / HBOND_TIMESERIES_FILE
        if not filepath.exists():
            print(f"[WARN] Missing {filepath}")
            continue

        df = load_hbond_timeseries(filepath)
        if df is None or df.empty:
            print(f"[WARN] Could not parse {filepath}")
            continue

        per_replicate[rep.index] = df
        series.append(df["HBonds"].to_numpy(dtype=float))
        labels.append(rep.label)

    if not series:
        return None

    matrix = truncate_series_to_min_length(series)
    time_ns = build_time_axis(matrix.shape[1])

    return {
        "time_ns": time_ns,
        "matrix": matrix,
        "labels": labels,
        "per_replicate": per_replicate,
        "mean": matrix.mean(axis=0),
        "std": matrix.std(axis=0, ddof=1) if matrix.shape[0] > 1 else np.zeros(matrix.shape[1]),
    }


def plot_hbond_total_per_replicate(replicates: list[Replicate], hbond_data: dict[str, object], outdir: Path) -> None:
    per_replicate = hbond_data["per_replicate"]

    for rep in replicates:
        if rep.index not in per_replicate:
            continue

        df = per_replicate[rep.index]
        time_ns = build_time_axis(len(df), frames=df["Frame"].to_numpy())
        hbonds = df["HBonds"].to_numpy(dtype=float)

        fig, ax = plt.subplots(figsize=(5.0, 3.5))
        ax.plot(time_ns, hbonds, lw=1.1)

        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Number of H-bonds")
        ax.set_title(f"Ligand–protein H-bonds - {rep.label}")
        nature_axes(ax)
        save_png(fig, outdir, f"hbonds_total_{rep.suffix}")


def plot_hbonds_replicates(hbond_data: dict[str, object], outdir: Path) -> None:
    time_ns = hbond_data["time_ns"]
    matrix = hbond_data["matrix"]
    labels = hbond_data["labels"]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for i, label in enumerate(labels):
        ax.plot(time_ns, matrix[i], lw=1.0, alpha=0.8, label=label)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title("Ligand–protein H-bonds per replicate")
    ax.legend(frameon=False)
    nature_axes(ax)
    save_png(fig, outdir, "hbonds_replicates")


def plot_hbonds_mean_std(hbond_data: dict[str, object], outdir: Path) -> None:
    time_ns = hbond_data["time_ns"]
    mean = hbond_data["mean"]
    std = hbond_data["std"]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(time_ns, mean, lw=1.8, label="Mean")
    ax.fill_between(time_ns, mean - std, mean + std, alpha=0.25, label="± SD")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title("Ligand–protein H-bonds")
    ax.legend(frameon=False)
    nature_axes(ax)
    save_png(fig, outdir, "hbonds_mean_std")


def plot_hbonds_mean(hbond_data: dict[str, object], outdir: Path) -> None:
    time_ns = hbond_data["time_ns"]
    mean = hbond_data["mean"]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(time_ns, mean, lw=1.8)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title("Mean ligand–protein H-bonds")
    nature_axes(ax)
    save_png(fig, outdir, "hbonds_mean")


def plot_hbonds_distribution(hbond_data: dict[str, object], outdir: Path) -> None:
    matrix = hbond_data["matrix"]
    labels = hbond_data["labels"]

    fig, ax = plt.subplots(figsize=(4.8, 3.4))

    boxplot_values = [matrix[i, :] for i in range(matrix.shape[0])]
    try:
        ax.boxplot(
            boxplot_values,
            tick_labels=labels,
            showfliers=False,
            patch_artist=False,
        )
    except TypeError:
        # Compatibility with older Matplotlib versions.
        ax.boxplot(
            boxplot_values,
            labels=labels,
            showfliers=False,
            patch_artist=False,
        )

    ax.set_xlabel("Replicate")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title("H-bond count distribution")
    nature_axes(ax)
    save_png(fig, outdir, "hbonds_distribution_boxplot")


# =============================================================================
# HBOND AVG PLOTS AND XLSX
# =============================================================================

def plot_hbond_avg_top_contacts(df_avg: pd.DataFrame, rep: Replicate, outdir: Path, top_n: int = HBOND_TOP_N) -> None:
    if df_avg.empty:
        return

    top = df_avg.head(min(top_n, len(df_avg))).copy()
    top_plot = top.iloc[::-1].copy()

    fig_height = max(4.2, 0.30 * len(top_plot) + 1.4)
    fig, ax = plt.subplots(figsize=(8.6, fig_height))

    y = np.arange(len(top_plot))
    ax.barh(y, top_plot["Frac"])
    ax.set_yticks(y)
    ax.set_yticklabels(top_plot["Interaction"], fontsize=8)

    ax.set_xlabel("Occupancy fraction")
    ax.set_ylabel("H-bond interaction")
    ax.set_title(f"Top ligand–protein H-bonds by occupancy - {rep.label}")
    ax.set_xlim(0, max(1.0, float(top_plot["Frac"].max()) * 1.05))
    nature_axes(ax)

    save_png(fig, outdir, f"hbond_avg_top_contacts_frac_{rep.suffix}")


def collect_hbond_avg_data(replicates: list[Replicate]) -> dict[int, pd.DataFrame]:
    avg_data: dict[int, pd.DataFrame] = {}

    for rep in replicates:
        filepath = rep.path / HBOND_AVG_FILE
        if not filepath.exists():
            print(f"[WARN] Missing {filepath}")
            continue

        df = load_hbond_avg(filepath, replicate=rep)
        if df.empty:
            print(f"[WARN] Could not parse {filepath}")
            continue

        avg_data[rep.index] = df

    return avg_data


def export_hbond_avg_detailed_xlsx(avg_data: dict[int, pd.DataFrame], outdir: Path) -> None:
    if not avg_data:
        print("[WARN] No hbond_lig_prot_avg.dat data available for XLSX export.")
        return

    xlsx_path = outdir / HBOND_AVG_DETAILED_XLSX

    preferred_cols = [
        "Rank",
        "Replicate",
        "Interaction",
        "Acceptor_residue",
        "Acceptor_residue_name",
        "Acceptor_residue_id",
        "Acceptor_atom",
        "Acceptor_atom_index",
        "DonorH_residue",
        "DonorH_residue_name",
        "DonorH_residue_id",
        "DonorH_atom",
        "DonorH_atom_index",
        "Donor_residue",
        "Donor_residue_name",
        "Donor_residue_id",
        "Donor_atom",
        "Donor_atom_index",
        "Frames",
        "Frac",
        "AvgDist_A",
        "AvgAng_deg",
        "Acceptor_raw",
        "DonorH_raw",
        "Donor_raw",
    ]

    combined = pd.concat(avg_data.values(), ignore_index=True)
    preferred_cols = [c for c in preferred_cols if c in combined.columns]
    combined = combined[preferred_cols]

    try:
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            combined.to_excel(writer, sheet_name="all_replicates", index=False)

            for rep_index, df in avg_data.items():
                sheet_name = f"replicate_{rep_index}"
                df_out = df[[c for c in preferred_cols if c in df.columns]]
                df_out.to_excel(writer, sheet_name=sheet_name, index=False)

            # Simple top contacts sheet for quick article screening.
            top_contacts = (
                combined.sort_values(["Replicate", "Frac", "Frames"], ascending=[True, False, False])
                .groupby("Replicate", as_index=False)
                .head(10)
            )
            top_contacts.to_excel(writer, sheet_name="top10_each_replicate", index=False)

            # Formatting.
            workbook = writer.book
            for worksheet in workbook.worksheets:
                worksheet.freeze_panes = "A2"
                for col_cells in worksheet.columns:
                    max_len = 0
                    column_letter = col_cells[0].column_letter
                    for cell in col_cells:
                        value = "" if cell.value is None else str(cell.value)
                        max_len = max(max_len, len(value))
                    worksheet.column_dimensions[column_letter].width = min(max(max_len + 2, 10), 42)

    except ImportError:
        print("[ERROR] openpyxl is required for XLSX export.")
        print("        Install it with: pip install openpyxl")
        return

    print(f"[OK] Saved: {xlsx_path}")


def plot_hbond_avg_subplot_3x1(
    replicates: list[Replicate],
    avg_data: dict[int, pd.DataFrame],
    outdir: Path,
    top_n: int = HBOND_TOP_N_SUBPLOT,
) -> None:
    valid_reps = [rep for rep in replicates if rep.index in avg_data and not avg_data[rep.index].empty]
    if not valid_reps:
        print("[WARN] No hbond_lig_prot_avg.dat data available for 3x1 subplot.")
        return

    nrows = len(valid_reps)
    fig_height = max(8.0, 0.33 * top_n * nrows + 1.5)
    fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(9.2, fig_height))

    if nrows == 1:
        axes = [axes]

    for ax, rep in zip(axes, valid_reps):
        df = avg_data[rep.index]
        top = df.head(min(top_n, len(df))).iloc[::-1].copy()

        y = np.arange(len(top))
        ax.barh(y, top["Frac"])
        ax.set_yticks(y)
        ax.set_yticklabels(top["Interaction"], fontsize=7)
        ax.set_xlim(0, max(1.0, float(top["Frac"].max()) * 1.05))
        ax.set_xlabel("Occupancy fraction")
        ax.set_title(f"{rep.label}")
        nature_axes(ax)

    fig.suptitle("Top ligand–protein H-bonds by occupancy", y=0.995, fontsize=13)
    fig.tight_layout()
    save_png(fig, outdir, HBOND_AVG_SUBPLOT_NAME)


# =============================================================================
# MD METRICS 2x2 SUBPLOT
# =============================================================================

def plot_mean_std_on_axis(
    ax: plt.Axes,
    x: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    ax.plot(x, mean, lw=1.6, label="Mean")
    ax.fill_between(x, mean - std, mean + std, alpha=0.25, label="± SD")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False)
    nature_axes(ax)


def plot_md_metrics_subplot_2x2(
    hbond_data: dict[str, object] | None,
    rmsd_data: dict[str, object] | None,
    rg_data: dict[str, object] | None,
    rmsf_data: dict[str, object] | None,
    outdir: Path,
) -> None:
    datasets = [hbond_data, rmsd_data, rg_data, rmsf_data]
    if all(d is None for d in datasets):
        print("[WARN] No data available for 2x2 mean ± SD subplot.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4))
    axes = axes.ravel()

    if hbond_data is not None:
        plot_mean_std_on_axis(
            axes[0],
            hbond_data["time_ns"],
            hbond_data["mean"],
            hbond_data["std"],
            "Time (ns)",
            "Number of H-bonds",
            "Ligand–protein H-bonds",
        )
    else:
        axes[0].text(0.5, 0.5, "H-bond data not found", ha="center", va="center")
        axes[0].axis("off")

    if rmsd_data is not None:
        plot_mean_std_on_axis(
            axes[1],
            rmsd_data["time_ns"],
            rmsd_data["mean"],
            rmsd_data["std"],
            "Time (ns)",
            "Ligand RMSD (Å)",
            "Ligand RMSD",
        )
    else:
        axes[1].text(0.5, 0.5, "RMSD data not found", ha="center", va="center")
        axes[1].axis("off")

    if rg_data is not None:
        plot_mean_std_on_axis(
            axes[2],
            rg_data["time_ns"],
            rg_data["mean"],
            rg_data["std"],
            "Time (ns)",
            "Radius of gyration (Å)",
            "Protein radius of gyration",
        )
    else:
        axes[2].text(0.5, 0.5, "Rg data not found", ha="center", va="center")
        axes[2].axis("off")

    if rmsf_data is not None:
        plot_mean_std_on_axis(
            axes[3],
            rmsf_data["residues"],
            rmsf_data["mean"],
            rmsf_data["std"],
            "Residue index",
            "Cα RMSF (Å)",
            "Protein Cα RMSF",
        )
    else:
        axes[3].text(0.5, 0.5, "RMSF data not found", ha="center", va="center")
        axes[3].axis("off")

    fig.tight_layout()
    save_png(fig, outdir, MD_METRICS_SUBPLOT_NAME)


# =============================================================================
# SUMMARY TABLE
# =============================================================================

def generate_summary_table(
    replicates: list[Replicate],
    avg_data: dict[int, pd.DataFrame],
    outdir: Path,
) -> None:
    rows: list[dict[str, object]] = []

    for rep in replicates:
        row: dict[str, object] = {
            "Replicate_index": rep.index,
            "Replicate": rep.label,
            "Directory": str(rep.path),
        }

        filepath = rep.path / RMSD_FILE
        if filepath.exists():
            df = load_rmsd_file(filepath)
            if df is not None and not df.empty:
                row["RMSD_mean"] = df["RMSD"].mean()
                row["RMSD_sd"] = df["RMSD"].std()
                row["RMSD_min"] = df["RMSD"].min()
                row["RMSD_max"] = df["RMSD"].max()

        filepath = rep.path / RMSF_FILE
        if filepath.exists():
            df = load_rmsf_file(filepath)
            if df is not None and not df.empty:
                row["RMSF_mean"] = df["RMSF"].mean()
                row["RMSF_sd"] = df["RMSF"].std()
                row["RMSF_min"] = df["RMSF"].min()
                row["RMSF_max"] = df["RMSF"].max()

        filepath = rep.path / RG_FILE
        if filepath.exists():
            df = load_rg_file(filepath)
            if df is not None and not df.empty:
                row["Rg_mean"] = df["Rg"].mean()
                row["Rg_sd"] = df["Rg"].std()
                row["Rg_min"] = df["Rg"].min()
                row["Rg_max"] = df["Rg"].max()

        filepath = rep.path / HBOND_TIMESERIES_FILE
        if filepath.exists():
            df = load_hbond_timeseries(filepath)
            if df is not None and not df.empty:
                row["HBonds_mean"] = df["HBonds"].mean()
                row["HBonds_sd"] = df["HBonds"].std()
                row["HBonds_min"] = df["HBonds"].min()
                row["HBonds_max"] = df["HBonds"].max()

        if rep.index in avg_data and not avg_data[rep.index].empty:
            df_avg = avg_data[rep.index]
            top = df_avg.iloc[0]
            row["HBond_avg_contacts_n"] = len(df_avg)
            row["HBond_avg_top_interaction"] = top["Interaction"]
            row["HBond_avg_top_acceptor_residue"] = top["Acceptor_residue"]
            row["HBond_avg_top_acceptor_atom"] = top["Acceptor_atom"]
            row["HBond_avg_top_donor_residue"] = top["Donor_residue"]
            row["HBond_avg_top_donor_atom"] = top["Donor_atom"]
            row["HBond_avg_top_frames"] = top["Frames"]
            row["HBond_avg_top_frac"] = top["Frac"]
            row["HBond_avg_top_dist_A"] = top["AvgDist_A"]
            row["HBond_avg_top_angle_deg"] = top["AvgAng_deg"]

        rows.append(row)

    if not rows:
        print("[WARN] No data available for summary table.")
        return

    summary = pd.DataFrame(rows)
    filepath = outdir / "md_replicates_summary.csv"
    summary.to_csv(filepath, index=False)
    print(f"[OK] Saved: {filepath}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    base = Path(".")
    outdir = ensure_output_dir(OUTPUT_DIR)

    replicates = find_replicates(base)
    if not replicates:
        print("[ERROR] No replicate directories found.")
        print("Expected one of these patterns:")
        for i in REPLICATE_INDICES:
            for pattern in REPLICATE_DIR_CANDIDATES:
                print(f"  - {pattern.format(i=i)}")
        sys.exit(1)

    print("[INFO] Replicates used:")
    for rep in replicates:
        print(f"  - {rep.label}: {rep.path}")

    # Collect data once, then reuse for individual plots, subplots, and tables.
    rmsd_data = collect_rmsd_data(replicates)
    rmsf_data = collect_rmsf_data(replicates)
    rg_data = collect_rg_data(replicates)
    hbond_data = collect_hbond_timeseries_data(replicates)
    hbond_avg_data = collect_hbond_avg_data(replicates)

    # RMSD outputs.
    if rmsd_data is not None:
        plot_rmsd_replicates(rmsd_data, outdir)
        plot_rmsd_mean_std(rmsd_data, outdir)
    else:
        print(f"[WARN] No valid {RMSD_FILE} data found.")

    # RMSF outputs.
    if rmsf_data is not None:
        plot_rmsf_replicates(rmsf_data, outdir)
        plot_rmsf_mean_std(rmsf_data, outdir)
    else:
        print(f"[WARN] No valid {RMSF_FILE} data found.")

    # Rg outputs.
    if rg_data is not None:
        plot_rg_replicates(rg_data, outdir)
        plot_rg_mean_std(rg_data, outdir)
    else:
        print(f"[WARN] No valid {RG_FILE} data found.")

    # H-bond time-series outputs.
    if hbond_data is not None:
        plot_hbond_total_per_replicate(replicates, hbond_data, outdir)
        plot_hbonds_replicates(hbond_data, outdir)
        plot_hbonds_mean_std(hbond_data, outdir)
        plot_hbonds_mean(hbond_data, outdir)
        plot_hbonds_distribution(hbond_data, outdir)
    else:
        print(f"[WARN] No valid {HBOND_TIMESERIES_FILE} data found.")

    # H-bond avg/occupancy individual plots.
    for rep in replicates:
        if rep.index in hbond_avg_data:
            plot_hbond_avg_top_contacts(hbond_avg_data[rep.index], rep, outdir, top_n=HBOND_TOP_N)

    # H-bond avg detailed XLSX.
    export_hbond_avg_detailed_xlsx(hbond_avg_data, outdir)

    # Summary CSV.
    generate_summary_table(replicates, hbond_avg_data, outdir)

    # Additional subplot figures requested.
    plot_hbond_avg_subplot_3x1(replicates, hbond_avg_data, outdir, top_n=HBOND_TOP_N_SUBPLOT)
    plot_md_metrics_subplot_2x2(hbond_data, rmsd_data, rg_data, rmsf_data, outdir)

    print(f"[DONE] Selected outputs generated in: {outdir.resolve()}")


if __name__ == "__main__":
    main()
