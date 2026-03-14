#!/usr/bin/env python3

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Total simulated time of each replicate in ns
# Example: 100.0 for a 100 ns production run
TOTAL_TIME_NS = 100.0

# If you prefer to define time per saved frame instead of total time,
# set this to a float and TOTAL_TIME_NS to None.
# Example: TIME_PER_FRAME_NS = 0.1
TIME_PER_FRAME_NS = None

# Output directory
OUTPUT_DIR = "figures_md"

# H-bond line color
HBOND_COLOR = "#d62728"

# Matplotlib global configuration
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

# Utilities

def nature_axes(ax):
    """Apply clean publication-style formatting."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4, width=1)


def ensure_output_dir(path=OUTPUT_DIR):
    outdir = Path(path)
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def find_replicate_dirs(base="."):
    """Find directories like analise_replica_1, analise_replica_2, ..."""
    base = Path(base)
    dirs = sorted(
        [d for d in base.glob("analise_replica_*") if d.is_dir()],
        key=lambda x: int(re.search(r"(\d+)$", x.name).group(1))
        if re.search(r"(\d+)$", x.name) else x.name
    )
    return dirs


def build_time_axis(n_points):
    """
    Convert point index to time in ns.
    Priority:
    1) TIME_PER_FRAME_NS
    2) TOTAL_TIME_NS
    """
    if n_points <= 0:
        return np.array([])

    if TIME_PER_FRAME_NS is not None:
        return np.arange(n_points) * TIME_PER_FRAME_NS

    if TOTAL_TIME_NS is not None:
        return np.linspace(0, TOTAL_TIME_NS, n_points)

    raise ValueError(
        "Define TIME_PER_FRAME_NS or TOTAL_TIME_NS in the configuration section."
    )


def load_numeric_dat(filepath):
    """
    Load a generic whitespace-separated numeric .dat file,
    ignoring comment lines (#) and blank lines.
    """
    rows = []
    with open(filepath, "r") as f:
        for line in f:
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

# RMSD

def load_rmsd_file(filepath):
    """
    Expect at least 2 numeric columns:
    col0 = frame/index
    col1 = RMSD
    """
    df = load_numeric_dat(filepath)
    if df.empty or df.shape[1] < 2:
        return None
    df = df.iloc[:, :2].copy()
    df.columns = ["Frame", "RMSD"]
    return df


def plot_rmsd_replicates(rep_dirs, outdir):
    series = []
    labels = []

    for rep in rep_dirs:
        f = rep / "rmsd_ligand.dat"
        if not f.exists():
            continue
        df = load_rmsd_file(f)
        if df is None or df.empty:
            continue
        series.append(df["RMSD"].to_numpy())
        labels.append(rep.name)

    if not series:
        print("[WARN] No valid rmsd_ligand.dat files were found.")
        return

    min_len = min(len(x) for x in series)
    series = np.array([x[:min_len] for x in series])

    mean = series.mean(axis=0)
    std = series.std(axis=0)
    time_ns = build_time_axis(min_len)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for i, y in enumerate(series):
        ax.plot(time_ns, y, lw=1.0, alpha=0.8, label=labels[i])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Ligand RMSD per replicate")
    ax.legend(frameon=False)
    nature_axes(ax)
    fig.savefig(outdir / "rmsd_ligand_replicates.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.plot(time_ns, mean, lw=1.8, label="Mean")
    ax.fill_between(time_ns, mean - std, mean + std, alpha=0.25, label="± SD")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Ligand RMSD (mean across replicates)")
    ax.legend(frameon=False)
    nature_axes(ax)
    fig.savefig(outdir / "rmsd_ligand_mean_std.png")
    plt.close(fig)

# RMSF

def load_rmsf_file(filepath):
    """
    Expect at least 2 numeric columns:
    col0 = residue
    col1 = RMSF
    """
    df = load_numeric_dat(filepath)
    if df.empty or df.shape[1] < 2:
        return None
    df = df.iloc[:, :2].copy()
    df.columns = ["Residue", "RMSF"]
    df["Residue"] = df["Residue"].astype(int)
    return df


def plot_rmsf_replicates(rep_dirs, outdir):
    dfs = []
    labels = []

    for rep in rep_dirs:
        f = rep / "rmsf_ca.dat"
        if not f.exists():
            continue
        df = load_rmsf_file(f)
        if df is None or df.empty:
            continue
        dfs.append(df)
        labels.append(rep.name)

    if not dfs:
        print("[WARN] No valid rmsf_ca.dat files were found.")
        return

    merged = dfs[0].rename(columns={"RMSF": labels[0]})
    for df, label in zip(dfs[1:], labels[1:]):
        merged = merged.merge(
            df.rename(columns={"RMSF": label}),
            on="Residue",
            how="inner"
        )

    value_cols = [c for c in merged.columns if c != "Residue"]
    rmsf_matrix = merged[value_cols].to_numpy()
    mean = rmsf_matrix.mean(axis=1)
    std = rmsf_matrix.std(axis=1)
    residues = merged["Residue"].to_numpy()

    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    for label in value_cols:
        ax.plot(residues, merged[label], lw=1.0, alpha=0.8, label=label)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Protein RMSF (Cα) per replicate")
    ax.legend(frameon=False, ncol=2)
    nature_axes(ax)
    fig.savefig(outdir / "rmsf_ca_replicates.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    ax.plot(residues, mean, lw=1.8, label="Mean")
    ax.fill_between(residues, mean - std, mean + std, alpha=0.25, label="± SD")
    ax.set_xlabel("Residue index")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Protein RMSF (Cα, mean across replicates)")
    ax.legend(frameon=False)
    nature_axes(ax)
    fig.savefig(outdir / "rmsf_ca_mean_std.png")
    plt.close(fig)

# Radius of gyration

def load_rg_file(filepath):
    """
    Expect at least 2 numeric columns:
    col0 = frame/index
    col1 = Rg
    """
    df = load_numeric_dat(filepath)
    if df.empty or df.shape[1] < 2:
        return None
    df = df.iloc[:, :2].copy()
    df.columns = ["Frame", "Rg"]
    return df


def plot_rg_replicates(rep_dirs, outdir):
    series = []
    labels = []

    for rep in rep_dirs:
        f = rep / "rg_protein.dat"
        if not f.exists():
            continue
        df = load_rg_file(f)
        if df is None or df.empty:
            continue
        series.append(df["Rg"].to_numpy())
        labels.append(rep.name)

    if not series:
        print("[WARN] No valid rg_protein.dat files were found.")
        return

    min_len = min(len(x) for x in series)
    series = np.array([x[:min_len] for x in series])

    mean = series.mean(axis=0)
    std = series.std(axis=0)
    time_ns = build_time_axis(min_len)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for i, y in enumerate(series):
        ax.plot(time_ns, y, lw=1.0, alpha=0.8, label=labels[i])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Radius of gyration (Å)")
    ax.set_title("Protein Rg per replicate")
    ax.legend(frameon=False)
    nature_axes(ax)
    fig.savefig(outdir / "rg_protein_replicates.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.plot(time_ns, mean, lw=1.8, label="Mean")
    ax.fill_between(time_ns, mean - std, mean + std, alpha=0.25, label="± SD")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Radius of gyration (Å)")
    ax.set_title("Protein Rg (mean across replicates)")
    ax.legend(frameon=False)
    nature_axes(ax)
    fig.savefig(outdir / "rg_protein_mean_std.png")
    plt.close(fig)

# H-bonds

def load_hbond_timeseries(filepath):
    """
    Parse H-bond time-series output.

    Expected behavior:
    - If file has >= 2 numeric columns, use the 2nd column as the number of H-bonds
    - If file has only 1 numeric column, use it directly
    - X axis is reconstructed as Time (ns) from the number of points
    """
    df = load_numeric_dat(filepath)
    if df.empty:
        return None

    if df.shape[1] >= 2:
        hb = df.iloc[:, 1].to_numpy()
    else:
        hb = df.iloc[:, 0].to_numpy()

    out = pd.DataFrame({
        "Index": np.arange(len(hb)),
        "HBonds": hb
    })
    return out


def plot_hbond_timeseries(rep_dirs, outdir):
    series = []
    labels = []

    for rep in rep_dirs:
        f = rep / "hbond_lig_prot.dat"
        if not f.exists():
            continue
        df = load_hbond_timeseries(f)
        if df is None or df.empty:
            continue
        series.append(df["HBonds"].to_numpy())
        labels.append(rep.name)

    if not series:
        print("[WARN] No valid hbond_lig_prot.dat files were found.")
        return

    min_len = min(len(x) for x in series)
    series = np.array([x[:min_len] for x in series])

    mean = series.mean(axis=0)
    std = series.std(axis=0)
    time_ns = build_time_axis(min_len)

    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    for i, y in enumerate(series):
        ax.plot(time_ns, y, lw=1.0, alpha=0.75, label=labels[i])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title("Ligand–Protein H-bonds per replicate")
    ax.legend(frameon=False)
    nature_axes(ax)
    fig.savefig(outdir / "hbonds_replicates.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.plot(time_ns, mean, lw=1.5, color=HBOND_COLOR)
    ax.fill_between(time_ns, mean - std, mean + std, color=HBOND_COLOR, alpha=0.18)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title("Ligand–Protein H-bonds")
    nature_axes(ax)
    fig.savefig(outdir / "hbonds_mean_std.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.plot(time_ns, mean, lw=1.6, color=HBOND_COLOR)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title("Ligand–Protein H-bonds")
    nature_axes(ax)
    fig.savefig(outdir / "hbonds_mean.png")
    plt.close(fig)

# Summary table

def generate_summary_table(rep_dirs, outdir):
    rows = []

    for rep in rep_dirs:
        row = {"Replicate": rep.name}

        f = rep / "rmsd_ligand.dat"
        if f.exists():
            df = load_rmsd_file(f)
            if df is not None and not df.empty:
                row["RMSD_mean"] = df["RMSD"].mean()
                row["RMSD_sd"] = df["RMSD"].std()

        f = rep / "rmsf_ca.dat"
        if f.exists():
            df = load_rmsf_file(f)
            if df is not None and not df.empty:
                row["RMSF_mean"] = df["RMSF"].mean()
                row["RMSF_sd"] = df["RMSF"].std()

        f = rep / "rg_protein.dat"
        if f.exists():
            df = load_rg_file(f)
            if df is not None and not df.empty:
                row["Rg_mean"] = df["Rg"].mean()
                row["Rg_sd"] = df["Rg"].std()

        f = rep / "hbond_lig_prot.dat"
        if f.exists():
            df = load_hbond_timeseries(f)
            if df is not None and not df.empty:
                row["HBonds_mean"] = df["HBonds"].mean()
                row["HBonds_sd"] = df["HBonds"].std()

        rows.append(row)

    if rows:
        summary = pd.DataFrame(rows)
        summary.to_csv(outdir / "md_replicates_summary.csv", index=False)

# Main

def main():
    base = Path(".")
    outdir = ensure_output_dir(OUTPUT_DIR)
    rep_dirs = find_replicate_dirs(base)

    if not rep_dirs:
        print("No analise_replica_* directories were found.")
        return

    print(f"Replicates found: {[d.name for d in rep_dirs]}")

    plot_rmsd_replicates(rep_dirs, outdir)
    plot_rmsf_replicates(rep_dirs, outdir)
    plot_rg_replicates(rep_dirs, outdir)
    plot_hbond_timeseries(rep_dirs, outdir)
    generate_summary_table(rep_dirs, outdir)

    print(f"Figures saved in: {outdir.resolve()}")


if __name__ == "__main__":
    main()
