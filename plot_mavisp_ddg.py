#!/usr/bin/env python3

from pathlib import Path
import argparse
import math
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D


# Constants and configuration
EXCLUDE_TOKENS = [
    "classification", "local", "interaction", "interactions",
    "st. dev", "stdev", "std", "sd", "count", "rank"
]
PREFERRED_TOKENS = ["stability", "kcal", "ddg", "ΔΔG", "delta", "dG"]
METHOD_KEYWORDS = {
    "FoldX": ["foldx"],
    "Rosetta": ["rosetta"],
    "RaSP": ["rasp"]
}
DICT_COLORS = {
    "FoldX": "#1f77b4",
    "Rosetta": "#ff7f0e",
    "RaSP": "#2ca02c"
}

# Helper functions
def column_matches_tokens(col: str, tokens: list) -> bool:
    """
    Check if any token in `tokens` is present in the column name `col_name`.
    Case-insensitive.
    
    Args:
        col_name (str): Column name to check.
        tokens (list): List of tokens to search for.
        
    Returns:
        bool: True if any token matches, False otherwise.
    """
    col_lower = col_name.strip().lower()
    return any(tok in col_lower for tok in tokens)

def find_all_method_columns(df: pd.DataFrame) -> dict:
    """
    Identify columns in the DataFrame for each method (FoldX, Rosetta, RaSP),
    separating ΔΔG values from standard deviation columns.

    Args:
        df (pd.DataFrame): Input DataFrame with mutation and method data.

    Returns:
        dict: Mapping of method to {"values": [...], "stdevs": [...]} columns.
    """
    mapping = {}
    for method, keys in METHOD_KEYWORDS.items():
        # Candidate columns containing method keywords, excluding unwanted tokens
        candidates = [c for c in df.columns if any(k in c.lower() for k in keys)]
        candidates = [c for c in candidates if not column_matches_tokens(c, EXCLUDE_TOKENS)]
        
        # Separate ΔΔG value columns from standard deviation columns
        value_cols = [c for c in candidates
                      if not column_matches_tokens(c, ["st. dev", "stdev", "std", "sd"])]
        stdev_cols = [c for c in candidates
                      if column_matches_tokens(c, ["st. dev", "stdev", "std", "sd"])]

        # Prefer columns that contain stability-related tokens
        preferred_values = [c for c in value_cols if column_matches_tokens(c, PREFERRED_TOKENS)]
        if preferred_values:
            value_cols = preferred_values

        mapping[method] = {"values": value_cols, "stdevs": stdev_cols}

    return mapping

# plotting function
def plot_chunk(df_chunk, method_cols, out_path: Path, chunk_idx: int,
               base_name: str, colors_map: dict):
    """
    Plot a chunk of mutations as a grouped bar chart.
    Args:
        df_chunk (pd.DataFrame): Chunk of mutation data.
        method_cols (dict): Dictionary of method columns with 'values' and 'stdevs'.
        out_path (Path): Directory to save output plots.
        chunk_idx (int): Index of this chunk.
        base_name (str): Base name for output files.
        colors_map (dict): Mapping of method names to colors.
    """    
    mutations = df_chunk.index.astype(str).tolist()
    n_mutations = len(mutations)
    methods = list(method_cols.keys())
    n_methods = len(methods)
    max_bars = max(len(method_cols[m]["values"]) for m in methods)
    bar_width = 0.8 / max(max_bars * n_methods, 1)
    x_positions = np.arange(n_mutations)

    fig, ax = plt.subplots(figsize=(max(6, n_mutations * 0.7), 6))

    for i, method in enumerate(methods):
        value_columns = method_cols[method]["values"]
        stdev_columns = method_cols[method]["stdevs"]

        for j, val_col in enumerate(value_columns):
            stdev_col = stdev_columns[j] if j < len(stdev_columns) else None
            if val_col not in df_chunk.columns:
                continue

            y_values = pd.to_numeric(df_chunk[val_col], errors="coerce").to_numpy()
            y_err = pd.to_numeric(df_chunk[stdev_col], errors="coerce").to_numpy() \
                if stdev_col and stdev_col in df_chunk.columns else None

            ax.bar(x_positions + (i * max_bars + j) * bar_width,
                   y_values,
                   width=bar_width,
                   yerr=y_err,
                   capsize=3 if y_err is not None else 0,
                   color=colors_map.get(method),
                   label=method if j == 0 else "")

    ax.axhline(0, linestyle="--", linewidth=1)
    total_bars = sum(len(method_cols[m]["values"]) for m in methods)
    ax.set_xticks(x_positions + bar_width * total_bars / 2)
    ax.set_xticklabels(mutations, rotation=45, ha="right")
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_xlabel("Mutation")
    ax.set_title(f"{base_name} — mutations {chunk_idx:02d}")

    ax.legend(frameon=True)
    plt.tight_layout()

    pdf_file = out_path / f"{base_name}_{chunk_idx:02d}.pdf"
    png_file = out_path / f"{base_name}_{chunk_idx:02d}.png"
    fig.savefig(pdf_file)
    fig.savefig(png_file, dpi=300)
    plt.close(fig)

    print(f"Saved: {pdf_file} and {png_file}")

# Main
def main():
     """Parse arguments, read CSV, split into chunks, and plot ΔΔG values."""
    parser = argparse.ArgumentParser(description="Plot stability ΔΔG columns in chunks.")
    parser.add_argument("-c", "--csv", required=True, help="Input CSV file.")
    parser.add_argument("-o", "--out", required=False,
                        default=None, help="Output directory or filename prefix.")
    parser.add_argument("-n", "--chunk-size", type=int, default=10,
                        help="Number of mutations per plot (default: 10).")
    
    args = parser.parse_args()
    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"ERROR: CSV file not found: {csv_path}", file=sys.stderr)
        sys.exit(2)

    # Set output directory
    out_dir = Path(args.out) if args.out else Path.cwd() / csv_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read CSV and set mutation column as index
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]
    df.set_index("Mutation", inplace=True)


    # Identify method columns
    method_columns = find_all_method_columns(df)
    colors_map = {m: METHOD_COLORS.get(m, "#7f7f7f") for m in method_columns.keys()}

    # Split DataFrame into chunks and plot each
    chunk_size = args.chunk_size
    total_mutations = len(df)
    n_chunks = math.ceil(total_mutations / chunk_size)
    print(f"Total mutations: {total_mutations} → {n_chunks} chunks")

    for idx in range(n_chunks):
        start, stop = idx * chunk_size, min((idx + 1) * chunk_size, total_mutations)
        df_chunk = df.iloc[start:stop]
        plot_chunk(df_chunk, method_columns, out_dir, idx + 1, csv_path.stem, colors_map)

    print("Done.")

if __name__ == "__main__":
    main()