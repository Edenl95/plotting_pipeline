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

# Utility and detection functions
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

# Unified column matching
def column_matches_tokens(col: str, tokens: list) -> bool:
    low = col.strip().lower()
    return any(tok in low for tok in tokens)


def find_method_columns(df: pd.DataFrame):
    mapping = {}
    for method, keys in METHOD_KEYWORDS.items():
        candidates = [c for c in df.columns if any(k in c.strip().lower() for k in keys)]
        candidates = [c for c in candidates if not column_matches_tokens(c, EXCLUDE_TOKENS)]

        value_candidates = [c for c in candidates if not column_matches_tokens(c, ["st. dev", "stdev", "std", "sd"])]
        stdev_candidates = [c for c in candidates if column_matches_tokens(c, ["st. dev", "stdev", "std", "sd"])]

        # Prefer columns with preferred tokens
        preferred_values = [c for c in value_candidates if column_matches_tokens(c, PREFERRED_TOKENS)]
        if preferred_values:
            value_candidates = preferred_values

        # Fallback if no value columns
        if not value_candidates and candidates:
            value_candidates = [c for c in candidates if not column_matches_tokens(c, ["classification"])]

        # **Fix**: store None only if the column really doesn't exist
        val_col = value_candidates[0] if value_candidates else None
        sd_col = stdev_candidates[0] if stdev_candidates else None

        # Extra fallback: sometimes stdev column name includes the method name
        if not sd_col:
            for c in df.columns:
                if any(k in c.strip().lower() for k in keys) and column_matches_tokens(c, ["st. dev", "stdev", "std", "sd"]):
                    sd_col = c
                    break

        mapping[method] = {"value": val_col, "stdev": sd_col}
    return mapping

# plotting function
def plot_chunk(df_chunk, mutation_col, method_cols, outpath: Path, chunk_idx: int, base_name: str, colors_map: dict, legend_loc="upper right"):
    mutations = df_chunk[mutation_col].astype(str).tolist()
    n = len(mutations)
    methods = list(method_cols.keys())
    n_methods = len(methods)
    bar_width = 0.8 / max(n_methods, 1)
    x = np.arange(n)

    fig, ax = plt.subplots(figsize=(max(6, n * 0.5), 6))

    for i, method in enumerate(methods):
        valcol = method_cols[method]["value"]
        sdcol = method_cols[method]["stdev"]
        if valcol is None or valcol not in df_chunk.columns:
            continue
        y = pd.to_numeric(df_chunk[valcol], errors="coerce").to_numpy()
        yerr = pd.to_numeric(df_chunk[sdcol], errors="coerce").to_numpy() if sdcol in df_chunk.columns else None
        ax.bar(x + i * bar_width, y, width=bar_width, label=method,
               yerr=yerr, capsize=3 if yerr is not None else 0,
               color=colors_map.get(method))

    ax.axhline(0, linestyle="--", linewidth=1)
    ax.set_xticks(x + bar_width * (n_methods - 1) / 2)
    ax.set_xticklabels(mutations, rotation=45, ha="right")
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_xlabel("Mutation")
    ax.set_title(f"{base_name} — mutations {chunk_idx:02d}")

    handles, labels = ax.get_legend_handles_labels()
    err_handle = Line2D([0], [0], color='black', lw=1, linestyle='', marker='|', markersize=10)
    handles.append(err_handle)
    labels.append("St. dev")
    ax.legend(handles, labels, loc=legend_loc, frameon=True)

    plt.tight_layout()
    out_file_pdf = outpath / f"{base_name}_{chunk_idx:02d}.pdf"
    out_file_png = outpath / f"{base_name}_{chunk_idx:02d}.png"
    fig.savefig(out_file_pdf)
    fig.savefig(out_file_png, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_file_pdf} and {out_file_png}")

# Main
def main():
    parser = argparse.ArgumentParser(description="Plot stability ΔΔG columns in chunks.")
    parser.add_argument("-c", "--csv", required=True, help="Input CSV file.")
    parser.add_argument("-o", "--out", required=False, default=None, help="Output directory or filename prefix.")
    parser.add_argument("-n", "--chunk-size", type=int, default=10, help="Number of mutations per plot (default: 10).")
    parser.add_argument("--mutation-col", default="Mutation", help="Name of mutation column.")
    parser.add_argument("--legend-loc", default="upper right", help="Legend location.")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"ERROR: CSV file not found: {csv_path}", file=sys.stderr)
        sys.exit(2)

    out_dir = Path(args.out) if args.out else Path.cwd() / csv_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]

    if args.mutation_col not in df.columns:
        print(f"ERROR: Mutation column '{args.mutation_col}' not found.", file=sys.stderr)
        sys.exit(2)

    mutation_col = args.mutation_col
    chunk_size = args.chunk_size
    base_name = csv_path.stem
    legend_loc = args.legend_loc

    method_cols = find_method_columns(df)
    colors_map = {m: DICT_COLORS.get(m, "#7f7f7f") for m in method_cols.keys()}

    total = len(df)
    n_chunks = math.ceil(total / chunk_size)
    print(f"Total mutations: {total} → {n_chunks} chunks")

    for idx in range(n_chunks):
        start, stop = idx * chunk_size, min((idx + 1) * chunk_size, total)
        df_chunk = df.iloc[start:stop]
        plot_chunk(df_chunk, mutation_col, method_cols, out_dir, idx + 1, base_name, colors_map, legend_loc)

    print("Done.")

    # run
if __name__ == "__main__":
    main()