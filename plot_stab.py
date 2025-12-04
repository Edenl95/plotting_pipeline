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
    "classification", "classif", "local", "interaction", "interactions",
    "st. dev", "stdev", "std", "sd", "count", "rank"
]
PREFERRED_TOKENS = ["stability", "kcal", "ddg", "ΔΔG", "delta", "ΔΔG", "dG"]

METHOD_KEYWORDS = {
    "FoldX": ["foldx", "fold x"],
    "Rosetta": ["rosetta"],
    "RaSP": ["rasp", "rosp"]  # allow common variants
}

def normalise_col(col: str) -> str:
    return col.strip().lower()

def column_is_excluded(col: str) -> bool:
    low = normalise_col(col)
    return any(tok in low for tok in EXCLUDE_TOKENS)

def column_is_preferred(col: str) -> bool:
    low = normalise_col(col)
    return any(tok in low for tok in PREFERRED_TOKENS)

def find_method_columns(df: pd.DataFrame):
    """
    For each method in METHOD_KEYWORDS, return a dict:
      method -> {"value_cols": [col, ...], "stdev_cols": [col, ...]}
    Preference is given to columns that contain PREFERRED_TOKENS.
    We exclude columns that contain EXCLUDE_TOKENS.
    """
    cols = list(df.columns)
    mapping = {}
    for method, keys in METHOD_KEYWORDS.items():
        candidates = []
        for col in cols:
            low = normalise_col(col)
            if any(k in low for k in keys):
                candidates.append(col)

        # remove excluded candidates
        candidates = [c for c in candidates if not column_is_excluded(c)]

        # Separate value columns vs stdev columns heuristically
        value_candidates = []
        stdev_candidates = []
        for c in candidates:
            low = normalise_col(c)
            if any(tok in low for tok in ["st. dev", "stdev", "std", "sd"]):
                stdev_candidates.append(c)
            else:
                value_candidates.append(c)

        # If we have many value candidates, prefer those with "stability" or "kcal" tokens
        if value_candidates:
            preferred = [c for c in value_candidates if column_is_preferred(c)]
            if preferred:
                value_candidates = preferred

        # Fallback: if we detected only stdev-like columns (rare), try to recover value columns from original set
        if not value_candidates and candidates:
            fallback = [c for c in candidates if not any(tok in normalise_col(c) for tok in ["classification"])]
            value_candidates = fallback

        mapping[method] = {
            "value_cols": value_candidates,
            "stdev_cols": stdev_candidates
        }
    return mapping

# Plotting
def plot_chunk(df_chunk: pd.DataFrame, mutation_col: str, method_cols: dict,
               outpath: Path, chunk_idx: int, base_name: str, colors_map: dict,
               legend_loc="upper right"):
  
    mutations = df_chunk[mutation_col].astype(str).tolist()
    n = len(mutations)
    methods = list(method_cols.keys())
    n_methods = len(methods)
    bar_width = 0.8 / max(n_methods, 1)  # keep group width ~0.8

    x = np.arange(n)
    fig, ax = plt.subplots(figsize=(max(6, n * 0.5), 6))

    for i, method in enumerate(methods):
        col_info = method_cols[method]
        valcol = col_info.get("value")
        sdcol = col_info.get("stdev")
        if valcol is None:
            continue
        y = df_chunk[valcol].to_numpy(dtype=float)
        if sdcol is not None and sdcol in df_chunk.columns:
            yerr = df_chunk[sdcol].to_numpy(dtype=float)
            ax.bar(x + i * bar_width, y, width=bar_width, label=method,
                   yerr=yerr, capsize=3, color=colors_map.get(method))
        else:
            ax.bar(x + i * bar_width, y, width=bar_width, label=method,
                   color=colors_map.get(method))

    ax.axhline(0, linestyle="--", linewidth=1)  # dashed zero line
    ax.set_xticks(x + bar_width * (n_methods - 1) / 2)
    ax.set_xticklabels(mutations, rotation=45, ha="right")
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_xlabel("Mutation")
    ax.set_title(f"{base_name} — mutations {chunk_idx:02d}")

    # Add st.dev to legend
    handles, labels = ax.get_legend_handles_labels()
    err_handle = Line2D([0], [0], color='black', lw=1, linestyle='', marker='|', markersize=10)
    handles.append(err_handle)
    labels.append("St. dev")
    ax.legend(handles, labels, loc=legend_loc, frameon=True)

    plt.tight_layout()
    out_file_pdf = outpath / f"{base_name}_{chunk_idx:02d}.pdf"
    out_file_png = outpath / f"{base_name}_{chunk_idx:02d}.png"

    fig.savefig(out_file_pdf)
    fig.savefig(out_file_png, dpi=300)  # high-resolution PNG for presentations
    plt.close(fig)
    print(f"Saved: {out_file_pdf} and {out_file_png}")


# Main
def main():
    parser = argparse.ArgumentParser(description="Plot stability ΔΔG columns in chunks.")
    parser.add_argument("-c", "--csv", required=True, help="Input CSV file (simple/ensemble mode CSV).")
    parser.add_argument("-o", "--out", required=False, default=None,
                        help="Output directory or filename prefix. If omitted, uses CSV basename in current dir.")
    parser.add_argument("-n", "--chunk-size", type=int, default=10,
                        help="Number of mutations per plot (default: 10).")
    parser.add_argument("--mutation-col", default="Mutation", help="Name of mutation column (default 'Mutation').")
    parser.add_argument("--legend-loc", default="upper right", help="Legend location (default 'upper right').")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"ERROR: CSV file not found: {csv_path}", file=sys.stderr)
        sys.exit(2)

    out_path = Path(args.out) if args.out else Path.cwd() / csv_path.stem
    if out_path.suffix == "":
        out_dir = out_path
    else:
        # if user provided something with a suffix, treat as dir prefix
        out_dir = out_path.parent / out_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]

    # Validate mutation column
    mutation_col = args.mutation_col
    if mutation_col not in df.columns:
        print(f"ERROR: Mutation column '{mutation_col}' not found in CSV columns.", file=sys.stderr)
        print("Available columns:", ", ".join(df.columns), file=sys.stderr)
        sys.exit(2)

    # Find method columns
    detected = find_method_columns(df)

    # For each method choose a single value column (and possibly stdev) following heuristics:
    method_cols = {}
    for method, info in detected.items():
        value_cols = info["value_cols"]
        stdev_cols = info["stdev_cols"]

        if not value_cols:
            # try a relaxed approach: look for any column with method name but exclude classification/local
            fallback = [c for c in df.columns if any(k in normalise_col(c) for k in METHOD_KEYWORDS[method])]
            fallback = [c for c in fallback if not column_is_excluded(c)]
            value_cols = fallback

        if not value_cols:
            # no column found for this method
            continue

        # Prefer columns that contain "stability" or "kcal" or "ddg"
        preferred = [c for c in value_cols if column_is_preferred(c)]
        chosen_value = preferred[0] if preferred else value_cols[0]

        chosen_stdev = None
        if stdev_cols:
            chosen_stdev = stdev_cols[0]
        else:
            # try to find an stdev-like column near the chosen value col name
            pattern = re.escape(chosen_value)
            # look for columns that contain "st. dev" or "std" together with method name
            for c in df.columns:
                low = normalise_col(c)
                if any(tok in low for tok in ["st. dev", "stdev", "std", "sd"]) and any(k in low for k in METHOD_KEYWORDS[method]):
                    chosen_stdev = c
                    break

        method_cols[method] = {"value": chosen_value, "stdev": chosen_stdev}

    if not method_cols:
        print("ERROR: No FoldX/Rosetta/RaSP stability columns identified.", file=sys.stderr)
        print("Columns present:", ", ".join(df.columns), file=sys.stderr)
        sys.exit(2)

    # Coerce chosen columns to numeric and drop rows where all chosen methods are NaN
    cols_to_check = [info["value"] for info in method_cols.values() if info["value"]]
    for c in cols_to_check:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Optionally coerce stdev columns
    for info in method_cols.values():
        s = info.get("stdev")
        if s and s in df.columns:
            df[s] = pd.to_numeric(df[s], errors="coerce")

    # Drop rows where all value columns are NaN (nothing to plot)
    value_cols_present = [info["value"] for info in method_cols.values() if info["value"] in df.columns]
    before = len(df)
    df = df.dropna(subset=value_cols_present, how="all").reset_index(drop=True)
    after = len(df)
    if after < before:
        print(f"Info: Dropped {before-after} rows with no numeric values in detected stability columns.")

    # Determine consistent colors for methods (stable mapping)
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", None)
    colors_map = {}
    if color_cycle:
        for i, method in enumerate(sorted(method_cols.keys())):
            colors_map[method] = color_cycle[i % len(color_cycle)]
    else:
        # fallback: let matplotlib choose
        for method in method_cols.keys():
            colors_map[method] = None

    # Create chunks and plot
    chunk_size = args.chunk_size
    total = len(df)
    n_chunks = math.ceil(total / chunk_size)
    base_name = csv_path.stem
    print(f"Found methods and selected columns:")
    for method, info in method_cols.items():
        print(f"  - {method}: value='{info['value']}' stdev='{info.get('stdev')}'")
    print(f"Total mutations to plot: {total}. Chunk size: {chunk_size}. Will create {n_chunks} files in '{out_dir}'")

    for idx in range(n_chunks):
        start = idx * chunk_size
        stop = min((idx + 1) * chunk_size, total)
        df_chunk = df.iloc[start:stop].copy()
        # Use only mutation + chosen columns
        keep_cols = [mutation_col] + [info["value"] for info in method_cols.values() if info["value"] in df_chunk.columns]
        # Include stdev cols if present
        for info in method_cols.values():
            s = info.get("stdev")
            if s and s in df_chunk.columns and s not in keep_cols:
                keep_cols.append(s)
        df_chunk = df_chunk[keep_cols]
        # Rename chosen value cols to method names for plotting convenience
        rename_map = {info["value"]: method for method, info in method_cols.items() if info["value"] in df_chunk.columns}
        df_chunk = df_chunk.rename(columns=rename_map)

        # Update method_cols for this chunk with actual column names in the chunk
        method_cols_for_plot = {}
        for method, info in method_cols.items():
            valcol_name = rename_map.get(info["value"])  # method name if present
            if valcol_name is None and info["value"] in df_chunk.columns:
                valcol_name = info["value"]
            if valcol_name is None:
                # method not present for this chunk (all NaN or missing)
                method_cols_for_plot[method] = {"value": None, "stdev": None}
                continue
            # stdev mapping: use original stdev column name
            stdev_col = info.get("stdev")
            if stdev_col and stdev_col in df_chunk.columns:
                method_cols_for_plot[method] = {"value": valcol_name, "stdev": stdev_col}
            else:
                method_cols_for_plot[method] = {"value": valcol_name, "stdev": None}

        plot_chunk(df_chunk, mutation_col, method_cols_for_plot, out_dir, idx+1, base_name, colors_map, legend_loc=args.legend_loc)

    print("Done.")

if __name__ == "__main__":
    main()