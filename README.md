# Stability Plot
## Description

The plot_stab script automates the visualization of protein mutation stability data from MAVISp-generated CSV files. Its primary goal is to generate bar plots of ΔΔG values for multiple methods (FoldX, Rosetta, RaSP), optionally including standard deviations. The script processes CSV files, detects the relevant columns for each method, and splits the data into chunks for plotting if there are many mutations.
**Important**: Make sure the input CSV file is either copied into the working folder where you run the script, or provide the full path to the CSV file with the -c option.

The script automatically:
* Detects FoldX, Rosetta, and RaSP stability columns in the CSV.
* Excludes irrelevant columns (e.g., classification or count columns).
* Converts chosen columns to numeric and drops rows without valid values.
* Generates PDF and high-resolution PNG plots of ΔΔG per mutation.

## Output

For each chunk of mutations, the script generates two files in the output directory:

| File | Description |
|------|-------------|
| `plot_01.pdf` | PDF plot of ΔΔG values |
| `plot_01.png` | PNG plot of ΔΔG values |

Each plot includes:
* ΔΔG values per mutation for each method
* Error bars if standard deviation columns are present
* Legend with method names and st. dev indicator

## Options
```bash
-h, --help            show this help message and exit
-c CSV, --csv CSV     Input CSV file (simple/ensemble mode CSV)
-o OUT, --out OUT     Output directory or filename prefix (default: CSV basename in current directory)
-n CHUNK_SIZE, --chunk-size CHUNK_SIZE
                      Number of mutations per plot (default: 10)
--mutation-col MUTATION_COL
                      Name of mutation column (default 'Mutation')
--legend-loc LEGEND_LOC
                      Legend location in the plot (default 'upper right')
```

## Requirements
Python 3.10+
Packages: numpy, pandas, matplotlib

Example module load (if using a module system):
module load python/3.10/modulefile

## Usage
1. Copy the CSV file into the folder where you will run the script, or use the full path to the CSV with -c.
2. Run the script:
```bash
# CSV in same folder
python3 plot_stab.py -c my_mutations.csv -o plots/ -n 10

# CSV in a different folder
python3 plot_stab.py -c /full/path/to/my_mutations.csv -o plots/ -n 10
```

**This will:**
* Read the CSV file (my_mutations.csv)
* Detect stability columns for FoldX, Rosetta, RaSP
* Split mutations into chunks of 10
* Save PDF and PNG plots into the plots/ directory
