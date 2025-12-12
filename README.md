# Stability Plot (plot_mavisp_ddg.py)
## Description

plot_mavisp_ddg.py automates the visualization of protein mutation stability data from MAVISp-generated CSV files. It generates grouped bar plots of ΔΔG (kcal/mol) values for multiple computational methods: FoldX, Rosetta, and RaSP. Standard deviations are plotted automatically when the corresponding columns (st. dev, stdev, std, or sd) are present.

The script handles large datasets efficiently by splitting mutations into chunks, generating both PDF and high-resolution PNG plots for each chunk.

The script automatically:
* Detects ΔΔG columns for FoldX, Rosetta, and RaSP.
* Excludes irrelevant columns (e.g., classification, count, rank).
* Converts values to numeric and ignores rows without valid data.
* Splits large datasets into chunks for plotting.
* Produces PDF and PNG plots for each chunk.

## Requirements
Python 3.10+
Packages: numpy, pandas, matplotlib

Example module load (if using a module system):
```bash
module load python/3.10/modulefile
```

**Important**: Make sure the input CSV file is either in the working folder or provide the full path using the -c option.

By default, if no output folder is specified with -o, the script creates a new folder in the current directory using the CSV filename as the folder name (without extension). For example:
```bash
python3 plot_mavisp_ddg.py -c ABI1-simple_mode.csv
# output saved in ./ABI1-simple_mode/
```

You can override this by specifying an output folder:
```bash
python3 plot_mavisp_ddg.py -c ABI1-simple_mode.csv -o plots/
# output saved in ./plots/
```

## Usage
1. Copy the CSV file into the folder where you will run the script, or use the full path to the CSV with -c.

2. Run the script:
```bash
# CSV in the same folder, default chunk size 10
python3 plot_mavisp_ddg.py -c my_mutations.csv 

# CSV in a different folder
python3 plot_mavisp_ddg.py -c /full/path/to/my_mutations.csv 

# Specify output folder and chunk size
python3 plot_mavisp_ddg.py -c my_mutations.csv -o plots/ -n 15
```

If no output folder is specified with -o, a folder will be created in the current directory using the CSV filename (without extension).

Example:
```bash
python3 plot_mavisp_ddg.py -c ABI1-simple_mode.csv
# Output saved in ./ABI1-simple_mode/
```
## Options
| Flag                                     | Description                                                   |
| ---------------------------------------- | ------------------------------------------------------------- |
| `-h, --help`                             | Show help message and exit                                    |
| `-c CSV, --csv CSV`                      | Input CSV file (simple or ensemble mode)                      |
| `-o OUT, --out OUT`                      | Output directory (default: CSV basename in current directory) |
| `-n CHUNK_SIZE, --chunk-size CHUNK_SIZE` | Number of mutations per plot (default: 10)                    |

## Output

For each chunk of mutations, the script generates two files in the output directory:

| File | Description |
|------|-------------|
| `CSVNAME_01.pdf` | PDF plot of ΔΔG values for the first chunk |
| `CSVNAME_01.png` | PNG plot of ΔΔG values for the first chunk |

Each plot includes:
* ΔΔG values per mutation for each method.
* Error bars for standard deviation (if present).
* Legend showing method names.

## Notes
* Make sure your CSV has a column named Mutation, which will be used as the index.
* Columns for ΔΔG values should include method keywords (foldx, rosetta, rasp) and * optionally preferred tokens (stability, kcal, ddg, ΔΔG, delta, dG).
* Standard deviation columns (st. dev, stdev, std, sd) are automatically detected if present.
* Large datasets are split into chunks for readability and clarity in the plots.

**This will:**
* Read the CSV file (my_mutations.csv)
* Detect stability columns for FoldX, Rosetta, RaSP
* Split mutations into chunks of 10
* Save PDF and PNG plots into the plots/ directory
