# Stability Plot (plot_mavisp_ddg.py)
## Description

The plot_mavisp_ddg.py script automates the visualization of protein mutation stability data from MAVISp-generated CSV files. It generates bar plots of ΔΔG (kcal/mol) values for multiple computational methods: FoldX, Rosetta, and RaSP. Standard deviations are automatically plotted whenever the corresponding columns (st. dev, stdev, std, or sd) are present in the CSV.

The script is optimized for large datasets by splitting mutations into chunks, generating PDF and high-resolution PNG plots for each chunk.


The script automatically:
* Detects stability columns for FoldX, Rosetta, and RaSP.
* Excludes irrelevant columns (e.g., classification, count, rank).
* Converts chosen columns to numeric values and drops rows without valid data.
* Splits mutations into chunks for plotting if the dataset is large.
* Generates PDF and high-resolution PNG plots for each chunk.

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

## Output

For each chunk of mutations, the script generates two files in the output directory:

| File | Description |
|------|-------------|
| `CSVNAME_01.pdf` | PDF plot of ΔΔG values for the first chunk of mutations |
| `CSVNAME_01.png` | PNG plot of ΔΔG values for the first chunk of mutations |

Each plot includes:
* ΔΔG values per mutation for each method.
* Error bars if standard deviation columns are present.
* Legend showing method names and a marker for standard deviation.

## Options
```bash
-h, --help            show this help message and exit
-c CSV, --csv CSV     Input CSV file (simple or ensemble mode)
-o OUT, --out OUT     Output directory (default: CSV basename in current directory)
-n CHUNK_SIZE, --chunk-size CHUNK_SIZE
                      Number of mutations per plot (default: 10)
--mutation-col MUTATION_COL
                      Name of the mutation column (default 'Mutation')
--legend-loc LEGEND_LOC
                      Legend location in the plot (default 'upper right')
```

## Requirements
Python 3.10+
Packages: numpy, pandas, matplotlib

Example module load (if using a module system):
```bash
module load python/3.10/modulefile
```

## Usage
1. Copy the CSV file into the folder where you will run the script, or use the full path to the CSV with -c.

2. Run the script:
```bash
# CSV in same folder, default chunk size 10
python3 plot_mavisp_ddg.py -c my_mutations.csv 

# CSV in a different folder
python3 plot_mavisp_ddg.py -c /full/path/to/my_mutations.csv 

# Specify output folder and chunk size
python3 plot_mavisp_ddg.py -c my_mutations.csv -o plots/ -n 15

```

**This will:**
* Read the CSV file (my_mutations.csv)
* Detect stability columns for FoldX, Rosetta, RaSP
* Split mutations into chunks of 10
* Save PDF and PNG plots into the plots/ directory
