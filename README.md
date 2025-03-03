# Assays

## Instructions

1. Install `uv`, a python package and environment manager:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2. Get data files. You will need at least two files: `clean_data.csv` and `target_quality.tsv`. Optionally, you can start with the processed version of this file, `scores.npy`, and skip to step #4.

3. Run the following command to process the `clean_data.csv` file and create the `scores.npy` file:

```bash
uv run main.py process
```

4. Run the following command to produce the analysis files from `scores.npy` and `target_quality.tsv`:

```bash
uv run main.py analyze
```

This should produce four Excel files:

- `highest_scores.xlsx`
- `allgood_lowest_scores.xlsx`
- `allgood_highest_scores.xlsx`
- `allgood_zero_scores.xlsx`