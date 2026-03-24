#!/usr/bin/env python3
"""
Export a cell-by-gene count matrix from AnnData for pySCENIC text input.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export a count matrix for pySCENIC.")
    parser.add_argument("--input", required=True, help="Input .h5ad file.")
    parser.add_argument("--output", required=True, help="Output .csv or .tsv file.")
    parser.add_argument(
        "--layer",
        default=None,
        help="Optional AnnData layer to export. Default uses .X.",
    )
    parser.add_argument(
        "--transpose",
        action="store_true",
        help="Export as gene-by-cell instead of the pySCENIC-recommended cell-by-gene format.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(input_path)
    matrix = adata.layers[args.layer] if args.layer else adata.X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()

    df = pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
    if args.transpose:
        df = df.T

    sep = "\t" if output_path.suffix.lower() in {".tsv", ".txt"} else ","
    df.to_csv(output_path, sep=sep)

    print(f"Matrix written to: {output_path}")
    print(f"Shape: {df.shape[0]} x {df.shape[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
