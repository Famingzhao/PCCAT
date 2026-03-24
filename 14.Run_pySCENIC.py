#!/usr/bin/env python3
"""
Run a reproducible pySCENIC workflow from a single command.

This script wraps the official pySCENIC CLI:
1. GRN inference
2. cisTarget motif pruning
3. AUCell scoring

Input expression matrix must be cell-by-gene for text input files.
Supported input formats:
- .csv
- .tsv / .txt
- .loom
- .h5ad
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path


def run_command(cmd: list[str]) -> None:
    print("Running:", " ".join(shlex.quote(x) for x in cmd), flush=True)
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pySCENIC from a single wrapper script.")
    parser.add_argument(
        "--input",
        required=True,
        help="Expression matrix file (.csv/.tsv/.txt/.loom/.h5ad). For text files, rows=cells and columns=genes.",
    )
    parser.add_argument(
        "--tfs",
        required=True,
        help="Transcription factor list file used by pySCENIC grn.",
    )
    parser.add_argument(
        "--motif-dbs",
        nargs="+",
        required=True,
        help="One or more cisTarget ranking databases in Feather v2 format.",
    )
    parser.add_argument(
        "--motif-annotations",
        required=True,
        help="Motif annotation table, for example motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for all pySCENIC outputs.",
    )
    parser.add_argument(
        "--method",
        default="grnboost2",
        choices=["grnboost2", "genie3"],
        help="Gene regulatory network inference method.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=10,
        help="Number of CPU workers for each pySCENIC step.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=777,
        help="Random seed used by pySCENIC where supported.",
    )
    parser.add_argument(
        "--sparse",
        action="store_true",
        help="Use sparse mode in `pyscenic grn`.",
    )
    parser.add_argument(
        "--mask-dropouts",
        action="store_true",
        help="Use dropout masking in `pyscenic ctx`.",
    )
    parser.add_argument(
        "--cell-id-attribute",
        default="CellID",
        help="Cell attribute name for loom input/output.",
    )
    parser.add_argument(
        "--gene-attribute",
        default="Gene",
        help="Gene attribute name for loom input/output.",
    )
    return parser.parse_args()


def build_paths(output_dir: Path) -> dict[str, Path]:
    return {
        "adj": output_dir / "pyscenic_adj.tsv",
        "reg": output_dir / "pyscenic_regulons.csv",
        "auc": output_dir / "pyscenic_auc.csv",
    }


def main() -> int:
    args = parse_args()

    input_path = Path(args.input).expanduser().resolve()
    tfs_path = Path(args.tfs).expanduser().resolve()
    motif_dbs = [str(Path(x).expanduser().resolve()) for x in args.motif_dbs]
    motif_annotations = Path(args.motif_annotations).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = build_paths(output_dir)

    grn_cmd = [
        "pyscenic",
        "grn",
        str(input_path),
        str(tfs_path),
        "-o",
        str(paths["adj"]),
        "--method",
        args.method,
        "--num_workers",
        str(args.num_workers),
        "--seed",
        str(args.seed),
    ]
    if args.sparse:
        grn_cmd.append("--sparse")

    ctx_cmd = [
        "pyscenic",
        "ctx",
        str(paths["adj"]),
        *motif_dbs,
        "--annotations_fname",
        str(motif_annotations),
        "--expression_mtx_fname",
        str(input_path),
        "--output",
        str(paths["reg"]),
        "--num_workers",
        str(args.num_workers),
    ]
    if args.mask_dropouts:
        ctx_cmd.append("--mask_dropouts")
    if input_path.suffix.lower() == ".loom":
        ctx_cmd.extend(
            [
                "--cell_id_attribute",
                args.cell_id_attribute,
                "--gene_attribute",
                args.gene_attribute,
            ]
        )

    auc_cmd = [
        "pyscenic",
        "aucell",
        str(input_path),
        str(paths["reg"]),
        "-o",
        str(paths["auc"]),
        "--num_workers",
        str(args.num_workers),
    ]
    if input_path.suffix.lower() == ".loom":
        auc_cmd.extend(
            [
                "--cell_id_attribute",
                args.cell_id_attribute,
                "--gene_attribute",
                args.gene_attribute,
            ]
        )

    run_command(grn_cmd)
    run_command(ctx_cmd)
    run_command(auc_cmd)

    print("\npySCENIC workflow completed successfully.")
    print(f"Adjacencies: {paths['adj']}")
    print(f"Regulons:    {paths['reg']}")
    print(f"AUC matrix:  {paths['auc']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
