#!/usr/bin/env python3
import argparse
import sys
import re

import yaml


def revcomp(seq: str) -> str:
    """Return reverse complement of an A/C/G/T sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def main():
    parser = argparse.ArgumentParser(
        description="Generate barcodes.txt from new debarcoding YAML config (reverse-complementing barcodes)."
    )
    parser.add_argument("--config", help="Path to config_new_debarcode.yaml")
    parser.add_argument(
        "-o",
        "--output",
        default="barcodes.txt",
        help="Output barcodes file (default: barcodes.txt)",
    )
    parser.add_argument("--reverse-complement", action="store_true", help="Output reverse-complemented barcodes")
    args = parser.parse_args()

    # Load YAML
    try:
        with open(args.config) as fh:
            cfg = yaml.safe_load(fh)
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: Config file not found: {args.config}\n")
        sys.exit(1)

    # Get sample name (optional, just for logging)
    sample_name = cfg.get("sample_name", "UNKNOWN_SAMPLE")

    # Get barcodes mapping (now at top level)
    barcodes = cfg.get("barcodes")
    if not isinstance(barcodes, dict):
        sys.stderr.write("ERROR: config must contain a top-level 'barcodes' mapping.\n")
        sys.exit(1)

    # Prepare output
    with open(args.output, "w") as out_fh:
        for name, seq in barcodes.items():
            seq_str = str(seq).strip()
            # Only keep barcodes that are pure A/C/G/T
            if not re.fullmatch(r"[ACGTacgt]+", seq_str):
                sys.stderr.write(
                    f"Skipping barcode '{name}' with non-ACGT characters: '{seq_str}'\n"
                )
                continue
            rc = revcomp(seq_str)
            out_fh.write(f"{seq_str}\t{rc}\n")

    sys.stderr.write(
        f"Wrote barcodes for sample '{sample_name}' to {args.output}\n"
    )


if __name__ == "__main__":
    main()