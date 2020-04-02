#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import sys

import click

from allfreqs import AlleleFreqs


@click.command()
@click.argument("input_file")
@click.option("--out", "-o", default="all_freqs.csv", show_default=True,
              help="Output filename")
@click.option("--reference", "-r", default=None,
              help="Optional reference file (if not present in INPUT_FILE)")
@click.option("--ambiguous", "-a", default=False, is_flag=True,
              show_default=True,
              help="Show frequencies for ambiguous nucleotides too")
@click.version_option()
def main(input_file, out, reference, ambiguous):
    """Calculate allele frequencies from the given input multialignment.

    Input can be either a fasta or csv file with multialigned sequences,
    which may or may not contain the reference sequence in the first
    position. In the latter case, an additional reference sequence file
    is needed, either in fasta or csv format.
    """
    input_ext = input_file.split(".")[-1]
    if input_ext == "fasta":
        a = AlleleFreqs.from_fasta(input_file, reference, ambiguous)
    elif input_ext == "csv":
        a = AlleleFreqs.from_csv(input_file, reference, ambiguous)
    else:
        click.echo("Input not recognised. "
                   "Please provide either a fasta or csv file.")
        return 1
    a.to_csv(out)
    click.echo(f"Allele frequencies saved to {out}.")

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
