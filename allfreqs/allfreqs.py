#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from typing import Optional

from cached_property import cached_property
import pandas as pd
from skbio import TabularMSA, DNA, Sequence

from allfreqs.classes import Reference, MultiAlignment


class AlleleFreqs:
    """Class used to calculate allele frequencies from a multialignment.

    Input can be either a fasta or csv file with multialigned sequences,
    which may or may not contain the reference sequence in the first
    position. In the latter case, an additional reference sequence file
    is needed, either in fasta or csv format.
    """

    def __init__(self, multialg: MultiAlignment, reference: Reference):
        self.multialg = multialg
        self.reference = reference
        if len(self.multialg.tabmsa.sequence[0]) != len(self.reference):
            raise ValueError("Reference and aligned sequences must have "
                             "the same length.")

    @classmethod
    def from_fasta(cls,
                   sequences: str,
                   reference: Optional[str] = None):
        """Read a multialignment from a fasta file.

        If `reference` is not provided, it is assumed that the first
        sequence of the multialignment is the reference sequence.
        Otherwise, an additional fasta file with the reference sequence
        is needed.

        Args:
            sequences: input fasta file with multialignment
            reference: optional fasta file with reference sequence
        """
        msa = TabularMSA.read(sequences, constructor=DNA)

        if not reference:
            refer = msa[0]
            multialg = {seq.metadata.get("id"): str(seq) for seq in msa[1:]}
        else:
            refer = Sequence.read(reference)
            multialg = {seq.metadata.get("id"): str(seq) for seq in msa}

        ref = Reference(refer)
        alg = MultiAlignment(multialg)

        return cls(multialg=alg, reference=ref)

    @classmethod
    def from_csv(cls,
                 sequences: str,
                 reference: Optional[str] = None,
                 **kwargs):
        """Read a multialignment from a csv file.

        If `reference` is not provided, it is assumed that the first
        sequence of the multialignment is the reference sequence.
        Otherwise, an additional csv file with the reference sequence is
        needed. In both cases, the input csv file must be composed of
        two columns only, one for sequences ids and the other for the
        actual sequences; if not, you can provide additional options for
        pandas to restrict the number of columns read.

        Args:
            sequences: input csv file with multialignment
            reference: optional csv file with reference sequence
            **kwargs: additional options for pandas.read_csv()
        """
        msa = pd.read_csv(sequences, **kwargs)
        if msa.shape[1] != 2:
            raise ValueError("Please make sure the input only contains two "
                             "columns.")

        if not reference:
            refer = msa.iloc[0, 1]
            msa = msa.iloc[1:, :]
            multialg = dict(zip(msa.iloc[:, 0], msa.iloc[:, 1]))
        else:
            refer = pd.read_csv(reference, **kwargs)
            if refer.shape[1] != 2:
                raise ValueError("Please make sure the input only contains "
                                 "two columns.")
            refer = refer.iloc[0, 1]
            multialg = dict(zip(msa.iloc[:, 0], msa.iloc[:, 1]))

        ref = Reference(refer)
        alg = MultiAlignment(multialg)

        return cls(multialg=alg, reference=ref)

    @cached_property
    def df(self) -> pd.DataFrame:
        """Convert sequences to the proper dataframe for further allele
        frequency calculations."""
        df = self.multialg.tabmsa["sequence"].str.split("", expand=True)
        df.drop([0, df.columns[len(df.columns) - 1]], axis=1, inplace=True)
        df.fillna("-", inplace=True)  # avoidable
        df.index = self.multialg.tabmsa["id"]
        df.columns = self.reference.indexes

        return df

    @cached_property
    def frequencies(self) -> pd.DataFrame:
        """Calculate allele frequencies for the 4 basic nucleotides,
        gaps and other (non-canonical) nucleotides."""
        rows = []
        cols = ["A", "C", "G", "T", "gap", "oth"]
        idxs = []
        for col in self.df.columns:
            freqs = self.df[col].value_counts(normalize=True)
            freq_A = freqs.get("A", 0.0)
            freq_C = freqs.get("C", 0.0)
            freq_G = freqs.get("G", 0.0)
            freq_T = freqs.get("T", 0.0)
            freq_gap = freqs.get("-", 0.0)
            freq_oth = 1.0 - (freq_A + freq_C + freq_G + freq_T + freq_gap)
            rows.append({"A": freq_A, "C": freq_C, "G": freq_G, "T": freq_T,
                         "gap": freq_gap, "oth": freq_oth})
            idxs.append(col)
        var_df = pd.DataFrame(rows, columns=cols, index=idxs)
        var_df.reset_index(inplace=True)
        var_df.rename({"index": "position"}, axis=1, inplace=True)

        return var_df

    def to_csv(self, output_file: str = "all_freqs.csv"):
        """Write the resulting allele frequency dataframe to disk.

        Args:
            output_file: output file name
        """
        self.frequencies.to_csv(output_file, index=False)

    def __repr__(self):
        return "<{} ({} sequences, {} positions)>".format(
            self.__class__.__name__, len(self.multialg), len(self.reference)
        )
