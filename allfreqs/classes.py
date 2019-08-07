#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from skbio import TabularMSA, DNA, Sequence
import pandas as pd
# import modin.pandas as pd


class MultiAlignment:

    def __init__(self):
        self.msa = None
        self.ref = None

    def from_fasta(self, input_file, with_ref=False):
        self.msa = TabularMSA.read(input_file, constructor=DNA)
        if with_ref:
            self.ref = self.msa[0]
            self.msa = self.msa[1:]

    def from_csv(self, input_file, with_ref=False, **kwargs):
        self.msa = pd.read_csv(input_file, **kwargs)
        if with_ref:
            self.ref = self.msa.iloc[0, :]
            self.msa = self.msa.iloc[1:, :]




class Multialignment:
    """Class to read and process a multialignment.

    It reads in a multisequence fasta file, then converts it to a dataframe
    with (id, sequence) columns.
    """

    def __init__(self, input_file):
        self.msa = input_file
        self.tabmsa = self.read()

    def read(self):
        """Read a multifasta using the TabularMSA class from skbio."""
        return TabularMSA.read(self.msa, constructor=DNA)

    @property
    def df(self):
        """Convert the multialignment to a dataframe with (id, sequence)
        columns."""
        df_dict = {seq.metadata.get("id"): str(seq) for seq in self.tabmsa}
        df = pd.DataFrame.from_dict(df_dict,
                                    orient="index", columns=["sequence"])
        df.reset_index(inplace=True)
        df.rename({"index": "id"}, axis=1, inplace=True)
        return df

    def write(self, output_file: str = "sequences_df.csv"):
        """Write the resulting dataframe to disk.

        Args:
            output_file: output file name
        """
        self.df.to_csv(output_file, index=False)


class Reference:
    """Class to read and process a reference genome.

    It reads a reference genome in fasta format, which is used to create
    the reference positions.
    """

    def __init__(self, input_file):
        self.input_file = input_file
        self.sequence = self.read()

    def read(self):
        """Read the input reference sequence."""
        return Sequence.read(self.input_file)

    @property
    def indexes(self):
        """Create a set of indexes based on the reference genome.

        Each reference position will be as follows:
            1.0_G       G in position 1
            1.1_-       insertion in position 1 (determines a gap in ref)
        and so on.
        """
        indexes = []
        n = 0
        i = 0
        for nt in str(self.sequence):
            if nt != "-":  # non-gap position
                i = 0
                n += 1
            else:  # gap position
                i += 1
            indexes.append("{}.{}_{}".format(n, i, nt))
        return indexes


class AlleleFreqs:
    """Class used to calculate allele frequencies from a multialignment.

    Input: Reference(), Multialignment().
    Output: dataframe with allele frequencies.
    """

    def __init__(self, reference, sequences):
        self.reference = Reference(reference)
        self.multialg = Multialignment(sequences)

    @property
    def df(self):
        """Convert sequences to the proper dataframe for further allele
        frequency calculations."""
        df = self.multialg.df["sequence"].str.split("", expand=True)
        df.drop([0, df.columns[len(df.columns) - 1]], axis=1, inplace=True)
        df.fillna("-", inplace=True)  # avoidable
        df.index = self.multialg.df["id"]
        df.columns = self.reference.indexes
        return df

    @property
    def frequencies(self):
        """Calculate allele frequencies for the 4 basic nucleotides,
        gaps and other (non-canonical) nucleotides."""
        var_df = pd.DataFrame(columns=["A", "C", "G", "T", "gap", "oth"])
        for col in self.df.columns:
            freqs = self.df[col].value_counts(normalize=True)
            freq_A = freqs.get("A", 0.0)
            freq_C = freqs.get("C", 0.0)
            freq_G = freqs.get("G", 0.0)
            freq_T = freqs.get("T", 0.0)
            freq_gap = freqs.get("-", 0.0)
            freq_oth = 1.0 - (freq_A + freq_C + freq_G + freq_T + freq_gap)
            var_df = pd.concat([var_df,
                                pd.DataFrame({"A": freq_A, "C": freq_C,
                                              "G": freq_G, "T": freq_T,
                                              "gap": freq_gap, "oth": freq_oth},
                                             index=[col])])
        var_df.reset_index(inplace=True)
        var_df.rename({"index": "position"}, axis=1, inplace=True)
        return var_df

    def write(self, output_file: str = "all_freqs.csv"):
        """Write the resulting allele frequency dataframe to disk.

        Args:
            output_file: output file name
        """
        self.frequencies.to_csv(output_file, index=False)
