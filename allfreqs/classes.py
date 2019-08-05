#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from skbio import TabularMSA, DNA, Sequence
import pandas as pd


class MultiAlignment:
    """Class that will store a multialignment.

    It will read in a multisequence fasta file, then convert it to a dataframe
    with (id, sequence) columns.
    """
    def __init__(self, input_msa):
        self.msa = input_msa
        self.tabmsa = self.read()
        # self.df = None

    def read(self):
        """Read a multifasta."""
        return TabularMSA.read(self.msa, constructor=DNA)

    @property
    def df(self):
        """Dump the multialignment to a dataframe with (id, sequence)
        columns."""
        df_dict = {seq.metadata.get("id"): str(seq) for seq in self.tabmsa}
        df = pd.DataFrame.from_dict(df_dict,
                                    orient="index", columns=["sequence"])
        df.reset_index(inplace=True)
        df.rename({"index": "id"}, axis=1, inplace=True)
        # self.df = df
        return df

    def write(self):
        self.df.to_csv("seq_df.csv", index=False)


class Reference:
    """Reference genome, will be used to create the reference positions."""
    def __init__(self, input_file):
        self.input_file = input_file
        self.sequence = self.read()
        # self.indexes = None

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
        # self.indexes = indexes
        return indexes


class AlignmentDataframe:
    """Input: Reference(), Multialignment().
    Output: dataframe with allele frequencies."""
    def __init__(self, reference, sequences):
        self.reference = Reference(reference)
        self.sequences = MultiAlignment(sequences)
        self.df = None
        self.var_df = None

    def to_df(self):
        df = self.sequences.df["sequence"].str.split("", expand=True)
        df.drop([0, df.columns[len(df.columns) - 1]], axis=1, inplace=True)
        df.fillna("-", inplace=True)  # avoidable
        df.index = self.sequences.df["id"]
        df.columns = self.reference.indexes
        self.df = df

    def calculate(self):
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
        self.var_df = var_df
