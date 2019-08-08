#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
import os
import pandas as pd
from pandas.testing import assert_frame_equal
from allfreqs import AlleleFreqs
from allfreqs.classes import Reference, MultiAlignment

DATADIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
TEST_CSV = os.path.join(DATADIR, "test.csv")


class TestAlleleFreqsBasic:
    ref = Reference("AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT")
    d = {"seq1": "AAGGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT",
         "seq2": "AAGGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGATAT",
         "seq3": "ACC-GGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGCTAA",
         "seq4": "AAA-CCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACT-CAT",
         "seq5": "AAAGCCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCT-TAA"}
    alg = MultiAlignment(d)
    a = AlleleFreqs(reference=ref, multialg=alg)

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_frequencies(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)


class TestAlleleFreqsFasta:
    a = AlleleFreqs.from_fasta(sequences=os.path.join(DATADIR, "sample_multialg.fasta"),
                               reference=os.path.join(DATADIR, "sample_ref.fasta"))

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)


class TestAlleleFreqsCsv:
    a = AlleleFreqs.from_csv(sequences=os.path.join(DATADIR, "sample_multialg.csv"),
                             reference=os.path.join(DATADIR, "sample_ref.csv"))

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)
