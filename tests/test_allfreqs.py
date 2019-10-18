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


class TestBasic:
    ref = Reference("AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT")
    d = {"seq1": "AAGGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT",
         "seq2": "AAGGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGATAT",
         "seq3": "ACC-GGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGCTAA",
         "seq4": "AAA-CCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACT-CAT",
         "seq5": "AAAGCCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCT-TAA"}
    alg = MultiAlignment(d)
    a = AlleleFreqs(multialg=alg, reference=ref)

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)


# From Fasta

class TestFromFasta:
    a = AlleleFreqs.from_fasta(
        sequences=os.path.join(DATADIR, "sample_multialg.fasta")
    )

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)


class TestFromFastaNoRef:
    a = AlleleFreqs.from_fasta(
        sequences=os.path.join(DATADIR, "sample_multialg_noref.fasta"),
        reference=os.path.join(DATADIR, "sample_ref.fasta")
    )

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)


# From Csv

class TestFromCsv:
    a = AlleleFreqs.from_csv(
        sequences=os.path.join(DATADIR, "sample_multialg.csv")
    )

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)


class TestFromCsvNoRef:
    a = AlleleFreqs.from_csv(
        sequences=os.path.join(DATADIR, "sample_multialg_noref.csv"),
        reference=os.path.join(DATADIR, "sample_ref.csv")
    )

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, sample_frequencies_csv)


# Real Datasets

class TestRealDatasetsX:
    a = AlleleFreqs.from_fasta(
        sequences=os.path.join(DATADIR, "real_datasets", "alg_X.fasta")
    )

    def test_df(self, alg_X_df):
        assert_frame_equal(self.a.df, alg_X_df)

    def test_freqs(self, alg_X_frequencies_csv):
        assert_frame_equal(self.a.frequencies, alg_X_frequencies_csv)

    def test_to_csv(self, alg_X_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, alg_X_frequencies_csv)


class TestRealDatasetsXNoRef:
    a = AlleleFreqs.from_fasta(
        sequences=os.path.join(DATADIR, "real_datasets", "alg_X_noref.fasta"),
        reference=os.path.join(DATADIR, "real_datasets", "RSRS.fasta")
    )

    def test_df(self, alg_X_df):
        assert_frame_equal(self.a.df, alg_X_df)

    def test_freqs(self, alg_X_frequencies_csv):
        assert_frame_equal(self.a.frequencies, alg_X_frequencies_csv)

    def test_to_csv(self, alg_X_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, alg_X_frequencies_csv)


class TestRealDatasetsL6:
    a = AlleleFreqs.from_fasta(
        sequences=os.path.join(DATADIR, "real_datasets", "alg_L6.fasta")
    )

    def test_df(self, alg_L6_df):
        assert_frame_equal(self.a.df, alg_L6_df)

    def test_freqs(self, alg_L6_frequencies_csv):
        assert_frame_equal(self.a.frequencies, alg_L6_frequencies_csv)

    def test_to_csv(self, alg_L6_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, alg_L6_frequencies_csv)


class TestRealDatasetsL6NoRef:
    a = AlleleFreqs.from_fasta(
        sequences=os.path.join(DATADIR, "real_datasets", "alg_L6_noref.fasta"),
        reference=os.path.join(DATADIR, "real_datasets", "RSRS.fasta"))

    def test_df(self, alg_L6_df):
        assert_frame_equal(self.a.df, alg_L6_df)

    def test_freqs(self, alg_L6_frequencies_csv):
        assert_frame_equal(self.a.frequencies, alg_L6_frequencies_csv)

    def test_to_csv(self, alg_L6_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        assert_frame_equal(result, alg_L6_frequencies_csv)

