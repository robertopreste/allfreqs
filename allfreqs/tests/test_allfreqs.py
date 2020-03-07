#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os
import pytest

import pandas as pd
import pandas.testing as pdtest

from allfreqs import AlleleFreqs
from allfreqs.classes import Reference, MultiAlignment
from allfreqs.tests.constants import (
    REAL_ALG_X_FASTA, REAL_ALG_X_NOREF_FASTA, REAL_RSRS_FASTA,
    REAL_ALG_L6_FASTA, REAL_ALG_L6_NOREF_FASTA,
    SAMPLE_MULTIALG_FASTA, SAMPLE_MULTIALG_NOREF_FASTA, SAMPLE_REF_FASTA,
    SAMPLE_MULTIALG_CSV, SAMPLE_MULTIALG_NOREF_CSV, SAMPLE_REF_CSV,
    TEST_CSV
)


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
        pdtest.assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        pdtest.assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, sample_frequencies_csv)


# From Fasta

class TestFromFasta:
    a = AlleleFreqs.from_fasta(sequences=SAMPLE_MULTIALG_FASTA)

    def test_df(self, sample_sequences_df):
        pdtest.assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        pdtest.assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, sample_frequencies_csv)


class TestFromFastaNoRef:
    a = AlleleFreqs.from_fasta(sequences=SAMPLE_MULTIALG_NOREF_FASTA,
                               reference=SAMPLE_REF_FASTA)

    def test_df(self, sample_sequences_df):
        pdtest.assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        pdtest.assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, sample_frequencies_csv)


# From Csv

class TestFromCsv:
    a = AlleleFreqs.from_csv(sequences=SAMPLE_MULTIALG_CSV)

    def test_df(self, sample_sequences_df):
        pdtest.assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        pdtest.assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, sample_frequencies_csv)


class TestFromCsvNoRef:
    a = AlleleFreqs.from_csv(sequences=SAMPLE_MULTIALG_NOREF_CSV,
                             reference=SAMPLE_REF_CSV)

    def test_df(self, sample_sequences_df):
        pdtest.assert_frame_equal(self.a.df, sample_sequences_df)

    def test_freqs(self, sample_sequences_freqs):
        pdtest.assert_frame_equal(self.a.frequencies, sample_sequences_freqs)

    def test_to_csv(self, sample_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, sample_frequencies_csv)


# Real Datasets

class TestRealDatasetsX:
    a = AlleleFreqs.from_fasta(sequences=REAL_ALG_X_FASTA)

    def test_df(self, alg_X_df):
        pdtest.assert_frame_equal(self.a.df, alg_X_df)

    def test_freqs(self, alg_X_frequencies_csv):
        pdtest.assert_frame_equal(self.a.frequencies, alg_X_frequencies_csv)

    def test_to_csv(self, alg_X_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, alg_X_frequencies_csv)


class TestRealDatasetsXNoRef:
    a = AlleleFreqs.from_fasta(sequences=REAL_ALG_X_NOREF_FASTA,
                               reference=REAL_RSRS_FASTA)

    def test_df(self, alg_X_df):
        pdtest.assert_frame_equal(self.a.df, alg_X_df)

    def test_freqs(self, alg_X_frequencies_csv):
        pdtest.assert_frame_equal(self.a.frequencies, alg_X_frequencies_csv)

    def test_to_csv(self, alg_X_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, alg_X_frequencies_csv)


class TestRealDatasetsL6:
    a = AlleleFreqs.from_fasta(sequences=REAL_ALG_L6_FASTA)

    def test_df(self, alg_L6_df):
        pdtest.assert_frame_equal(self.a.df, alg_L6_df)

    def test_freqs(self, alg_L6_frequencies_csv):
        pdtest.assert_frame_equal(self.a.frequencies, alg_L6_frequencies_csv)

    def test_to_csv(self, alg_L6_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, alg_L6_frequencies_csv)


class TestRealDatasetsL6NoRef:
    a = AlleleFreqs.from_fasta(sequences=REAL_ALG_L6_NOREF_FASTA,
                               reference=REAL_RSRS_FASTA)

    def test_df(self, alg_L6_df):
        pdtest.assert_frame_equal(self.a.df, alg_L6_df)

    def test_freqs(self, alg_L6_frequencies_csv):
        pdtest.assert_frame_equal(self.a.frequencies, alg_L6_frequencies_csv)

    def test_to_csv(self, alg_L6_frequencies_csv):
        self.a.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, alg_L6_frequencies_csv)
