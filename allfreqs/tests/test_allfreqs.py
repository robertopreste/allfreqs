#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import unittest

import pandas as pd
import pandas.testing as pdtest

from allfreqs import AlleleFreqs
from allfreqs.classes import Reference, MultiAlignment
from allfreqs.tests.constants import (
    REAL_ALG_X_FASTA, REAL_ALG_X_NOREF_FASTA, REAL_RSRS_FASTA,
    REAL_ALG_L6_FASTA, REAL_ALG_L6_NOREF_FASTA,
    SAMPLE_MULTIALG_FASTA, SAMPLE_MULTIALG_NOREF_FASTA, SAMPLE_REF_FASTA,
    SAMPLE_MULTIALG_CSV, SAMPLE_MULTIALG_NOREF_CSV, SAMPLE_REF_CSV,
    SAMPLE_SEQUENCES_DF, SAMPLE_SEQUENCES_FREQS, SAMPLE_FREQUENCIES_CSV,
    ALG_X_DF, ALG_X_FREQUENCIES_CSV, ALG_L6_DF, ALG_L6_FREQUENCIES_CSV,
    TEST_CSV
)


class TestBasic(unittest.TestCase):

    def setUp(self) -> None:
        ref = Reference("AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT")
        d = {"seq1": "AAGGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT",
             "seq2": "AAGGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGATAT",
             "seq3": "ACC-GGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGCTAA",
             "seq4": "AAA-CCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACT-CAT",
             "seq5": "AAAGCCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCT-TAA"}
        alg = MultiAlignment(d)

        self.af = AlleleFreqs(multialg=alg, reference=ref)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, SAMPLE_SEQUENCES_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, SAMPLE_SEQUENCES_FREQS)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, SAMPLE_FREQUENCIES_CSV)


# From Fasta

class TestFromFasta(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=SAMPLE_MULTIALG_FASTA)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, SAMPLE_SEQUENCES_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, SAMPLE_SEQUENCES_FREQS)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, SAMPLE_FREQUENCIES_CSV)


class TestFromFastaNoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=SAMPLE_MULTIALG_NOREF_FASTA,
                                         reference=SAMPLE_REF_FASTA)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, SAMPLE_SEQUENCES_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, SAMPLE_SEQUENCES_FREQS)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, SAMPLE_FREQUENCIES_CSV)


# From Csv

class TestFromCsv(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_csv(sequences=SAMPLE_MULTIALG_CSV)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, SAMPLE_SEQUENCES_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, SAMPLE_SEQUENCES_FREQS)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, SAMPLE_FREQUENCIES_CSV)


class TestFromCsvNoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_csv(sequences=SAMPLE_MULTIALG_NOREF_CSV,
                                       reference=SAMPLE_REF_CSV)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, SAMPLE_SEQUENCES_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, SAMPLE_SEQUENCES_FREQS)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, SAMPLE_FREQUENCIES_CSV)


# Real Datasets

class TestRealDatasetsX(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_X_FASTA)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, ALG_X_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, ALG_X_FREQUENCIES_CSV)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, ALG_X_FREQUENCIES_CSV)


class TestRealDatasetsXNoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_X_NOREF_FASTA,
                                         reference=REAL_RSRS_FASTA)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, ALG_X_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, ALG_X_FREQUENCIES_CSV)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, ALG_X_FREQUENCIES_CSV)


class TestRealDatasetsL6(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_L6_FASTA)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, ALG_L6_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, ALG_L6_FREQUENCIES_CSV)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, ALG_L6_FREQUENCIES_CSV)


class TestRealDatasetsL6NoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_L6_NOREF_FASTA,
                                         reference=REAL_RSRS_FASTA)

    def test_df(self):
        pdtest.assert_frame_equal(self.af.df, ALG_L6_DF)

    def test_frequencies(self):
        pdtest.assert_frame_equal(self.af.frequencies, ALG_L6_FREQUENCIES_CSV)

    def test_to_csv(self):
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        pdtest.assert_frame_equal(result, ALG_L6_FREQUENCIES_CSV)
