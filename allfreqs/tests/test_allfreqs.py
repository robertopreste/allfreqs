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
    sample_sequences_df, SAMPLE_SEQUENCES_DICT, sample_sequences_freqs,
    sample_sequences_freqs_amb, SAMPLE_FREQUENCIES,
    SAMPLE_FREQUENCIES_AMB, REAL_ALG_X_DF, REAL_X_FREQUENCIES, REAL_ALG_L6_DF,
    REAL_L6_FREQUENCIES, TEST_CSV
)


class TestBasic(unittest.TestCase):

    def setUp(self) -> None:
        ref = Reference("AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT")
        alg = MultiAlignment(SAMPLE_SEQUENCES_DICT)

        self.af = AlleleFreqs(multialg=alg, reference=ref)
        self.af_amb = AlleleFreqs(multialg=alg, reference=ref, ambiguous=True)

    def test_df(self):
        # Given/When
        exp_df = sample_sequences_df()

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = sample_sequences_freqs()

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_frequencies_ambiguous(self):
        # Given/When
        exp_freqs = sample_sequences_freqs_amb()

        # Then
        pdtest.assert_frame_equal(self.af_amb.frequencies, exp_freqs)

    def test__get_frequencies(self):
        # Given
        test_freq = pd.Series({'A': 0.2, 'C': 0.2, 'G': 0.1, 'T': 0.3,
                               '-': 0.1, 'N': 0.1})
        exp_freq = {'A': 0.2, 'C': 0.2, 'G': 0.1, 'T': 0.3, 'gap': 0.1,
                    'oth': 0.1}

        # When
        result = self.af._get_frequencies(test_freq)

        # Then
        self._dict_almost_equal(result, exp_freq)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(SAMPLE_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)

    def test_to_csv_ambiguous(self):
        # Given/When
        self.af_amb.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(SAMPLE_FREQUENCIES_AMB)

        # Then
        pdtest.assert_frame_equal(result, expected)

    @staticmethod
    def _dict_almost_equal(expected: dict, result: dict, acc=10**-8) -> bool:
        """Compare to dictionaries and ensure that all their values are the
        same, accounting for some fluctuation up to the given accuracy value.

        Args:
            expected: expected dictionary
            result: resulting dictionary
            acc: accuracy to use [default: 10**-8]
        """
        if expected.keys() == result.keys():
            for key in expected.keys():
                if abs(expected[key] - result[key]) < acc:
                    continue
            return True
        return False


# From Fasta

class TestFromFasta(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=SAMPLE_MULTIALG_FASTA)

    def test_df(self):
        # Given/When
        exp_df = sample_sequences_df()

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = sample_sequences_freqs()

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(SAMPLE_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)


class TestFromFastaNoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=SAMPLE_MULTIALG_NOREF_FASTA,
                                         reference=SAMPLE_REF_FASTA)

    def test_df(self):
        # Given/When
        exp_df = sample_sequences_df()

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = sample_sequences_freqs()

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(SAMPLE_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)


# From Csv

class TestFromCsv(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_csv(sequences=SAMPLE_MULTIALG_CSV)

    def test_df(self):
        # Given/When
        exp_df = sample_sequences_df()

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = sample_sequences_freqs()

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(SAMPLE_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)


class TestFromCsvNoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_csv(sequences=SAMPLE_MULTIALG_NOREF_CSV,
                                       reference=SAMPLE_REF_CSV)

    def test_df(self):
        # Given/When
        exp_df = sample_sequences_df()

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = sample_sequences_freqs()

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(SAMPLE_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)


# Real Datasets

class TestRealDatasetsX(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_X_FASTA)

    def test_df(self):
        # Given/When
        exp_df = pd.read_csv(REAL_ALG_X_DF, index_col=0)

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = pd.read_csv(REAL_X_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/when
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(REAL_X_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)


class TestRealDatasetsXNoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_X_NOREF_FASTA,
                                         reference=REAL_RSRS_FASTA)

    def test_df(self):
        # Given/When
        exp_df = pd.read_csv(REAL_ALG_X_DF, index_col=0)

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = pd.read_csv(REAL_X_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(REAL_X_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)


class TestRealDatasetsL6(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_L6_FASTA)

    def test_df(self):
        # Given/When
        exp_df = pd.read_csv(REAL_ALG_L6_DF, index_col=0)

        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = pd.read_csv(REAL_L6_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(REAL_L6_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)


class TestRealDatasetsL6NoRef(unittest.TestCase):

    def setUp(self) -> None:
        self.af = AlleleFreqs.from_fasta(sequences=REAL_ALG_L6_NOREF_FASTA,
                                         reference=REAL_RSRS_FASTA)

    def test_df(self):
        # Given/When
        exp_df = pd.read_csv(REAL_ALG_L6_DF, index_col=0)

        # Then
        pdtest.assert_frame_equal(self.af.df, exp_df)

    def test_frequencies(self):
        # Given/When
        exp_freqs = pd.read_csv(REAL_L6_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(self.af.frequencies, exp_freqs)

    def test_to_csv(self):
        # Given/When
        self.af.to_csv(TEST_CSV)
        result = pd.read_csv(TEST_CSV)
        expected = pd.read_csv(REAL_L6_FREQUENCIES)

        # Then
        pdtest.assert_frame_equal(result, expected)
