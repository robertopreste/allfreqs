#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os
import unittest

from click.testing import CliRunner
import pandas as pd
import pandas.testing as pdtest

from allfreqs import cli
from allfreqs.tests.constants import (
    SAMPLE_MULTIALG_CSV, SAMPLE_MULTIALG_NOREF_CSV, SAMPLE_REF_CSV,
    SAMPLE_MULTIALG_FASTA, SAMPLE_MULTIALG_NOREF_FASTA, SAMPLE_REF_FASTA,
    SAMPLE_FREQUENCIES, SAMPLE_FREQUENCIES_AMB, TEST_CSV
)


class TestFromFastaCLI(unittest.TestCase):

    def setUp(self) -> None:
        self.runner = CliRunner()

    def test_cli_help(self):
        # Given/When
        result = self.runner.invoke(cli.main, ["--help"])

        # Then
        self.assertEqual(0, result.exit_code)
        self.assertIn("Show this message and exit.", result.output)

    def test_basic(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_FASTA])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES)
        test_csv = pd.read_csv("all_freqs.csv")

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove("all_freqs.csv")

    def test_output_name(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_FASTA,
                                     "--out", TEST_CSV])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES)
        test_csv = pd.read_csv(TEST_CSV)

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove(TEST_CSV)

    def test_ambiguous(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_FASTA,
                                     "--ambiguous"])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES_AMB)
        test_csv = pd.read_csv("all_freqs.csv")

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove("all_freqs.csv")

    def test_no_ref(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_NOREF_FASTA,
                                     "--reference", SAMPLE_REF_FASTA])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES)
        test_csv = pd.read_csv("all_freqs.csv")

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove("all_freqs.csv")


class TestFromCsvCLI(unittest.TestCase):

    def setUp(self) -> None:
        self.runner = CliRunner()

    def test_basic(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_CSV])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES)
        test_csv = pd.read_csv("all_freqs.csv")

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove("all_freqs.csv")

    def test_output_name(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_CSV,
                                     "--out", TEST_CSV])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES)
        test_csv = pd.read_csv(TEST_CSV)

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove(TEST_CSV)

    def test_ambiguous(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_CSV,
                                     "--ambiguous"])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES_AMB)
        test_csv = pd.read_csv("all_freqs.csv")

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove("all_freqs.csv")

    def test_no_ref(self):
        # Given/When
        result = self.runner.invoke(cli.main,
                                    [SAMPLE_MULTIALG_NOREF_CSV,
                                     "--reference", SAMPLE_REF_CSV])
        exp_csv = pd.read_csv(SAMPLE_FREQUENCIES)
        test_csv = pd.read_csv("all_freqs.csv")

        # Then
        self.assertEqual(0, result.exit_code)
        pdtest.assert_frame_equal(test_csv, exp_csv)
        # Cleanup
        os.remove("all_freqs.csv")
