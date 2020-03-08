#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import unittest

from skbio import Sequence
import pandas.testing as pdtest

from allfreqs.classes import MultiAlignment, Reference
from allfreqs.tests.constants import (
    SAMPLE_SEQUENCES_DICT, SAMPLE_SEQUENCES_TABMSA,
    SAMPLE_REF_FASTA, SAMPLE_REFERENCE_INDEXES
)


class TestMultiAlignment(unittest.TestCase):

    def setUp(self) -> None:
        self.multialg = MultiAlignment(SAMPLE_SEQUENCES_DICT)

    def test_tabmsa(self):
        pdtest.assert_frame_equal(self.multialg.tabmsa,
                                  SAMPLE_SEQUENCES_TABMSA)

    def test_length(self):
        self.assertEqual(5, len(self.multialg))


class TestReference(unittest.TestCase):

    def test_indexes_from_str(self):
        ref = "AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT"
        reference = Reference(ref)
        self.assertEqual(SAMPLE_REFERENCE_INDEXES, reference.indexes)

    def test_indexes_from_sequence(self):
        ref = Sequence.read(SAMPLE_REF_FASTA)
        reference = Reference(ref)
        self.assertEqual(SAMPLE_REFERENCE_INDEXES, reference.indexes)

    def test_length(self):
        ref = "AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT"
        reference = Reference(ref)
        self.assertEqual(44, len(reference))
