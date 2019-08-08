#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
import os
from skbio import Sequence
from pandas.testing import assert_frame_equal
from allfreqs.classes import MultiAlignment, Reference

DATADIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


class TestMultiAlignment:

    def test_tabmsa(self, sample_sequences_dict, sample_sequences_tabmsa):
        m = MultiAlignment(sample_sequences_dict)
        assert_frame_equal(m.tabmsa, sample_sequences_tabmsa)


class TestReference:

    def test_indexes_from_str(self, sample_reference_indexes):
        ref = "AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT"
        r = Reference(ref)
        assert r.indexes == sample_reference_indexes

    def test_indexes_from_sequence(self, sample_reference_indexes):
        ref = Sequence.read(os.path.join(DATADIR, "sample_ref.fasta"))
        r = Reference(ref)
        assert r.indexes == sample_reference_indexes
