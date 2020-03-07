#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os
import pytest

from skbio import Sequence
import pandas.testing as pdtest

from allfreqs.classes import MultiAlignment, Reference
from allfreqs.tests.constants import DATADIR


class TestMultiAlignment:

    def test_tabmsa(self,
                    sample_sequences_dict,
                    sample_sequences_tabmsa):
        m = MultiAlignment(sample_sequences_dict)
        pdtest.assert_frame_equal(m.tabmsa,
                                  sample_sequences_tabmsa)


class TestReference:

    def test_indexes_from_str(self,
                              sample_reference_indexes):
        ref = "AAG-CTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT"
        r = Reference(ref)
        assert r.indexes == sample_reference_indexes

    def test_indexes_from_sequence(self,
                                   sample_reference_indexes):
        ref = Sequence.read(os.path.join(DATADIR,
                                         "sample_ref.fasta"))
        r = Reference(ref)
        assert r.indexes == sample_reference_indexes