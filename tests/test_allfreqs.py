#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
from allfreqs.classes import Multialignment, Reference
import os
from pandas.testing import assert_frame_equal

DATADIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


class TestMultialignment:
    seq_path = os.path.join(DATADIR, "sample_multialg.fasta")
    m = Multialignment(seq_path)

    def test_msa(self):
        assert self.m.msa == self.seq_path

    def test_df(self, sample_sequences_df):
        assert_frame_equal(self.m.df, sample_sequences_df)


class TestReference:
    ref_path = os.path.join(DATADIR, "sample_ref.fasta")
    r = Reference(ref_path)

    def test_input_file(self):
        assert self.r.input_file == self.ref_path

    def test_indexes(self, sample_reference_indexes):
        assert self.r.indexes == sample_reference_indexes
