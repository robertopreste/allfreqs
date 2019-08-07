#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
import pandas as pd


@pytest.fixture
def sample_sequences_df():
    df = pd.DataFrame({
        "id": ["seq1", "seq2", "seq3", "seq4", "seq5"],
        "sequence": ["AAGGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT",
                     "AAGGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGATAT",
                     "ACC-GGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGCTAA",
                     "AAA-CCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACT-CAT",
                     "AAAGCCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCT-TAA"]
    })
    return df


@pytest.fixture
def sample_reference_indexes():
    indexes = ["1.0_A", "2.0_A", "3.0_G", "3.1_-", "4.0_C", "5.0_T", "6.0_N",
               "7.0_G", "8.0_G", "9.0_G", "10.0_C", "11.0_A", "12.0_T",
               "13.0_T", "14.0_T", "15.0_C", "16.0_A", "17.0_G", "18.0_G",
               "19.0_G", "20.0_T", "21.0_G", "22.0_A", "23.0_G", "24.0_C",
               "25.0_C", "26.0_C", "27.0_G", "28.0_G", "29.0_G", "30.0_C",
               "31.0_A", "32.0_A", "33.0_T", "34.0_A", "35.0_C", "36.0_A",
               "37.0_G", "38.0_G", "39.0_G", "39.1_-", "40.0_T", "41.0_A",
               "42.0_T"]
    return indexes

