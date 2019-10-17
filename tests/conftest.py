#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
import os
import pandas as pd

DATADIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


@pytest.fixture
def sample_sequences_tabmsa():
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
def sample_sequences_dict():
    d = {"seq1": "AAGGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGG-TAT",
         "seq2": "AAGGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGATAT",
         "seq3": "ACC-GGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGCTAA",
         "seq4": "AAA-CCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACT-CAT",
         "seq5": "AAAGCCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCT-TAA"}
    return d


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


@pytest.fixture
def sample_sequences_df():
    df = pd.DataFrame.from_dict({
        "1.0_A": ["A", "A", "A", "A", "A"],
        "2.0_A": ["A", "A", "C", "A", "A"],
        "3.0_G": ["G", "G", "C", "A", "A"],
        "3.1_-": ["G", "G", "-", "-", "G"],
        "4.0_C": ["C", "C", "G", "C", "C"],
        "5.0_T": ["T", "C", "G", "C", "C"],
        "6.0_N": ["N", "T", "T", "C", "C"],
        "7.0_G": ["G", "T", "T", "T", "T"],
        "8.0_G": ["G", "G", "G", "T", "T"],
        "9.0_G": ["G", "G", "G", "G", "G"],
        "10.0_C": ["C", "C", "C", "C", "C"],
        "11.0_A": ["A", "A", "C", "C", "C"],
        "12.0_T": ["T", "G", "G", "G", "G"],
        "13.0_T": ["T", "T", "T", "T", "G"],
        "14.0_T": ["T", "G", "T", "T", "T"],
        "15.0_C": ["C", "C", "C", "A", "A"],
        "16.0_A": ["A", "A", "A", "C", "C"],
        "17.0_G": ["G", "G", "G", "G", "G"],
        "18.0_G": ["G", "G", "G", "C", "C"],
        "19.0_G": ["G", "G", "G", "T", "T"],
        "20.0_T": ["T", "T", "T", "T", "T"],
        "21.0_G": ["G", "G", "A", "A", "A"],
        "22.0_A": ["A", "A", "C", "A", "A"],
        "23.0_G": ["G", "G", "A", "A", "A"],
        "24.0_C": ["C", "C", "G", "C", "C"],
        "25.0_C": ["C", "C", "G", "C", "C"],
        "26.0_C": ["C", "G", "T", "G", "A"],
        "27.0_G": ["G", "T", "T", "A", "T"],
        "28.0_G": ["G", "G", "G", "G", "T"],
        "29.0_G": ["G", "G", "G", "G", "G"],
        "30.0_C": ["C", "C", "C", "C", "C"],
        "31.0_A": ["A", "C", "C", "C", "C"],
        "32.0_A": ["A", "G", "G", "G", "G"],
        "33.0_T": ["T", "G", "T", "G", "G"],
        "34.0_A": ["A", "G", "T", "G", "T"],
        "35.0_C": ["C", "C", "C", "A", "A"],
        "36.0_A": ["A", "A", "A", "C", "C"],
        "37.0_G": ["G", "C", "G", "A", "G"],
        "38.0_G": ["G", "G", "G", "C", "C"],
        "39.0_G": ["G", "G", "G", "T", "T"],
        "39.1_-": ["-", "A", "C", "-", "-"],
        "40.0_T": ["T", "T", "T", "C", "T"],
        "41.0_A": ["A", "A", "A", "A", "A"],
        "42.0_T": ["T", "T", "A", "T", "A"]
    })
    df.index = ["seq1", "seq2", "seq3", "seq4", "seq5"]
    df.index.name = "id"
    return df


@pytest.fixture
def sample_sequences_freqs():
    df = pd.DataFrame.from_dict(
        {0: {'position': '1.0_A', 'A': 1.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         1: {'position': '2.0_A', 'A': 0.8, 'C': 0.2, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         2: {'position': '3.0_G', 'A': 0.4, 'C': 0.2, 'G': 0.4, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         3: {'position': '3.1_-', 'A': 0.0, 'C': 0.0, 'G': 0.6, 'T': 0.0, 'gap': 0.4, 'oth': 0.0},
         4: {'position': '4.0_C', 'A': 0.0, 'C': 0.8, 'G': 0.2, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         5: {'position': '5.0_T', 'A': 0.0, 'C': 0.6, 'G': 0.2, 'T': 0.2, 'gap': 0.0, 'oth': 0.0},
         6: {'position': '6.0_N', 'A': 0.0, 'C': 0.4, 'G': 0.0, 'T': 0.4, 'gap': 0.0, 'oth': 0.2},
         7: {'position': '7.0_G', 'A': 0.0, 'C': 0.0, 'G': 0.2, 'T': 0.8, 'gap': 0.0, 'oth': 0.0},
         8: {'position': '8.0_G', 'A': 0.0, 'C': 0.0, 'G': 0.6, 'T': 0.4, 'gap': 0.0, 'oth': 0.0},
         9: {'position': '9.0_G', 'A': 0.0, 'C': 0.0, 'G': 1.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         10: {'position': '10.0_C', 'A': 0.0, 'C': 1.0, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         11: {'position': '11.0_A', 'A': 0.4, 'C': 0.6, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         12: {'position': '12.0_T', 'A': 0.0, 'C': 0.0, 'G': 0.8, 'T': 0.2, 'gap': 0.0, 'oth': 0.0},
         13: {'position': '13.0_T', 'A': 0.0, 'C': 0.0, 'G': 0.2, 'T': 0.8, 'gap': 0.0, 'oth': 0.0},
         14: {'position': '14.0_T', 'A': 0.0, 'C': 0.0, 'G': 0.2, 'T': 0.8, 'gap': 0.0, 'oth': 0.0},
         15: {'position': '15.0_C', 'A': 0.4, 'C': 0.6, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         16: {'position': '16.0_A', 'A': 0.6, 'C': 0.4, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         17: {'position': '17.0_G', 'A': 0.0, 'C': 0.0, 'G': 1.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         18: {'position': '18.0_G', 'A': 0.0, 'C': 0.4, 'G': 0.6, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         19: {'position': '19.0_G', 'A': 0.0, 'C': 0.0, 'G': 0.6, 'T': 0.4, 'gap': 0.0, 'oth': 0.0},
         20: {'position': '20.0_T', 'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 1.0, 'gap': 0.0, 'oth': 0.0},
         21: {'position': '21.0_G', 'A': 0.6, 'C': 0.0, 'G': 0.4, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         22: {'position': '22.0_A', 'A': 0.8, 'C': 0.2, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         23: {'position': '23.0_G', 'A': 0.6, 'C': 0.0, 'G': 0.4, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         24: {'position': '24.0_C', 'A': 0.0, 'C': 0.8, 'G': 0.2, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         25: {'position': '25.0_C', 'A': 0.0, 'C': 0.8, 'G': 0.2, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         26: {'position': '26.0_C', 'A': 0.2, 'C': 0.2, 'G': 0.4, 'T': 0.2, 'gap': 0.0, 'oth': 0.0},
         27: {'position': '27.0_G', 'A': 0.2, 'C': 0.0, 'G': 0.2, 'T': 0.6, 'gap': 0.0, 'oth': 0.0},
         28: {'position': '28.0_G', 'A': 0.0, 'C': 0.0, 'G': 0.8, 'T': 0.2, 'gap': 0.0, 'oth': 0.0},
         29: {'position': '29.0_G', 'A': 0.0, 'C': 0.0, 'G': 1.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         30: {'position': '30.0_C', 'A': 0.0, 'C': 1.0, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         31: {'position': '31.0_A', 'A': 0.2, 'C': 0.8, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         32: {'position': '32.0_A', 'A': 0.2, 'C': 0.0, 'G': 0.8, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         33: {'position': '33.0_T', 'A': 0.0, 'C': 0.0, 'G': 0.6, 'T': 0.4, 'gap': 0.0, 'oth': 0.0},
         34: {'position': '34.0_A', 'A': 0.2, 'C': 0.0, 'G': 0.4, 'T': 0.4, 'gap': 0.0, 'oth': 0.0},
         35: {'position': '35.0_C', 'A': 0.4, 'C': 0.6, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         36: {'position': '36.0_A', 'A': 0.6, 'C': 0.4, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         37: {'position': '37.0_G', 'A': 0.2, 'C': 0.2, 'G': 0.6, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         38: {'position': '38.0_G', 'A': 0.0, 'C': 0.4, 'G': 0.6, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         39: {'position': '39.0_G', 'A': 0.0, 'C': 0.0, 'G': 0.6, 'T': 0.4, 'gap': 0.0, 'oth': 0.0},
         40: {'position': '39.1_-', 'A': 0.2, 'C': 0.2, 'G': 0.0, 'T': 0.0, 'gap': 0.6, 'oth': 0.0},
         41: {'position': '40.0_T', 'A': 0.0, 'C': 0.2, 'G': 0.0, 'T': 0.8, 'gap': 0.0, 'oth': 0.0},
         42: {'position': '41.0_A', 'A': 1.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, 'gap': 0.0, 'oth': 0.0},
         43: {'position': '42.0_T', 'A': 0.4, 'C': 0.0, 'G': 0.0, 'T': 0.6, 'gap': 0.0, 'oth': 0.0}},
        orient="index")
    return df


@pytest.fixture
def sample_frequencies_csv():
    df = pd.read_csv(os.path.join(DATADIR, "sample_frequencies.csv"))
    return df


@pytest.fixture
def alg_X_frequencies_csv():
    df = pd.read_csv(os.path.join(DATADIR, "real_datasets", "frequencies_X.csv"))
    return df


@pytest.fixture
def alg_X_df():
    df = pd.read_csv(os.path.join(DATADIR, "real_datasets", "alg_X_df.csv"),
                     index_col=0)
    return df


@pytest.fixture
def alg_L6_frequencies_csv():
    df = pd.read_csv(os.path.join(DATADIR, "real_datasets", "frequencies_L6.csv"))
    return df


@pytest.fixture
def alg_L6_df():
    df = pd.read_csv(os.path.join(DATADIR, "real_datasets", "alg_L6_df.csv"),
                     index_col=0)
    return df
