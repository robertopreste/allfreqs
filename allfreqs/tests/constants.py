#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os


DATADIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
TEST_CSV = os.path.join(DATADIR, "test.csv")

SAMPLE_MULTIALG_FASTA = os.path.join(DATADIR, "sample_multialg.fasta")
SAMPLE_MULTIALG_NOREF_FASTA = os.path.join(DATADIR,
                                           "sample_multialg_noref.fasta")
SAMPLE_REF_FASTA = os.path.join(DATADIR, "sample_ref.fasta")

SAMPLE_MULTIALG_CSV = os.path.join(DATADIR, "sample_multialg.csv")
SAMPLE_MULTIALG_NOREF_CSV = os.path.join(DATADIR, "sample_multialg_noref.csv")
SAMPLE_REF_CSV = os.path.join(DATADIR, "sample_ref.csv")

SAMPLE_FREQUENCIES = os.path.join(DATADIR, "sample_frequencies.csv")

REAL_ALG_X_FASTA = os.path.join(DATADIR, "real_datasets", "alg_X.fasta")
REAL_ALG_X_NOREF_FASTA = os.path.join(DATADIR,
                                      "real_datasets", "alg_X_noref.fasta")
REAL_ALG_X_DF = os.path.join(DATADIR, "real_datasets", "alg_X_df.csv")

REAL_ALG_L6_FASTA = os.path.join(DATADIR, "real_datasets", "alg_L6.fasta")
REAL_ALG_L6_NOREF_FASTA = os.path.join(DATADIR,
                                       "real_datasets", "alg_L6_noref.fasta")
REAL_ALG_L6_DF = os.path.join(DATADIR, "real_datasets", "alg_L6_df.csv")

REAL_RSRS_FASTA = os.path.join(DATADIR, "real_datasets", "RSRS.fasta")

REAL_X_FREQUENCIES = os.path.join(DATADIR,
                                  "real_datasets", "frequencies_X.csv")

REAL_L6_FREQUENCIES = os.path.join(DATADIR,
                                   "real_datasets", "frequencies_L6.csv")
