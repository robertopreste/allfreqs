#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from typing import Dict, List, Union

from cached_property import cached_property
import pandas as pd
from skbio import Sequence


class MultiAlignment:
    """Class that creates a dataframe from a multialignment.

    Input is a dictionary of the form {"sequence id": "sequence"}, which
    is then used to create the resulting dataframe.
    """

    def __init__(self, msa: Dict[str, str]):
        self._msa = msa

    @cached_property
    def tabmsa(self) -> pd.DataFrame:
        """Create a dataframe with id and sequence columns.

        The resulting dataframe will be used to calculate allele frequencies.
        """
        df = pd.DataFrame.from_dict(self._msa,
                                    orient="index", columns=["sequence"])
        df.reset_index(inplace=True)
        df.rename({"index": "id"}, axis=1, inplace=True)

        return df

    def __len__(self):
        return self.tabmsa.shape[0]

    def __repr__(self):
        return repr(self.tabmsa)


class Reference:
    """Class to read and process a reference genome.

    Input is a reference genome in string format, which is used to create
    the reference positions.
    """

    def __init__(self, ref: Union[str, Sequence]):
        self._ref = str(ref)

    @cached_property
    def indexes(self) -> List[str]:
        """Create a set of indexes based on the reference genome.

        Each reference position will be as follows:
            1.0_G       G in position 1
            1.1_-       insertion in position 1 (determines a gap in ref)
        and so on.
        """
        indexes = []
        n = 0
        i = 0
        for nt in self._ref:
            if nt != "-":  # non-gap position
                i = 0
                n += 1
            else:  # gap position
                i += 1
            indexes.append("{}.{}_{}".format(n, i, nt))

        return indexes

    def __len__(self):
        return len(self._ref)

    def __repr__(self):
        return self._ref
