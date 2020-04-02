=====
Usage
=====

Command Line Interface
======================

allfreqs can be used as a command line tool, using the ``allfreqs`` command and providing the input
fasta or csv file with multialigned sequences:

.. code-block:: console

    # multialignment in fasta format
    $ allfreqs multialg_seqs.fasta
    # if reference is stored separately:
    $ allfreqs multialg_seqs.fasta --reference my_ref.fasta

    # multialignment in csv format
    # e.g. seq1,ACGTACGT
    #      seq2,A-CTAGGT
    $ allfreqs multialg_seqs.csv
    # if reference is stored separately:
    $ allfreqs multialg_seqs.csv --reference my_ref.csv

The program will use the first sequence in the multialignment as the reference sequence; if this is
not the case, you can supply a reference sequence using the ``--reference|-r`` option followed by
the fasta or csv file with the desired reference sequence to use. **Please note that in this case
both multialigned sequences and reference sequence must be in the same format (both fasta or csv
files).**

By default, allfreqs will add frequencies of non-standard (ambiguous) nucleotides together, showing
them in the ``oth`` column of the output; it is possible to show them in separate columns, specific
for each of them, using the ``--ambiguous`` flag.

allfreqs will calculate allele frequencies for each position in the multialignment and save them as
a csv file called ``all_freqs.csv`` in the current working directory. It is possible to specify a
different output location using the ``--out|-o`` option followed by the desired path/filename.

____

Python Module
=============

allfreqs can be used in a Python script by importing its ``AlleleFreqs`` class.

This class has two methods, ``.from_fasta()`` and ``.from_csv()``, which can be used to load
multialignments from either fasta or csv files respectively. Both methods accept a mandatory
``sequences`` argument, which specifies the file containing multialigned sequences, and an optional
``reference`` argument, which can be used to specify the reference sequence in case it is not
reported as the first sequence in the provided ``sequences`` file:

.. code-block:: python

    from allfreqs import AlleleFreqs

    # multialignment in fasta format
    a = AlleleFreqs.from_fasta(sequences="multialg_seqs.fasta")
    # if reference is stored separately:
    a = AlleleFreqs.from_fasta(sequences="multialg_seqs.fasta", reference="my_ref.fasta")

    # multialignment in csv format
    a = AlleleFreqs.from_csv(sequences="multialg_seqs.csv")
    # if reference is stored separately:
    a = AlleleFreqs.from_csv(sequences="multialg_seqs.csv", reference="my_ref.csv")

The ``AlleleFreqs`` class has two useful properties:

- ``df``, which returns a dataframe with sequences as rows and single positions as columns;
- ``frequencies``, which returns a dataframe with the actual allele frequencies for each position.

By default, allfreqs will add frequencies of non-standard (ambiguous) nucleotides together, showing
them in the ``oth`` column of the output; it is possible to show them in separate columns, specific
for each of them, using the ``ambiguous=True`` option.

These allele frequencies can be saved to a csv file using the ``.to_csv()`` method; by default they
will be saved to a file called ``all_freqs.csv`` in the current working directory, but this can be
overridden by providing the ``output_file`` argument to ``.to_csv()``.
