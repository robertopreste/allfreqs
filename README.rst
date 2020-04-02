========
allfreqs
========


.. image:: https://img.shields.io/pypi/v/allfreqs.svg
        :target: https://pypi.python.org/pypi/allfreqs

.. image:: https://www.repostatus.org/badges/latest/wip.svg
    :alt: Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.
    :target: https://www.repostatus.org/#wip

.. image:: https://travis-ci.com/robertopreste/allfreqs.svg?branch=master
        :target: https://travis-ci.com/robertopreste/allfreqs

.. image:: https://readthedocs.org/projects/allfreqs/badge/?version=latest
        :target: https://allfreqs.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Calculate allele frequencies from a sequence multialignment.


* Free software: MIT license
* Documentation: https://allfreqs.readthedocs.io
* GitHub repo: https://github.com/robertopreste/allfreqs


Features
========

Calculate allele frequencies from a nucleotide multialignment in fasta or csv format.

Allele frequencies will be returned as a table in which each row is a nucleotide position (based on
the provided reference sequence) and columns are A, C, G, T frequencies as well as gaps and other
non-canonical nucleotides.

For example, given the following multialignment:

+------+----------+
| ID   | Sequence |
+======+==========+
| ref  | ACGTACGT |
+------+----------+
| seq1 | A-GTAGGN |
+------+----------+
| seq2 | ACCAGCGT |
+------+----------+

the resulting allele frequencies will be:

+----------+-----+-----+-----+-----+-----+-----+
| position | A   | C   | G   | T   | gap | oth |
+==========+=====+=====+=====+=====+=====+=====+
| 1.0_A    | 1.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
+----------+-----+-----+-----+-----+-----+-----+
| 2.0_C    | 0.0 | 0.5 | 0.0 | 0.0 | 0.5 | 0.0 |
+----------+-----+-----+-----+-----+-----+-----+
| 3.0_G    | 0.0 | 0.5 | 0.5 | 0.0 | 0.0 | 0.0 |
+----------+-----+-----+-----+-----+-----+-----+
| 4.0_T    | 0.5 | 0.0 | 0.0 | 0.5 | 0.0 | 0.0 |
+----------+-----+-----+-----+-----+-----+-----+
| 5.0_A    | 0.5 | 0.0 | 0.5 | 0.0 | 0.0 | 0.0 |
+----------+-----+-----+-----+-----+-----+-----+
| 6.0_C    | 0.0 | 0.5 | 0.5 | 0.0 | 0.0 | 0.0 |
+----------+-----+-----+-----+-----+-----+-----+
| 7.0_G    | 0.0 | 0.0 | 1.0 | 0.0 | 0.0 | 0.0 |
+----------+-----+-----+-----+-----+-----+-----+
| 8.0_T    | 0.0 | 0.0 | 0.0 | 0.5 | 0.0 | 0.5 |
+----------+-----+-----+-----+-----+-----+-----+

Frequencies of non-canonical (ambiguous) nucleotides are by default squashed into the ``oth``
column, but they can also be shown separately using a simple flag.

allfreqs can be used either as a command line tool or through its Python API.

For more information, please refer to the Usage_ section of the documentation.

Installation
============

**PLEASE NOTE: allfreqs only supports Python >= 3.6!**

The preferred installation method for allfreqs is using ``pip``:

.. code-block:: console

    $ pip install allfreqs

For more information, please refer to the Installation_ section of the documentation.

Credits
=======

This package was created with Cookiecutter_ and the `cc-pypackage`_ project template.

.. _Usage: https://allfreqs.readthedocs.io/en/latest/usage.html
.. _Installation: https://allfreqs.readthedocs.io/en/latest/installation.html
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cc-pypackage`: https://github.com/robertopreste/cc-pypackage
