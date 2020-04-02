=======
History
=======

0.1.0 (2019-07-08)
==================

* First release.

0.1.1 (2019-08-08)
------------------

* Read and process multialignments from fasta and csv files (Python module only).

0.1.2 (2019-10-17)
------------------

* Add tests with and without reference included in multialignments;
* Add tests with real datasets (coming from haplogroup-specific multialignments).

0.1.3 (2019-10-18)
------------------

* Add more detailed tests for real datasets;
* Implement more efficient frequency calculation;
* Add dunder methods and sanity checks;
* Fix requirements and testing framework;
* Clean code.

0.2.0 (2020-03-07)
==================

* Remove `numpy` and `pandas` from requirements as they are installed by `scikit-bio`;
* Move `tests` module inside `allfreqs`;
* Add `ci` module for internal management;
* Clean code.

0.3.0 (2020-04-02)
==================

* Add option to allow ambiguous nucleotides shown separately.
