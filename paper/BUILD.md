# Building the paper

The manuscript uses AMS-LaTeX, Latin Modern Type 1 fonts, geometry, hyperref, and the amsplain BibTeX style. From this directory, run the following with a locally available Tectonic bundle.

~~~sh
SOURCE_DATE_EPOCH=1795305600 tectonic --only-cached --reruns 2 manuscript.tex
~~~

For installations using a separate existing cache, set TECTONIC_CACHE_DIR to that cache before running the command. Tectonic runs BibTeX automatically. Alternatively, with a complete TeX Live installation, run pdflatex manuscript, bibtex manuscript, then pdflatex manuscript twice.

The release PDF was built with Tectonic 0.17.0, double-checked in separate output directories, rendered page by page, and inspected visually. The original reviewed paper and the integer programs are preserved separately in [proof](../proof/). The instructions for replaying the Lean arithmetic and coefficient lemmas are in the [formalization README](../formalization/README.md).
