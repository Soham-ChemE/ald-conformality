# Paper

`main.tex` is the arXiv-ready LaTeX source of the paper (article class, 11 pt, letter paper, manual `thebibliography`); it reads the figures from `../figures/`.

Build it with Tectonic from this directory: `cd paper && tectonic main.tex`. Tectonic downloads the required packages on the first run and writes `main.pdf` here.

`draft.md` is the Markdown source of record for the text; `main.tex` is a faithful typesetting of it and should be regenerated from it if the draft changes.

`assumptions.md` is the run log: every number in the paper traces to a run recorded there and to a script in the repository (see also `../CLAIMS.md`).

The figure PNGs are produced by the scripts in `../reactor/`, `../coupling/` and `../analysis/`; rerun those before rebuilding the PDF if any case is rerun.
