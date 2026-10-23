SHARP ASYMPTOTICS FOR QUADRATIC LAGRANGE INTERPOLATION

main.tex is the complete standalone amsart source. main.pdf is its compiled
A4 version, with one-inch margins, no author entry, and a five-sentence
abstract. The argument does not require supplementary code or unpublished
working notes to be included by LaTeX.

The theorem is M_n = 2 - C_*/n + o(1/n), where
C_* = 2 Gamma(5/4)^2 / Gamma(3/4)^2 = 1.0942198076132386...
The reference array consists of both endpoints and the zeros of the Jacobi
polynomial P_{n-2}^{(5/4,5/4)}. The paper proves an unrestricted lower
comparison for competitive arrays and an interior rigidity conclusion.

The coefficient-one conjecture was already disproved in the public forum.
The paper credits that development and the earlier numerical work. See
literature_notes.txt for the sources checked on 23 October 2026. The sharp
matching lower bound is the claimed contribution here. No statement of
catalogue acceptance, external referee approval, or formal verification is
made. Numerical outputs are consistency checks; the analytic proof is the
mathematical claim being presented for review.

To compile with a standard TeX installation, run twice from this folder:

  pdflatex -interaction=nonstopmode -halt-on-error main.tex

Alternatively run `tectonic --keep-logs main.tex`, which resolves cross-
references automatically. build.ps1 uses Tectonic and keeps its cache and
temporary files under this folder. Packages are declared in main.tex.

Supplementary calculations and a single-command runner are in verification/.
No external file, bibliography database, font download, or numerical output
is required after the declared LaTeX packages have been installed.

This directory is suitable for copying to a GitHub repository. Keep
main.tex, main.pdf, README.txt, literature_notes.txt, build.ps1, and the
verification directory. SHA256SUMS.txt records the prepared file hashes.
The TeX compiler binaries, caches, and the historical research archive are
not part of the publication bundle. Keep the existing relative paths when
uploading these files so that the verification commands and hashes work.
