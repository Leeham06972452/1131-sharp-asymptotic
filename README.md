# Sharp asymptotics for Erdős problem 1131

Research manuscript and reproducible numerical checks for the sharp
first-order asymptotic in the quadratic Lagrange interpolation problem.

For distinct nodes $`Y=\{y_1,\ldots,y_n\}\subset[-1,1]`$, let
$`\ell_j^Y`$ be their Lagrange fundamental polynomials and define

```math
M_n=\inf_Y\int_{-1}^{1}\sum_{j=1}^{n}\bigl(\ell_j^Y(x)\bigr)^2\,dx.
```

The manuscript presents a proof of

```math
M_n=2-\frac{C_{\ast}}{n}+o(n^{-1}),\qquad
C_{\ast}=2\frac{\Gamma(5/4)^2}{\Gamma(3/4)^2}
=1.0942198076132386\ldots.
```

The proposed asymptotically optimal nodes are the endpoints together with
the zeros of the Jacobi polynomial $`P_{n-2}^{(5/4,5/4)}`$. The main
contribution claimed by the manuscript is the matching unrestricted lower
bound. The result concerns the first-order asymptotic, rather than the exact
minimum or the minimizing nodes for each finite $`n`$.

## Paper and supporting files

- [Read the paper (PDF)](main.pdf)
- [LaTeX source](main.tex)
- [Verification programs and instructions](verification/README.txt)
- [Recorded verification results](verification/results/summary.json)
- [Prior work and source notes](literature_notes.txt)
- [Build validation](build_validation.json) and [SHA-256 checksums](SHA256SUMS.txt)

## Reproduce the numerical checks

With Python installed, run these commands from the repository root:

```sh
python -m pip install -r verification/requirements.txt
python -B verification/verify.py
```

The supplied results were generated with Python 3.13.1, NumPy 2.4.2 and
SciPy 1.18.0. The runner writes its outputs to `verification/results/` and
checks the finite interpolation identities, exterior factorization, and an
independent quadrature comparison with the predicted limiting integral.
See [verification/README.txt](verification/README.txt) for the individual
programs and numerical tolerances.

The calculations are floating-point consistency checks. The analytic proof
is contained in the manuscript and does not assume their numerical results.
No proof-assistant certification or external referee approval is claimed.

## Build the PDF

With a standard LaTeX installation, run twice:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Alternatively, use `tectonic --keep-logs main.tex`. On Windows,
[`build.ps1`](build.ps1) runs Tectonic and keeps its cache and temporary files
inside the repository. The source uses `amsart`, A4 paper and one-inch
margins, with no external bibliography database or generated proof inputs.

## Prior work

The [problem discussion](https://www.erdosproblems.com/forum/thread/1131)
already records a disproof of the coefficient-one conjecture and earlier
numerical evidence for the value above. The manuscript credits those
developments and distinguishes them from its claim of a sharp unrestricted
coefficient. References and the scope of the literature check are recorded
in [literature_notes.txt](literature_notes.txt).
