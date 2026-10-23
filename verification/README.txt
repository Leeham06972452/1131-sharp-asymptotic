SUPPLEMENTARY VERIFICATION FOR THE SHARP ASYMPTOTIC

From this directory run:

  python -B verify.py

The requirements file records the versions used for the supplied run.
Install these in your preferred Python environment. The programs have no
network dependencies or external data files. verify.py writes only inside
verification/results/ and disables Python bytecode generation.

erdos_1131_regularity_reproduce.py checks exact finite cardinal-mass,
complex rank-one commutator, exterior logarithmic transport, reciprocal
amplitude, and derivative-weighted interpolation identities at n=16,32,64,128.
It also compares associated polynomial boundary amplitudes at larger n.

erdos_1131_global_reproduce.py checks the reference mass matrices, relative
mass-deficit spectra, endpoint rows, near-collision costs, exterior positive
factorization, smooth-phase predictions, and stationary loop identities.
The stationary diagnostics are supplementary and are not used in main.tex.

independent_energy.py computes cardinal polynomials directly from node
differences and integrates their squares with n-point Gauss-Legendre
quadrature, whose order is exact for the degree in exact arithmetic. Its
node perturbation has phase 0.2 sin(theta), with predicted modulus change
eta(x)=0.2x. It compares n[I(Y)-I(X)] against the two-dimensional fractional
integral at quadrature orders 128,256,512 and n up to 1024. This evaluation
does not reuse the Legendre evaluation-matrix method in the other scripts.

The first two programs are retained from the supplied research package.
The third program and the runner were prepared independently during review.
Results are floating-point consistency checks. They do not certify global
inequalities, limiting arguments, or the unrestricted minimum. The proof
is contained in main.tex and has no computational hypothesis.

summary.json reports the largest observed identity errors and the final
independent smooth-phase comparison. Individual JSON and log files retain
the full data. The script runner fails if the stated numerical tolerances
are exceeded; tolerances are numerical diagnostics, not theorem constants.
