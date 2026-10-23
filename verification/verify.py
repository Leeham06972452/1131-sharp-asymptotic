"""Run the supplementary checks. All generated files stay in this directory."""
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["TEMP"] = os.environ["TMP"] = str(RESULTS)

for name, output in [
    ("erdos_1131_regularity_reproduce.py", "regularity.json"),
    ("erdos_1131_global_reproduce.py", "global.json"),
    ("independent_energy.py", None),
]:
    command = [sys.executable, "-B", str(ROOT / name)]
    if output:
        command += ["--output", str(RESULTS / output)]
    completed = subprocess.run(command, cwd=ROOT, check=True, text=True,
                               encoding="utf-8", capture_output=True)
    (RESULTS / (Path(name).stem + ".log")).write_text(
        completed.stdout + completed.stderr, encoding="utf-8")
    print(name + ": completed", flush=True)

regularity = json.loads((RESULTS / "regularity.json").read_text())
global_checks = json.loads((RESULTS / "global.json").read_text())
energy = json.loads((RESULTS / "independent_energy.json").read_text())
cases = regularity["finite_identities"]
errors = {
    "mass_identity": max(r["mass_identity_max_error"] for r in cases),
    "rank_one_commutator": max(r["rank_one_commutator_max_error"] for r in cases),
    "exterior_transport": max(abs(r["transport_error"]) for r in cases),
    "reciprocal_identity": max(abs(r["reciprocal_prufer_error"]) for r in cases),
    "weighted_interpolation": max(r["derivative_weighted_interpolation_max_error"]
                                  for r in cases),
    "exterior_factorization": max(abs(r["difference"])
        for r in global_checks["exterior_factorization_checks"]),
}
for label, error in errors.items():
    tolerance = 2e-8 if label == "exterior_factorization" else 1e-9
    if error >= tolerance:
        raise AssertionError(f"{label}: {error} exceeds {tolerance}")

prediction = energy["continuum_defect"][-1]["defect"]
finite = energy["finite_objectives"][-1]["scaled_excess"]
if abs(finite - prediction) >= 5e-5:
    raise AssertionError("The independent smooth-phase convergence check failed")
import numpy
import scipy
summary = {
    "status": "All numerical consistency checks passed; no formal proof certification.",
    "python": platform.python_version(),
    "numpy": numpy.__version__,
    "scipy": scipy.__version__,
    "maximum_absolute_errors": errors,
    "independent_defect_prediction": prediction,
    "independent_scaled_excess_at_n_1024": finite,
}
(RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                       encoding="utf-8")
print(json.dumps(summary, indent=2))
