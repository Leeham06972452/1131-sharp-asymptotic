#!/usr/bin/env python3
"""Finite numerical checks for the global Erdős 1131 continuation.

The accompanying proofs establish the theorems. These floating point checks
test identities, normalizations, and deliberately ill-conditioned examples.
Requires Python 3, numpy, scipy. No network access or input files are needed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legvander
from scipy.linalg import eigvalsh, solve, svdvals
from scipy.optimize import brentq
from scipy.special import eval_jacobi, gamma, roots_jacobi, roots_legendre


ALPHA = 1.25
CSTAR = float(2 * gamma(1.25)**2 / gamma(0.75)**2)


def basis(nodes):
    n = len(nodes)
    return legvander(nodes, n - 1) * np.sqrt((2*np.arange(n)+1)/2)


def baseline(n, alpha=ALPHA):
    interior, weights = roots_jacobi(n-2, alpha, alpha)
    x = np.r_[1.0, interior[::-1], -1.0]
    weights = weights[::-1]
    inverse = solve(basis(x), np.eye(n))
    mass = inverse.T @ inverse
    omega = weights * (1-x[1:-1]**2)**(-alpha)
    return x, inverse, mass, omega


def barycentric_objective(nodes):
    n = len(nodes)
    distances = nodes[:, None] - nodes[None, :]
    np.fill_diagonal(distances, 1.0)
    signs = np.prod(np.sign(distances), axis=1)
    logs = -np.sum(np.log(np.abs(distances)), axis=1)
    bary = signs * np.exp(logs-np.max(logs))
    t, weights = roots_legendre(n)
    values = np.empty((n,n))
    for i, ti in enumerate(t):
        diff = ti-nodes
        exact = np.flatnonzero(np.abs(diff)<3e-15)
        if len(exact):
            values[i] = 0
            values[i, exact[0]] = 1
        else:
            raw = bary/diff
            values[i] = raw/np.sum(raw)
    return float(weights @ np.sum(values**2, axis=1))


def test_shift(n, name, x, inverse, mass, delta):
    theta = np.arccos(x)
    y = np.cos(theta+delta/n)
    matrix = basis(y) @ inverse
    z = solve(matrix, np.eye(n))
    new_coefficients = inverse @ z
    old_objective = float(np.trace(mass))
    new_objective = float(np.sum(new_coefficients**2))
    w0 = np.pi*np.sin(theta)
    # J T^-1 J=(A* - A^-1)(A - A^-*).
    h = matrix.T-z
    positive_energy = float(w0 @ np.sum(h*h, axis=1))
    direct_new_inverse = solve(basis(y), np.eye(n))
    direct_objective = float(np.sum(direct_new_inverse**2))
    singular = svdvals(matrix)
    record = {
        "n": n, "case": name,
        "changed_nodes": int(np.count_nonzero(delta)),
        "total_squared_shift": float(delta @ delta),
        "mean_squared_shift": float(delta @ delta/n),
        "largest_absolute_shift": float(np.max(np.abs(delta))),
        "minimum_scaled_angular_gap": float(n*np.min(np.diff(theta+delta/n))),
        "smallest_singular_value": float(singular[-1]),
        "inverse_operator_norm": float(1/singular[-1]),
        "baseline_objective": old_objective,
        "new_objective": new_objective,
        "scaled_objective_cost": float(n*(new_objective-old_objective)),
        "positive_column_gram_energy": positive_energy,
        "cost_minus_energy": float(n*(new_objective-old_objective)-positive_energy),
        "independent_direct_inverse_error": direct_objective-new_objective,
    }
    if n in (40,160) and name == "nearly_colliding_pair":
        check = barycentric_objective(y)
        record["independent_barycentric_error"] = check-new_objective
        assert abs(check-new_objective) < 2e-9
    assert abs(direct_objective-new_objective) < 2e-9
    assert positive_energy >= 0
    return record


def exterior_factorization_check(n, outer_order=480, inner_order=None):
    """Independent positive-integral evaluation of z*Rz.

    c=1/u reduces the exterior interval to (0,1). Jacobi quadrature
    handles u^(1-2r)(1-u)^r exactly as the endpoint weight.
    """
    x, inverse, mass, omega = baseline(n)
    m = n-2
    xi = x[1:-1]
    a = 1-xi*xi
    derivative = 0.5*(m+2*ALPHA+1)*eval_jacobi(m-1, ALPHA+1, ALPHA+1, xi)
    v = 1/(a*derivative)
    z = (-1.0)**np.arange(m) * np.sqrt(a)
    z /= np.linalg.norm(z)
    remainder = n*(np.diag(omega)-mass[1:-1,1:-1])
    matrix_value = float(z @ remainder @ z)
    r = 2-ALPHA
    t, wt = roots_jacobi(inner_order or (6*n+160), ALPHA, ALPHA)
    wp2 = wt * eval_jacobi(m, ALPHA, ALPHA, t)**2
    raw_u, wu = roots_jacobi(outer_order, r, 1-2*r)
    u = (raw_u+1)/2
    wu = wu / 2**(2-r)
    jtilde = np.sum(wp2[None,:]/(1-u[:,None]*t[None,:]), axis=1)
    positive_mode = np.sum((z*v)[None,:]/(1-u[:,None]*xi[None,:]),axis=1)
    negative_mode = np.sum((z*v)[None,:]/(1+u[:,None]*xi[None,:]),axis=1)
    integral = n*math.sin(math.pi*r)/math.pi * np.sum(
        wu * (1+u)**r * jtilde * (positive_mode**2+negative_mode**2))
    assert abs(integral-matrix_value) < 2e-8
    return {"n":n, "outer_order":outer_order, "matrix_quadratic_form":matrix_value,
            "positive_exterior_integral":float(integral),
            "difference":float(integral-matrix_value)}


def smooth_phase_prediction(amplitude, order=220):
    """Picone coefficient for f=g*exp(2a*T2), independently quadratured.

    This is a diagnostic prediction for the smooth phase experiment, not
    an asserted theorem about the perturbed arrays in this script.
    """
    x, w = roots_jacobi(order, 0.25, 0.25)
    diff = x[:,None]-x[None,:]
    ratio = 2*amplitude*(x[:,None]+x[None,:])
    dh = ratio*diff
    # (cosh(2dh)-1)/diff² = 2*ratio²*(sinh(dh)/dh)².
    sinhc = np.ones_like(dh)
    mask = np.abs(dh)>1e-10
    sinhc[mask] = np.sinh(dh[mask])/dh[mask]
    integrand = 2*ratio**2*sinhc**2
    deficit = float(w @ integrand @ w/(2*np.pi))
    return {"amplitude":amplitude, "predicted_deficit":deficit,
            "predicted_coefficient":CSTAR-deficit,
            "status":"diagnostic; smooth-array asymptotic not proved by this computation"}


def global_row_checks(n, x, inverse, mass, omega):
    """Probe the new uniform row estimates in bulk and shrinking edge layers.

    Independent Legendre evaluation constructs the cardinals. The grid
    includes distances in angle from 10^-3/n to 60/n at both endpoints.
    These diagnostics test the scaling, not the analytic supremum bound.
    """
    edge = np.geomspace(1e-3, 60.0, 400)/n
    theta = np.unique(np.r_[np.linspace(0, np.pi, 1201),
                            edge, np.pi-edge, np.arccos(x)])
    theta = theta[(theta >= 0) & (theta <= np.pi)]
    evaluation = legvander(np.cos(theta),n-1)*np.sqrt((2*np.arange(n)+1)/2)
    values = evaluation @ inverse
    w0 = np.pi*np.sqrt(np.maximum(0, 1-x*x))
    square_sum = np.sum(values*values, axis=1)
    weighted_error = n*(values*values @ w0-np.pi*np.sin(theta))
    remainder = n*(np.diag(omega)-mass[1:-1,1:-1])
    interior = values[:,1:-1]
    rows = n*np.sum((interior @ remainder)*interior, axis=1)
    assert float(np.min(rows)) > -2e-8
    j = int(np.argmax(rows))
    k = int(np.argmax(np.abs(weighted_error)))
    legendre = basis(x)
    christoffel_sum = float(np.sum(1/np.sum(legendre*legendre,axis=1)))
    return {"n":n, "grid_size":len(theta),
            "maximum_cardinal_square_sum":float(np.max(square_sum)),
            "maximum_n_mass_deficit_row":float(rows[j]),
            "angle_at_maximum_mass_deficit_row":float(theta[j]),
            "maximum_absolute_n_weighted_square_sum_error":float(abs(weighted_error[k])),
            "angle_at_maximum_weighted_square_sum_error":float(theta[k]),
            "maximum_n_gaussian_weight_error":float(np.max(
                n*np.abs(n*omega-w0[1:-1]))),
            "baseline_n_christoffel_excess":float(
                n*(np.trace(mass)-christoffel_sum))}


def loop_data(nodes, z=2.0):
    """Independent Cauchy integrals and exact-loop residual at a small array."""
    n = len(nodes)
    inverse = solve(basis(nodes), np.eye(n))
    mass = inverse.T @ inverse
    t, weights = roots_legendre(max(80, 2*n+20))
    polynomials = legvander(t, n-1) * np.sqrt((2*np.arange(n)+1)/2)
    cardinals = polynomials @ inverse
    squared_sum = np.sum(cardinals**2, axis=1)
    f = float(weights @ (squared_sum/(z-t)))
    c = float(np.diag(mass) @ (1/(z-nodes)))
    cp = float(-np.diag(mass) @ (1/(z-nodes)**2))
    q = np.prod(t[:,None]-nodes[None,:], axis=1)
    distances = nodes[:,None]-nodes[None,:]
    np.fill_diagonal(distances,1.0)
    reciprocal_derivatives = 1/np.prod(distances,axis=1)
    hp = float(-weights @ (q*q/(z-t)**2))
    r2 = float(np.sum(reciprocal_derivatives**2/(z-nodes)**2))
    log_derivative = float(np.sum(1/(z-nodes)))
    objective = float(np.trace(mass))
    residual = cp-objective/(1-z*z)-2*log_derivative*(f-c)+hp*r2
    return {"n":n, "nodes":nodes.tolist(), "z":z, "objective":objective,
            "F_minus_C":f-c, "stationary_loop_residual":residual}


def stationary_loop_checks():
    # Reflection symmetry reduces the four-node stationarity equations
    # to one derivative. Complex-step differentiation avoids cancellation.
    def derivative(a):
        ac = a+1e-25j
        x = np.array([-1.0,-ac,ac,1.0],dtype=complex)
        inverse = solve(basis(x),np.eye(4))
        return float(np.imag(np.sum(inverse*inverse))/1e-25)
    stationary_a = brentq(derivative,0.3,0.6,xtol=5e-15)
    records = [loop_data(np.array([-1.0,0.0,1.0])),
               loop_data(np.array([-1.0,-stationary_a,stationary_a,1.0]))]
    for record in records:
        assert abs(record["stationary_loop_residual"]) < 2e-12
    records[-1]["one_parameter_derivative"] = derivative(stationary_a)
    return records


def jacobi_loop_diagnostics():
    # The fixed Jacobi arrays are not asserted to be finite-n stationary.
    # This records their first-order loop, whose limit vanishes at alpha=5/4.
    records=[]
    for alpha in (1.0,1.15,1.25,1.4):
        for n in (40,80,160,320):
            x,inverse,mass,_ = baseline(n,alpha)
            t,weights=roots_legendre(2*n+30)
            values=(legvander(t,n-1)*np.sqrt((2*np.arange(n)+1)/2))@inverse
            f=float(weights @ (np.sum(values*values,axis=1)/(2-t)))
            c=float(np.diag(mass) @ (1/(2-x)))
            records.append({"alpha":alpha,"n":n,"z":2.0,
                            "F_minus_C":f-c,"n_F_minus_C":n*(f-c),
                            "n_squared_F_minus_C":n*n*(f-c)})
    return records


def run():
    mass_records, shifts, coherent, row_records = [], [], [], []
    for n in (40,80,160,320):
        x, inv, mass, omega = baseline(n)
        theta = np.arccos(x)
        rem = n*(np.diag(omega)-mass[1:-1,1:-1])
        whitened = rem/np.sqrt(n*omega[:,None]*n*omega[None,:])
        eig = eigvalsh(whitened)
        mass_records.append({
            "n":n, "baseline_coefficient":float(n*(2-np.trace(mass))),
            "quadrature_weight_sum":float(np.sum(omega)),
            "mass_deficit_trace":float(np.trace(rem)),
            "smallest_normalized_deficit_eigenvalue":float(eig[0]),
            "largest_normalized_deficit_eigenvalue":float(eig[-1]),
            "minimum_relative_mass_eigenvalue":float(1-eig[-1]),
        })
        row_records.append(global_row_checks(n,x,inv,mass,omega))
        bulk = (theta >= np.pi/6) & (theta <= 5*np.pi/6)
        alternating = np.zeros(n)
        alternating[bulk] = 0.8/np.sqrt(n)*(-1.0)**np.flatnonzero(bulk)
        shifts.append(test_shift(n,"spread_alternating_weak_shifts",x,inv,mass,alternating))

        smooth = np.zeros(n)
        phase = (theta[bulk]-np.pi/6)/(2*np.pi/3)
        smooth[bulk] = 0.6*n**(-0.25)*np.sin(np.pi*phase)**2
        shifts.append(test_shift(n,"spread_coherent_vanishing_shifts",x,inv,mass,smooth))

        collision = np.zeros(n)
        j = n//2-1
        half_gap = n*(theta[j+1]-theta[j])/2
        epsilon = n**(-0.25)
        collision[j] = half_gap-epsilon
        collision[j+1] = -(half_gap-epsilon)
        shifts.append(test_shift(n,"nearly_colliding_pair",x,inv,mass,collision))

        amp = 0.20
        fixed_phase = amp*np.sin(2*theta)
        fixed_phase[[0,-1]]=0
        rec = test_shift(n,"fixed_smooth_phase_outside_vanishing_mean_class",x,inv,mass,fixed_phase)
        rec["new_coefficient"] = float(n*(2-rec["new_objective"]))
        coherent.append(rec)

    factorization = [exterior_factorization_check(n) for n in (22,42,82)]
    return {"date":"2026-10-23", "candidate_C_star":CSTAR,
            "warning":"Numerical checks are not proofs or a solution of the unrestricted problem.",
            "baseline_mass_checks":mass_records,
            "global_endpoint_row_checks":row_records,
            "global_bulk_shift_checks":shifts,
            "exterior_factorization_checks":factorization,
            "smooth_phase_experiment":coherent,
            "smooth_phase_prediction":smooth_phase_prediction(0.20),
            "stationary_loop_checks":stationary_loop_checks(),
            "jacobi_loop_diagnostics":jacobi_loop_diagnostics()}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("erdos_1131_global_numerics.json"))
    args=parser.parse_args()
    results=run()
    args.output.write_text(json.dumps(results,indent=2)+"\n")
    print(json.dumps({"saved":str(args.output),"candidate_C_star":CSTAR,
        "factorization_errors":[r["difference"] for r in results["exterior_factorization_checks"]],
        "collision_cases":[{k:r[k] for k in ("n","inverse_operator_norm","scaled_objective_cost")}
            for r in results["global_bulk_shift_checks"] if r["case"]=="nearly_colliding_pair"],
        "smooth_phase_coefficient_prediction":results["smooth_phase_prediction"]["predicted_coefficient"]},indent=2))
