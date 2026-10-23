#!/usr/bin/env python3
"""Independent finite checks for the regularity/scattering continuation.

The accompanying proofs, not these floating point checks, establish the
claims. Requires Python 3, NumPy and SciPy. No network or input files.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import legvander
from scipy.integrate import quad
from scipy.linalg import solve, eigvalsh
from scipy.special import roots_jacobi, roots_legendre


def nodes(n, alpha=1.25):
    interior, _ = roots_jacobi(n-2, alpha, alpha)
    return np.r_[1., interior[::-1], -1.]


def basis(x, n):
    return legvander(x, n-1)*np.sqrt((2*np.arange(n)+1)/2)


def nodal_values(x, y):
    """Chebyshev-leading Q and Q' away from nodes, using signed products."""
    x = np.asarray(x)
    distances = x[:, None]-y[None, :]
    logs = (len(y)-1)*math.log(2)+np.sum(np.log(np.abs(distances)), axis=1)
    q = np.prod(np.sign(distances), axis=1)*np.exp(logs)
    qp = q*np.sum(1/distances, axis=1)
    return q, qp


def nodal_derivatives(y):
    d = y[:, None]-y[None, :]
    np.fill_diagonal(d, 1.)
    return np.prod(np.sign(d), axis=1)*np.exp(
        (len(y)-1)*math.log(2)+np.sum(np.log(np.abs(d)), axis=1))


def finite_checks(n, case):
    x = nodes(n)
    theta = np.arccos(x)
    if case == "baseline":
        y = x.copy()
    elif case == "smooth_phase":
        y = np.cos(theta+0.23*np.sin(2*theta)/n)
        y[0], y[-1] = 1., -1.
    else:
        raise ValueError(case)
    v = basis(y,n)
    vi = solve(v,np.eye(n))
    mass = vi.T@vi
    lam = 1/np.sum(v*v,axis=1)
    gram = np.sqrt(lam[:,None]*lam[None,:])*(v@v.T)
    comparison = np.sqrt(lam[:,None]*lam[None,:])*gram
    excess = float(np.trace(mass)-sum(lam))
    defect = mass+comparison-2*np.diag(lam)
    theoretical = (np.sqrt(lam[:,None]*lam[None,:]) *
                   (gram+solve(gram,np.eye(n))-2*np.eye(n)))
    mass_identity_error = float(np.max(np.abs(defect-theoretical)))
    mass_inequality_min = float(eigvalsh(2*np.diag(lam)+excess*np.eye(n)-mass)[0])

    # Rank-one identity used in the new discrete A2 proof.
    dp = nodal_derivatives(y)
    keep = np.flatnonzero(np.abs(y)<0.8)
    yy, dd = y[keep], dp[keep]
    epsilon = 0.7/n
    matrix = 1j*epsilon*dd[:,None]/dd[None,:]/(yy[:,None]-yy[None,:]+1j*epsilon)
    commutator = (yy[:,None]-yy[None,:])*matrix+1j*epsilon*matrix
    rank_one = 1j*epsilon*dd[:,None]/dd[None,:]
    rank_one_error = float(np.max(np.abs(commutator-rank_one)))
    a2 = 1.
    for start in range(len(dd)):
        a = np.cumsum(dd[start:]**2)
        b = np.cumsum(dd[start:]**-2)
        lengths = np.arange(1,len(dd)-start+1)
        a2 = max(a2,float(np.max(a*b/lengths**2)))

    # Exact exterior logarithmic identity. The stable integrand is a
    # sum of log1p ratios; common endpoint factors cancel identically.
    ax, ay = 1-x[1:-1]**2, 1-y[1:-1]**2
    def exterior_integrand(t):
        return float(np.sum(np.log1p((ay-ax)/(t*t+ax))))
    exterior, estimate = quad(exterior_integrand,0,np.inf,epsabs=2e-11,
                              epsrel=2e-11,limit=300)
    transport = math.pi*float(np.sum(np.sqrt(ay)-np.sqrt(ax)))

    # Exact reciprocal Prüfer mass, including nodal endpoints.
    raw, weights = roots_legendre(max(1024,12*n))
    angles = (raw+1)*math.pi/2
    weights = weights*math.pi/2
    q, qp = nodal_values(np.cos(angles),y)
    s = np.sin(angles)
    envelope = q*q+(s*qp/n)**2
    reciprocal_integral = float(weights@(s*s/envelope))
    reciprocal_discrete = math.pi*n*float(np.sum(dp**-2))

    # The exact derivative-weighted interpolation identity used to
    # select the weighted inverse without a Fredholm uniqueness claim.
    testx, _ = roots_legendre(n+7)
    q, qp = nodal_values(testx,y)
    h = lambda z: (1-z*z)*(z+0.17)
    data = dp*h(y)/n
    interpolated = basis(testx,n)@vi@data
    divided = (h(testx)[:,None]-h(y)[None,:])/(testx[:,None]-y[None,:])
    exact = h(testx)*qp/n-q*np.sum(divided,axis=1)/n
    interpolation_error = float(np.max(np.abs(interpolated-exact)))

    assert mass_identity_error < 2e-10
    assert mass_inequality_min > -2e-10
    assert rank_one_error < 2e-12
    assert abs(exterior-transport) < 2e-9
    assert abs(reciprocal_integral-reciprocal_discrete) < 2e-7
    assert interpolation_error < 2e-9
    return {
        "n": n, "case": case, "objective":float(np.trace(mass)),
        "n_christoffel_excess":n*excess,
        "n_mass_operator_norm":n*float(eigvalsh(mass)[-1]),
        "mass_identity_max_error":mass_identity_error,
        "mass_bound_min_eigenvalue":mass_inequality_min,
        "rank_one_commutator_max_error":rank_one_error,
        "maximum_bulk_discrete_A2":a2,
        "transport_direct":transport,"transport_exterior":exterior,
        "transport_error":exterior-transport,
        "quadrature_estimated_transport_error":estimate,
        "reciprocal_prufer_integral":reciprocal_integral,
        "reciprocal_prufer_discrete":reciprocal_discrete,
        "reciprocal_prufer_error":reciprocal_integral-reciprocal_discrete,
        "derivative_weighted_interpolation_max_error":interpolation_error,
    }


def amplitude_check(n,alpha):
    y=nodes(n,alpha)
    raw,weights=roots_legendre(max(1024,12*n))
    angles=0.25+(raw+1)*(math.pi-0.5)/2
    weights=weights*(math.pi-0.5)/2
    q,qp=nodal_values(np.cos(angles),y)
    stable=np.exp(1j*n*angles)*(q-1j*np.sin(angles)*qp/n)
    target=(1-np.exp(2j*angles))**(1.5-alpha)
    return {"n":n,"alpha":alpha,
            "bulk_stable_amplitude_L2_error":float(np.sqrt(weights@abs(stable-target)**2)),
            "bulk_scattering_L2_error":float(np.sqrt(weights@abs(
                stable/np.conj(stable)-target/np.conj(target))**2))}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="erdos_1131_regularity_numerics.json")
    args=parser.parse_args()
    result={
        "status":"Floating point checks of exact identities and limit predictions; not proof certificates.",
        "finite_identities":[finite_checks(n,c) for n in (16,32,64,128)
                             for c in ("baseline","smooth_phase")],
        "amplitude_limits":[amplitude_check(n,a) for a in (.3,1.25,2.2)
                            for n in (32,64,128,256)],
    }
    Path(args.output).write_text(json.dumps(result,indent=2)+"\n")
    checks=result["finite_identities"]
    print(json.dumps({"cases":len(checks),
        "max_transport_error":max(abs(r["transport_error"]) for r in checks),
        "max_reciprocal_error":max(abs(r["reciprocal_prufer_error"]) for r in checks),
        "max_interpolation_error":max(r["derivative_weighted_interpolation_max_error"] for r in checks),
        "output":str(Path(args.output).resolve())},indent=2))


if __name__=="__main__":
    main()
