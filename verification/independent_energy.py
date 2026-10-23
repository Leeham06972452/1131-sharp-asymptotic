"""Independent check of the proposed global continuum defect; outputs stay beside this script."""
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_jacobi, roots_legendre, gamma

ROOT = Path(__file__).resolve().parent
CSTAR = float(2*gamma(1.25)**2/gamma(.75)**2)

def objective(y):
    n=len(y)
    differences=y[:,None]-y[None,:]
    np.fill_diagonal(differences,1.)
    logs=-np.log(np.abs(differences)).sum(axis=1)
    weights=np.prod(np.sign(differences),axis=1)*np.exp(logs-logs.max())
    x,q=roots_legendre(n)
    terms=weights[None,:]/(x[:,None]-y[None,:])
    cardinals=terms/terms.sum(axis=1)[:,None]
    return float(q@(cardinals**2).sum(axis=1))

def reference(n,alpha=1.25):
    interior,_=roots_jacobi(n-2,alpha,alpha)
    return np.r_[1.,interior[::-1],-1.]

def defect_linear_modulus(amplitude,order):
    x,w=roots_jacobi(order,.25,.25)
    difference=x[:,None]-x[None,:]
    # (cosh(2*a*d)-1)/d**2 = 2*a**2*sinhc(a*d)**2.
    value=amplitude*difference
    sinhc=np.ones_like(value)
    np.divide(np.sinh(value),value,out=sinhc,where=value!=0)
    return float(w@(2*amplitude**2*sinhc**2)@w/(2*np.pi))

predictions=[{'quadrature_order':k,'defect':defect_linear_modulus(.2,k)} for k in (128,256,512)]
rows=[]
for n in (64,128,256,512,1024):
    x=reference(n)
    theta=np.arccos(x)
    y=np.cos(theta+.2*np.sin(theta)/n)
    y[0],y[-1]=1.,-1.
    ix,iy=objective(x),objective(y)
    rows.append({'n':n,'reference_coefficient':n*(2-ix),
                 'perturbed_coefficient':n*(2-iy),
                 'scaled_excess':n*(iy-ix)})
result={'candidate_constant':CSTAR,'prediction_modulus_eta':'0.2*x',
        'phase_delta':'0.2*sin(theta)','continuum_defect':predictions,
        'finite_objectives':rows,'status':'Independent floating point checks, not a proof.'}
(ROOT/'results'/'independent_energy.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
