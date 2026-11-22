"""Refinement of the independent kink response, including the small difference on a sphere."""
import json
import numpy as np
from order_check import energy
from scipy.special import expit
from scipy.integrate import quad
alpha=-2*quad(lambda t:t*(.5-t/3)**2/(1-t)**3,0,.5)[0]
beta=4*quad(lambda s:s*s*(expit(-s)*expit(s))**2,0,60)[0]
print(json.dumps(dict(analytic_H2=alpha,analytic_pair=beta)),flush=True)
for d in ((1,0,0),(1,-1,0),(1,1,1),(-.3,.7,.7)):
 a=np.array(d); P=(a.sum()**2-a@a)/2; theory=2*(alpha*a.sum()**2+beta*P)
 t=.0025
 Z=[energy(z,4096,d) for z in (-t,0,t)]
 E,M=np.array(Z).T
 k=(M[-1]-M[0])/(2*t); k2=(E[0]+E[-1]-2*E[1])/t**2
 print(json.dumps(dict(nodes=4097,d=d,step=t,M=M.tolist(),E=E.tolist(),k=k,k2=k2,theory=theory)),flush=True)
 assert abs(k-theory)<max(.0004*abs(theory),.0001) and abs(k2-theory)<max(.0004*abs(theory),.0001)
 assert abs(E[0]-E[-1])<1e-12 and abs(M[1])<1e-12
