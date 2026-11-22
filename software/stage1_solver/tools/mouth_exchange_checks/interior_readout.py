"""Independent triangular carrier check of the smaller held-end wall contribution."""
import json
import carrier
from mouth_exchange_interior import InteriorShape
carrier.ExchangeShape=InteriorShape
for n in (16,32):
 for beta in [0.,.04/6]:
    lo,hi=[carrier.solve(e,n,beta) for e in (-.005,.005)]
    print(json.dumps(dict(n=n,beta=beta,minus=lo,plus=hi,slope=(hi-lo)/.01)),flush=True)
