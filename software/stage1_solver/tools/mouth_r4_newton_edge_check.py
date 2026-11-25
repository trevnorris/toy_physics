"""Use exact triple products of eight stored labels; no field or trial evaluations."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json
import hashlib
import numpy as np

base = Path('_scratch/s11c/r4_newton_held01')
src = base/'primary/last_active.npz'
assert hashlib.sha256(src.read_bytes()).hexdigest() == 'a2210361aa5b38712de7db618b5220985e12f39f3a089351dfd3767bef0a7178'
record = json.loads(Path('research/pde_ledger_v3/_measurements/mouth_r4_newton_held_stop.json').read_bytes())['limiting_cell']
cell = tuple(record['cell'])
bits = list(product((0, 1), repeat=3))
with np.load(src, allow_pickle=False) as f:
    labels = {b: tuple(map(Q, f['labels'][tuple(i+j for i,j in zip(cell,b))])) for b in bits}

def value(point):
    cols = [[Q(0) for j in range(3)] for i in range(3)]
    for b in bits:
        w = [p if v else 1-p for p,v in zip(point,b)]
        for d in range(3):
            c = (1 if b[d] else -1)*w[(d+1)%3]*w[(d+2)%3]
            for j in range(3): cols[d][j] += c*labels[b][j]
    a,b,c=cols
    return sum(a[i]*(b[(i+1)%3]*c[(i-1)%3]-c[(i+1)%3]*b[(i-1)%3]) for i in range(3))

norm=Q(record['accepted_chart_det_exact'])/value(tuple(map(Q,record['local'])))
det=lambda p: norm * value(p)
corner=[float(det(v)) for v in bits]
p0=det((0,0,0)); p1=det((0,0,1)); pm=det((0,0,Q(1,2)))
a=2*(p1-2*pm+p0); b=p1-p0-a; t= -b/(2*a)
assert 0<t<1 and a>0
m=det((0,0,t))
assert m==p0-b*b/(4*a) and m>0
print(json.dumps(dict(vertex_order=bits, vertex_values=corner, min_vertex=min(corner), minimizing_fraction_on_blocking_edge=float(t), min_on_blocking_edge=float(m), witness=float(det(tuple(map(Q,record['local'])))), min_vertex_over_edge_min=float(min(corner)/float(m)), normalizing_cell_volume=float(1/norm))))
