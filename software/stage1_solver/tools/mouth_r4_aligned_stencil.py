"""Frozen R4 determinant means over chart corners or the aligned 39-point rule.

The tensor Simpson/Lobatto weight at each of the 12 edge midpoints is split
equally between that midpoint and a moving point on the same edge. Other
weights retain their tensor-rule values. For positive endpoint values A,C and
middle Bernstein coefficient B, evaluating the edge quadratic at
t=sqrt(A)/(sqrt(A)+sqrt(C)) gives a positive multiple of B+sqrt(A*C).
Its sign tests the entire edge. Faces and cell interiors still need the
independent whole-cell certificate. Only determinant sampling changes here.
"""
from itertools import product
import numpy as np


VERTICES = np.array(list(product((0, 1), repeat=3)))
GRID = np.array(list(product((0, 1, 2), repeat=3)))
MIDCOUNT = np.sum(GRID == 1, axis=1)
WEIGHTS = np.array([1, 2, 16, 64])[MIDCOUNT]  # common denominator 216
EDGE_WEIGHT = 2
EDGES = [(m-stride, m, m+stride) for m in np.flatnonzero(MIDCOUNT == 1)
         for stride in [np.array([9, 3, 1])[GRID[m] == 1].item()]]


def gradient_basis(local):
    factors = np.where(VERTICES, local, 1-np.asarray(local))
    return np.column_stack([(2*VERTICES[:,d]-1)*factors[:,(d+1)%3]*factors[:,(d+2)%3]
                            for d in range(3)])


BASIS = tuple(gradient_basis(q/2) for q in GRID)


def samples(points, widths):
    """Return (27 Lobatto determinants, 12 edge t, 12 moving determinants).

    points[c,v,component] are the eight transported nodal labels in 000,
    001,...,111 order. widths[c,axis] are the corresponding chart increments.
    Shape callers can pass complex arrays for holomorphic shape derivatives.
    """
    dtype = np.result_type(points, widths, np.longdouble)
    points = np.asarray(points, dtype=dtype)
    widths = np.asarray(widths, dtype=dtype)
    if points.ndim != 3 or points.shape[1:] != (8,3) or widths.shape != (len(points),3):
        raise ValueError('Expected one eight-node chart and three widths per cell')
    if np.any(widths.real <= 0) or not np.all(np.isfinite(widths)):
        raise ValueError('Nonpositive or unresolved chart width')
    offsets = (points-points[:,:1]).transpose(0,2,1)
    det = []
    for basis in BASIS:
        matrix = (offsets @ basis)/widths[:,None,:]
        a,b,c = matrix.transpose(2,0,1)
        det.append(np.sum(a*np.cross(b,c), axis=1))
    det = np.stack(det,axis=1)
    # Positive endpoints define the moving rule even for an inverted midpoint,
    # permitting diagnostics on rejected charts and off-center counterexamples.
    left,mid,right = np.array(EDGES).T
    A, M, C = det[:,left], det[:,mid], det[:,right]
    if np.any(A.real <= 0) or np.any(C.real <= 0) or not np.isfinite(det).all():
        raise ValueError('Nonpositive chart corner or nonfinite determinant')
    sa,sc = np.sqrt(A),np.sqrt(C)
    t,u = sa/(sa+sc),sc/(sa+sc)
    moved = u*(1-2*t)*A+4*t*u*M+t*(2*t-1)*C
    return det,t,moved


def barrier_means(points, widths):
    """Return (original 8-corner mean, aligned mean), one value per cell."""
    det, _, moved = samples(points, widths)
    if np.any(det.real <= 0) or np.any(moved.real <= 0) or not np.isfinite(moved).all():
        raise ValueError('Nonpositive or unresolved chart determinant on aligned stencil')
    def density(j):
        return j-1-np.log(j)
    fixed = density(det)
    old = np.mean(fixed[:,MIDCOUNT == 0],axis=1)
    aligned = (fixed@WEIGHTS+EDGE_WEIGHT*np.sum(density(moved),axis=1))/216
    return old,aligned
