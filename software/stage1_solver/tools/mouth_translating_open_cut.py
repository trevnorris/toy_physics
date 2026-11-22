"""Test the rigid carrier's weak translation source when its remote traces are freed.

A radial test lives only at large r and vanishes on all solid faces with
nonzero x normal. Integration by parts separates the artificial-cut term
from physical advection, which is measured independently in the volume.
"""
import argparse
from hashlib import file_digest
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_free_center import basis, motion_source
from mouth_carrier_lateral_response import emit, prepare, trace
from mouth_occupied_push_baths import load_bath


def volume_advection(dom, field, test):
    root, wt = leggauss(4)
    u, v = (s.ravel() for s in np.meshgrid(root, root, indexing='ij'))
    weights = np.outer(wt, wt).ravel()
    total = 0j
    for block, (pg, vg, dofs) in dom.blocks.items():
        for start in range(0, len(dofs), 128):
            i, j = np.unravel_index(np.arange(start, min(start+128, len(dofs))),
                                    (len(pg)-1, len(vg)-1))
            i, j = i[:, None], j[:, None]
            hp, hv = np.diff(pg)[i]/2, np.diff(vg)[j]/2
            p, s = pg[i] + hp*(1+u), vg[j] + hv*(1+v)
            idx, r, _, det, n, gr, _ = basis(dom, block, p, s, i, j)
            g = np.sum(n*test[idx], axis=-1)
            fr = np.sum(gr*field[idx], axis=-1)
            total += np.sum(4*np.pi*r*r*hp*hv*det*weights*g*fr)
    return -1j*total


def port(dom, field, test, block, along, fixed, radial_normal):
    pg, vg, _ = dom.blocks[block]
    edges = pg if along == 'p' else vg
    root, wt = leggauss(8)
    half = np.diff(edges)[:, None]/2
    q = (edges[:-1, None] + half*(1+root)).ravel()
    other = np.full_like(q, fixed)
    p, v = (q, other) if along == 'p' else (other, q)
    r, _, rp, rv, wp, wv, val, fr, _ = trace(dom, block, p, v, np.column_stack((field, test)))
    tr, tw = (rp, wp) if along == 'p' else (rv, wv)
    ds = np.hypot(tr, tw)
    measure = 4*np.pi*r*r*(half*wt).ravel()*ds
    # Core source conventions have divided through the 1/3 dipole angular average.
    weak = .5j*radial_normal*np.sum(measure*val[:, 0]*val[:, 1])
    mean_square = np.sum(measure*abs(fr[:, 0])**2)/measure.sum()
    return weak, float(np.sqrt(mean_square)), float(np.max(abs(fr[:, 0])))


def pair(z):
    return [float(z.real), float(z.imag)]


def check(dom, field, omega, label):
    radial = dom.points[:, 0]
    test = np.maximum(2*radial/dom.far-1, 0)**2
    raw = motion_source(dom, field)
    real, imag = np.split(raw, 2)
    old = test @ (real + 1j*imag)
    removed = 0j
    ports = {}
    for name, args in [('feed', (1, 'v', dom.far, 1)),
                       ('bulk_side', (4, 'p', 1., 1)),
                       ('roof_inner', (3, 'v', dom.roof, 0)),
                       ('roof_outer', (4, 'v', dom.roof, 0))]:
        weak, rms, maximum = port(dom, field, test, *args)
        removed += weak
        # RMS refers to the yz-radial coefficient multiplying the brane angular
        # dipole Y=x/r, before its 1/sqrt(3) angular RMS factor.
        ports[name] = dict(source_face=pair(weak), base_fr_rms=rms, base_fr_max=maximum,
                           lab_chart_offset_rms_per_unit_speed=rms/omega)
    volume = volume_advection(dom, field, test)
    emit(test='freed_cut_source', background=label, far=dom.far, roof=dom.roof, omega=omega,
         old_source_pairing=pair(old), artificial_cut_pairing=pair(removed),
         corrected_source_pairing=pair(old-removed), measured_volume_advection=pair(volume),
         absolute_residual=float(abs(old-removed-volume)), faces=ports)
    np.testing.assert_allclose(old-removed, volume, rtol=2e-8, atol=3e-9)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    prep = parser.add_mutually_exclusive_group(required=True)
    prep.add_argument('--source', type=Path)
    prep.add_argument('--bath', type=Path)
    parser.add_argument('--omega', type=float, required=True)
    args = parser.parse_args()
    if args.omega <= 0:
        parser.error('Positive omega required.')
    dom, field = load_bath(args.bath) if args.bath else prepare(args)
    path = (args.bath or args.source)/'field.npz'
    with path.open('rb') as stream:
        digest = file_digest(stream, 'sha256').hexdigest()
    emit(test='source_field', path=str(path), sha256=digest, nodes=len(field))
    check(dom, np.ones_like(field), args.omega, 'uniform zero intake')
    check(dom, field, args.omega, 'archived drain')


if __name__ == '__main__':
    main()
