"""Initial material order-rate probes of a prepared, frozen shear reference.

Use the archived flared-mouth carrier, order profiles and native order residual.
Add a declared local shear storage term, prescribing the shear component and
reference at the instant of each read. No reference dynamics or physical
relaxation constants are inferred by this forward test.
"""
import json
from pathlib import Path
from hashlib import file_digest

import numpy as np
from scipy.sparse.linalg import splu

from mouth_contact_age import PRIOR
from mouth_free_contact import Carrier
from mouth_free_order_diffuse import Order, read
from mouth_native_order_contrast import flat_window, rate_field


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def check_patch(case, operator, factor, carrier_field):
    dom = operator.dom
    order_file = Path(case['order_file'])
    with order_file.open('rb') as stream:
        if file_digest(stream, 'sha256').hexdigest() != case['order_sha256']:
            raise ValueError('Archived order mismatch')
    with np.load(order_file) as data:
        chi = data['chi']
    rc, wc = case['r'], case['w']
    q0, qr, qw = read(dom, chi, rc, wc)
    psi, pr, pw = read(dom, carrier_field, rc, wc)
    drift = (pr/psi).imag*qr + (pw/psi).imag*qw
    base = operator.measure(chi)[0]
    q = dom.N @ chi
    r, w = (dom.N @ dom.points).T
    def window(flat_r, flat_w):
        return flat_window(r-rc, flat_r, flat_r+.045)*flat_window(w-wc, flat_w, flat_w+.045)
    load = window(.11,.085)
    probe = window(.08,.065)
    anchor = .006*window(.07,.055)
    mobility = .6
    # The gradients of the elastic coefficients differ across these laws,
    # though their values at quarter order agree. The third law is a local
    # zero-derivative control for the link.
    for link_power, support_power in ((1,1), (2,3), (0,2)):
        def energy(c, strain):
            return .5*((c/q0)**link_power*strain**2 + (c/q0)**support_power*anchor**2)

        def rate(strain):
            # Differentiate the trial energy at fixed hidden displacement.
            # A complex step is exact to roundoff for the polynomials here.
            step = 1.e-20
            conjugate = energy(q+1j*step,strain).imag/step
            residual = base + operator.Nf.T @ (operator.wt * mobility * conjugate)
            euler = rate_field(dom, factor, residual)
            return float(read(dom, euler, rc, wc)[0] + drift)

        for hidden in (-.0017, 0., .0023):
            mid = rate(hidden*load)
            for pulse in (.001, .0005):
                minus = rate(hidden*load-pulse*probe)
                plus = rate(hidden*load+pulse*probe)
                second = plus+minus-2*mid
                estimate = None if abs(second)<1.e-12 else pulse*(plus-minus)/(2*second)
                emit(test="instantaneous_three_preparations", profile=case["profile"],
                     state=case["state"], mesh=case["mesh"], intake=case["intake"],
                     carrier=case["carrier_file"], order=str(order_file), q0=q0,
                     hessian_link_power=link_power, hessian_support_power=support_power,
                     prescribed_residual=hidden, probe=pulse, rates=[minus,mid,plus],
                     rate_curvature=second, inferred_residual=estimate)
                if link_power == 0 and estimate is not None:
                    raise ArithmeticError('Ungated local control returned a strain')
                if link_power and (estimate is None or abs(estimate-hidden)>3.e-5):
                    raise ArithmeticError('Signed strain readout unresolved on this mesh')


def run():
    patches = [json.loads(s) for s in json.loads(PRIOR.read_text())["runs"]
               ["mouth_mixed_order_ripple02"]["stdout_lines"]]
    patches = [p for p in patches if p["test"]=="archive_patch" and p["profile"]=="cubic"
               and p["mesh"] in (34,48)]
    for source in sorted({p['carrier_file'] for p in patches}):
        carrier = Carrier(Path(source).parent)
        model = Order(carrier.domain, carrier.field, .04)
        factor = splu(model.mass.tocsc())
        for patch in patches:
            if patch['carrier_file'] == source:
                check_patch(patch, model, factor, carrier.field)


if __name__ == '__main__':
    run()
