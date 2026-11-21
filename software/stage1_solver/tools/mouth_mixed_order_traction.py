"""Compare the boundary force in two independently assembled material chains."""
import json

import numpy as np
from scipy.sparse import bmat, eye, lil_matrix
from scipy.sparse.linalg import spsolve


def material_readings(q, omega, cells=180):
    # Equal local rho, mu and relaxation time in unit conventions. Clamp the
    # right end and impose a unit transverse displacement on the left. Read
    # just the material force in the first bond, with a massless end fixture.
    h = 3/cells
    gradient = lil_matrix((cells, cells-1))
    for i in range(cells-1):
        gradient[i, i] = 1/h
        gradient[i+1, i] = -1/h
    gradient = gradient.tocsr()
    end_strain = np.zeros(cells)
    end_strain[0] = -1/h
    n = cells-1
    ident = eye(n, format="csr")
    spring = (gradient.T @ gradient).tocsr()

    drag = q*(1-q)
    coupled = bmat([[q*spring-(omega**2*q+1j*omega*drag)*ident, 1j*omega*drag*ident],
                    [1j*omega*drag*ident, (-omega**2*(1-q)-1j*omega*drag)*ident]], format="csc")
    rhs = np.r_[-q*(gradient.T @ end_strain), np.zeros(n)]
    split = spsolve(coupled, rhs)[:n]
    split_stress = q*(gradient @ split+end_strain)

    # Ordinary spring in parallel with a positive spring and dashpot in
    # series. The internal variable is the dashpot strain on each bond.
    # Fade the extra branch below the tested order interval to ensure that
    # even instantaneous stiffness vanishes at zero order.
    t = min(4*q, 1.)
    activity = t*t*(3-2*t)
    extra = (1-q)*activity
    instantaneous = q+extra
    memory = bmat([[instantaneous*spring-omega**2*ident, -extra*gradient.T],
                   [-gradient, (1-1j*omega)*eye(cells)]], format="csc")
    rhs = np.r_[-instantaneous*(gradient.T @ end_strain), end_strain]
    one = spsolve(memory, rhs)
    displacement = one[:n]
    strain = gradient @ displacement + end_strain
    one_stress = q*strain + extra*(strain-one[n:])
    ratio = one_stress[0]/split_stress[0]
    spatial_difference = float(max(abs(displacement-split)))
    # A fixed-displacement comparison is valid with returning waves present.
    predicted_ratio = 1/(q+(1-q)/(1-1j*omega/q))
    print(json.dumps(dict(chi=q, omega=omega, cells=cells,
        spring_at_rest=q, extra_spring=extra, max_displacement_difference=spatial_difference,
        stress_split=[split_stress[0].real, split_stress[0].imag],
        stress_memory=[one_stress[0].real, one_stress[0].imag],
        stress_ratio_abs=abs(ratio), stress_ratio_phase_deg=np.angle(ratio, deg=True),
        unmodified_ratio_abs=abs(predicted_ratio), unmodified_ratio_phase_deg=np.angle(predicted_ratio, deg=True))),
        flush=True)
    if q >= .25:
        assert spatial_difference < 1e-9 and abs(ratio-predicted_ratio) < 1e-9
    return spatial_difference


if __name__ == "__main__":
    for chi in (.25, .75, .1):
        for x in (.1, 1., 10.):
            material_readings(chi, x*chi)
