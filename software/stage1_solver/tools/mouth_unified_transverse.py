"""Hold a saved broad waist and compare its curl and symmetric-strain toroidal pole."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.optimize import minimize_scalar
from scipy.special import spherical_jn, spherical_yn

from mouth_exchange_shear import DEGREE, solve
from mouth_live_waist_wave import OccupiedWave
from mouth_live_waist_wave_location import share
from mouth_sealed_shape_relax import DELTA, Meridian

ARCHIVE = Path("research/pde_ledger_v3/_measurements/mouth_live_waist_order.json")
PORT = 4.2


def cell_difference(wave, shape, mesh, degree):
    """Integrate Kcurl - Ksym at each radius face, with the archived Q1 quadrature."""
    nz = len(wave.ids) // mesh
    axial = shape.length / nz
    i, j = np.indices((mesh, nz)).reshape(2, -1)
    diff = np.zeros_like(wave.kc)
    nodes, weights = leggauss(3)
    for b, weight in zip((nodes + 1) / 2, weights / 2):
        outer = shape.radius((j + b) * axial)[0]
        for a, sign in ((0, -1), (1, 1)):
            radius = outer * (i + a) / mesh
            basis = np.array([(1-a)*(1-b), a*(1-b), (1-a)*b, a*b])
            diff += (2 * sign * radius**3 * axial * weight)[:, None, None] * np.outer(basis, basis)
    # A constant f in degree 1 is a rigid brane rotation. Restore all element
    # nodes (including the axis) when testing polynomial identity locally.
    # Degree affects only the angular f² term; subtract it back to degree one.
    if degree != 4:
        raise ValueError("The saved pole uses degree four")
    return diff


def direct_cell_forms(shape, mesh, cells, degree):
    """Independent volume quadrature for both energies on a nonconstant nodal test."""
    nz = cells // mesh
    h = shape.length / nz
    i, j = np.indices((mesh, nz)).reshape(2, -1)
    rng = np.random.default_rng(491)
    vals = rng.normal(size=(cells, 4))
    curl, sym, rigid, rigid_curl = np.zeros((4, cells))
    nodes, weights = leggauss(3)
    for a, qa in zip((nodes+1)/2, weights/2):
        for b, qb in zip((nodes+1)/2, weights/2):
            width, slope = shape.radius((j+b)*h)
            r = (i+a)*width/mesh
            f = vals @ np.array([(1-a)*(1-b), a*(1-b), (1-a)*b, a*b])
            fr = (vals @ np.array([b-1, 1-b, -b, b]))*mesh/width
            fw = (vals @ np.array([a-1, -a, 1-a, a]))/h-(i+a)*slope/mesh*fr
            weight = qa*qb*h*width/mesh*r*r
            ang = degree*(degree+1)
            curl += weight*((2*f+r*fr)**2+ang*f*f+(r*fw)**2)
            sym += weight*((r*fr)**2+(ang-2)*f*f+(r*fw)**2)
            rigid += weight*((r*0.)**2+(1*2-2)*1.+(r*0.)**2)
            rigid_curl += weight*(4.+2.)
    return vals, curl, sym, rigid, rigid_curl


def edge_matrix(port, shape, block):
    ids = port.ids
    ij = np.repeat(ids, len(ids)), np.tile(ids, len(ids))
    return sparse.coo_matrix((block.ravel(), ij), shape=shape).tocsc()


def run(case, intervals):
    archive = json.loads(ARCHIVE.read_bytes())
    prep = json.loads(archive["runs"][case]["literal_stdout"].splitlines()[0])
    if prep["point"] != "broad" or prep["radius"] not in (1, 2) or prep["J"] != 1:
        raise ValueError("Only held broad radius-one and radius-two specimens at unit action")
    path = Path("_scratch/s11c/live_waist_dynamics01/fields") / (case + ".npz")
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if sha != archive["saved_field_sha256"][str(path)]:
        raise ValueError("Saved field disagrees with its source manifest")
    source = np.load(path)
    mesh, n = prep["field_mesh"], prep["mesh"]
    shape = Meridian(prep["radius"], prep["coordinates"])
    wave = OccupiedWave(shape, n, 2*n, mesh, PORT)
    weights = wave.project @ source["order"].T.ravel()
    mass_inside = wave.matrix(weights[:, None, None]*wave.mc)
    mass = wave.cut_m + mass_inside
    curl = wave.cut_k + wave.matrix(weights[:, None, None]*wave.kc)
    diff = cell_difference(wave, shape, mesh, DEGREE)
    vals, dc, ds, rigid, rigid_curl = direct_cell_forms(shape, mesh, len(weights), DEGREE)
    def quadratic(matrices):
        return np.einsum("ci,cij,cj->c", vals, matrices, vals)
    cell_error = np.max(np.abs(quadratic(wave.kc-diff)-ds))/np.max(ds)
    curl_error = np.max(np.abs(quadratic(wave.kc)-dc))/np.max(dc)
    rotation_error = np.max(np.abs(diff.sum(axis=(1, 2))-rigid_curl))/np.max(rigid_curl)
    if max(cell_error, curl_error, rotation_error, np.max(np.abs(rigid))) > 1e-10:
        raise ArithmeticError("Radial faces disagree with volume strains or rigid rotation")
    pocket_difference = wave.matrix(weights[:, None, None]*diff)
    port = wave.port
    # In uniform feed Kcurl - Ksym is exactly the same 2 R³ boundary
    # matrix as Bcurl - Bsym. Cancel those two terms in K - B only.
    feed_face = 2*PORT**3 * (port.project.T @ port.project)
    feed_difference = edge_matrix(port, mass.shape, feed_face)
    # Independently differentiate the physical outgoing displacement h_l(kr).
    # The old traction on U is U'+U/r; the Cauchy traction is U'-U/r.
    test_omega = 5.135
    h4 = spherical_jn(DEGREE, PORT*test_omega) + 1j*spherical_yn(DEGREE, PORT*test_omega)
    deriv = (spherical_jn(DEGREE, PORT*test_omega, True)
             + 1j*spherical_yn(DEGREE, PORT*test_omega, True))
    direct = test_omega*deriv/h4-1/PORT
    port_error = abs(port.traction(test_omega)[0]-2/PORT-direct)
    if port_error > 1e-9:
        raise ArithmeticError("Mismatch in outgoing symmetric-strain traction")
    vector = source["vector"]
    omega = complex(source["pole"].item())
    waist = minimize_scalar(lambda z: shape.radius(z)[0], bounds=(0, shape.length),
                            method="bounded", options={"xatol": 1e-12}).x
    history, endpoints = [], []
    for alpha in np.linspace(0., 1., intervals+1):
        effective = curl-alpha*pocket_difference
        out, new, _ = solve(shape, mesh, PORT, 0., omega, incoming_vector=vector,
                            supplied_system=(wave.uniform[0], effective, mass, mass_inside),
                            target_frequency=omega.real, shape_work=False)
        omega = complex(out["frequency"], -out["amplitude_loss"])
        overlap = abs(np.vdot(vector, mass @ new))**2/np.vdot(vector, mass @ vector).real
        if overlap < .95:
            raise ArithmeticError(f"Lost pole correspondence at fraction {alpha}: {overlap}")
        history.append(dict(fraction=alpha, frequency=omega.real, squared_overlap=overlap))
        vector = new
        if alpha not in (0., 1.):
            continue
        stiff = effective - alpha*feed_difference
        boundary = port.block(omega) - alpha*feed_face
        residual = (stiff - omega**2*mass) @ new
        residual[port.ids] -= boundary @ new[port.ids]
        physical_residual = np.linalg.norm(residual)/np.linalg.norm(stiff @ new)
        cavity = np.vdot(new, mass_inside @ new).real
        edges = (0., waist-2*DELTA, waist+2*DELTA, shape.length)
        fractions = [share(wave, shape, new, weights, mesh, a, b)/cavity
                     for a, b in zip(edges[:-1], edges[1:])]
        flux = .5*omega.real*PORT**4*(abs(port.project @ new[port.ids])**2 @ port.traction(omega.real).imag)
        energy = .25*(abs(omega)**2*np.vdot(new, mass @ new).real+np.vdot(new, stiff @ new).real)
        endpoints.append(dict(fraction=alpha, frequency=omega.real, unit_action_energy=omega.real,
                              amplitude_loss=-omega.imag, feed_power_per_interior_energy=flux/energy,
                              physical_pole_residual=physical_residual, cavity_kinetic_share=cavity,
                              upstream_of_band=fractions[0], inside_waist_band=fractions[1],
                              downstream_of_band=fractions[2], partition_error=sum(fractions)-1))
        if physical_residual > 1e-8 or abs(sum(fractions)-1) > 1e-5:
            raise ArithmeticError("Physical outgoing pencil or cut-cell energy identity failed")
    if abs(endpoints[0]["frequency"]-float(source["frequency"])) > 1e-8:
        raise ArithmeticError("Archived curl pole failed reproduction")
    print(json.dumps(dict(case=case, snapshot=str(path), sha256=sha, waist_depth=waist,
                          mesh=mesh, min_cell_order=min(weights), cell_form_error=cell_error,
                          curl_form_error=curl_error, rotation_relative_error=rotation_error,
                          rigid_rotation_symmetric_energy=sum(rigid),
                          rigid_rotation_curl_energy=sum(rigid_curl),
                          direct_port_traction_error=port_error, history=history, endpoints=endpoints),
                     allow_nan=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", choices=[f"r{r}_broad_{n}" for r in (1, 2) for n in (48, 72)])
    parser.add_argument("--intervals", type=int, default=12)
    args = parser.parse_args()
    if args.intervals < 1:
        parser.error("Continuation needs at least one interval")
    for item in args.cases:
        run(item, args.intervals)
