"""Signed elastic work of an equilibrated flat puncture and two attached annular specimens."""
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_sealed_shape_relax import Meridian

ROOT = Path(__file__).resolve().parents[3]
RECORD = ROOT / "research/pde_ledger_v3/_measurements/mouth_live_waist_order.json"
INSET = .5  # A chosen specimen, the outer radial half of the archived ordered duct.
OFFSET = .001  # The same signed length in all preparations.
REST_AREA = .25**2  # Sets the amplitude of the volume-preserving reference field.
SHEAR = BULK = 1.  # Pure ordered material; both moduli carry one power of frozen order.


def flat_puncture(radius):
    """Infinite lateral extent, H=1, two free thickness faces, prescribed spherical rim."""
    # At the current aperture an outward shift q from the held state gives
    # u_r=(q-d)(R/r)^2 relative to a natural aperture offset d. In linear
    # elasticity the extension to the entire slab is an exact equilibrium.
    x, weight = leggauss(32)
    x, weight = (1+x)/2, weight/2  # r=R/x
    s, depth_weight = leggauss(12)
    s, depth_weight = s/2, depth_weight/2
    # A depth-dependent test has the same integrated shear modulus as the
    # unit material. It respects both free faces because sigma_sw=sigma_ww=0.
    profiles = {"uniform": np.ones_like(s), "graded": 1. + .6*s}
    r = radius/x
    measure = 4*np.pi*r**2*radius/x**2*weight
    v_hoop = (radius/r)**2/r
    v_radial = -2*v_hoop
    result = []
    for name, mu in profiles.items():
        for d in (-OFFSET, 0., OFFSET):
            def energy(q):
                return ((depth_weight@mu)*(measure @ (v_radial**2+2*v_hoop**2))*(q-d)**2)
            h = .0002
            central = (energy(-h)-energy(h))/(2*h)
            aperture_stress_rr = 4*mu*d/radius
            surface = 4*np.pi*radius**2*(depth_weight @ aperture_stress_rr)
            # Interior radial equilibrium, d sigma_rr/dr+2(sigma_rr-sigma_h)/r.
            radial_div = -12*d*radius**2/r**4 + 2*(4+2)*d*radius**2/r**4
            row = dict(profile=name, natural_radius_offset=d, applied_radius=radius,
                       elastic_opening_force=float(central), from_rim_traction=float(surface),
                       from_formula=float(16*np.pi*radius*d*(depth_weight@mu)),
                       face_traction=0., max_radial_divergence=float(max(abs(radial_div))),
                       energy_center=float(energy(0)))
            if abs(central-surface) > 1.e-12:
                raise ArithmeticError("The through-hole energy and rim traction disagree")
            result.append(row)
    return result


def archived(radius, mesh, manifest):
    name = f"r{radius}_broad_{mesh}"
    data = [json.loads(line) for line in manifest["runs"][name]["literal_stdout"].splitlines()]
    prep = next(r for r in data if r["kind"] == "preparation")
    center = next(r for r in data if r["kind"] == "coupled" and r["action"] == 1)
    path = Path("_scratch/s11c/live_waist_dynamics01/fields") / (name+".npz")
    sha = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    if sha != manifest["saved_field_sha256"][str(path)] or prep["radius"] != radius:
        raise ValueError("Specimen differs from registered preparation")
    with np.load(ROOT/path) as source:
        chi = source["order"]
        if abs(source["frequency"]-center["wave_frequency"]) > 1.e-10:
            raise ValueError("Occupation differs from archived state")
    if chi.shape != (2*mesh+1, mesh+1) or prep["mesh"] != mesh or prep["point"] != "broad":
        raise ValueError("Unexpected order chart")
    return Meridian(radius, np.array(prep["coordinates"])), chi, dict(
        snapshot=str(path), sha256=sha, archived_phase_bias=prep["phase_bias"],
        archived_action=1., archived_inward_current=center["current"],
        archived_waist=center["waist_radius"], coordinates=prep["coordinates"])


class Annulus:
    """Q1 meridional r,w displacements, two identical angular normal strains on S²."""

    def __init__(self, shape, order):
        self.shape = shape
        self.mesh = order.shape[1]-1
        t = np.linspace(0, 1, 2*self.mesh+1)
        x = np.arange(self.mesh//2, self.mesh+1)/self.mesh
        self.x, self.t = np.meshgrid(x, t)
        self.r = self.x*shape.radius(self.t*shape.length)[0]
        self.chi = order[:, self.mesh//2:]
        ids = np.arange(self.r.size).reshape(self.r.shape)
        nodes = np.stack((ids[:-1, :-1], ids[:-1, 1:], ids[1:, :-1], ids[1:, 1:]), -1)
        self.nodes = nodes.reshape(-1, 4)
        self.dofs = np.stack((2*self.nodes, 2*self.nodes+1), -1).reshape(-1, 8)
        self.free = np.repeat(((self.x > INSET) & (self.x < 1) &
                               (self.t > 0) & (self.t < 1)).ravel(), 2)

    def integration(self):
        q, weights = leggauss(3)
        x0 = self.x[:-1, :-1].ravel()
        t0 = self.t[:-1, :-1].ravel()
        dx, dt = 1/self.mesh, 1/(2*self.mesh)
        length = self.shape.length
        for a, wa in zip((1+q)/2, weights/2):
            x = x0+a*dx
            for b, wb in zip((1+q)/2, weights/2):
                t = t0+b*dt
                width, slope = self.shape.radius(t*length)
                r = width*x
                basis = np.array([(1-a)*(1-b), a*(1-b), (1-a)*b, a*b])
                dr = np.array([b-1, 1-b, -b, b])/(width[:, None]*dx)
                dw = np.array([a-1, -a, 1-a, a])/(length*dt)-(x*slope)[:, None]*dr
                hoop = basis[None, :]/r[:, None]
                vol = 4*np.pi*r**2*width*length*dx*dt*wa*wb
                chi = self.chi.ravel()[self.nodes] @ basis
                yield dr, dw, hoop, vol, chi

    def stiffness(self):
        local = np.zeros((len(self.nodes), 8, 8))
        volume = 0.
        # In four dimensions lambda=K-mu/2. Physical shear strain is
        # (du_r/dw+du_w/dr)/2; the fourth row below is engineering shear.
        weights = np.array([2*SHEAR, 2*SHEAR, 4*SHEAR, SHEAR])
        for dr, dw, hoop, vol, chi in self.integration():
            form = np.zeros((len(vol), 4, 8))
            form[:, 0, ::2], form[:, 1, 1::2] = dr, dw
            form[:, 2, ::2] = hoop
            form[:, 3, ::2], form[:, 3, 1::2] = dw, dr
            trace = form[:, 0] + form[:, 1] + 2*form[:, 2]
            gram = np.einsum("cik,i,cij->ckj", form, weights, form)
            local += (vol*chi)[:, None, None]*(
                gram+(BULK-SHEAR/2)*trace[:, :, None]*trace[:, None, :])
            volume += sum(vol)
        rows = np.broadcast_to(self.dofs[:, :, None], local.shape)
        cols = np.broadcast_to(self.dofs[:, None, :], local.shape)
        count = 2*self.r.size
        return sparse.coo_matrix((local.ravel(), (rows.ravel(), cols.ravel())),
                                 shape=(count, count)).tocsc(), volume

    def boundary_modes(self):
        # Mode 0 is the tangent of rho^3=r^3+3 d A f(t), w_ref=w.
        # This globally compatible reference map on the annulus has det=1
        # exactly and the identity value and gradient at the two joins.
        f = 16*(self.t*(1-self.t))**2
        fields = np.zeros((*self.r.shape, 2, 3))
        fields[:, :, 0, 0] = REST_AREA*f/self.r**2
        fields[:, :, 0, 1] = self.x*4*self.t*(1-self.t)  # saved widening per q
        fields[:, :, 1, 2] = self.t  # saved lengthening per actual L
        if np.min(self.r**3-3*OFFSET*REST_AREA*f) <= 0:
            raise ValueError("Reference map does not preserve a positive bore")
        return fields.reshape(-1, 3)

    def stress_work(self, relative, virtual):
        # Contract the physical stress on all four spatial directions against
        # the original test extension, independently of the reaction rows.
        u = relative.reshape(-1, 2)[self.nodes]
        v = virtual.reshape(-1, 2, 2)[self.nodes]
        total = np.zeros(2)
        for dr, dw, hoop, vol, chi in self.integration():
            def gradient(g, values):
                return np.einsum("ci,ci...->c...", g, values)
            rr = gradient(dr, u[:, :, 0])
            ww = gradient(dw, u[:, :, 1])
            hh = gradient(hoop, u[:, :, 0])
            gamma = gradient(dw, u[:, :, 0])+gradient(dr, u[:, :, 1])
            pressure = (BULK-SHEAR/2)*(rr+ww+2*hh)
            wrk = ((2*SHEAR*rr+pressure)[:, None]*gradient(dr, v[:, :, 0]) +
                   (2*SHEAR*ww+pressure)[:, None]*gradient(dw, v[:, :, 1]) +
                   2*(2*SHEAR*hh+pressure)[:, None]*gradient(hoop, v[:, :, 0]) +
                   SHEAR*gamma[:, None]*(gradient(dw, v[:, :, 0])+gradient(dr, v[:, :, 1])))
            total += np.sum((vol*chi)[:, None]*wrk, axis=0)
        return total


def curved(radius, mesh, manifest):
    shape, order, source = archived(radius, mesh, manifest)
    solid = Annulus(shape, order)
    matrix, volume = solid.stiffness()
    trial = solid.boundary_modes()
    # The unknowns are displacements relative to the natural placement. For
    # each boundary trace the interior is free, so KU=0 at all interior DOFs.
    response = trial.copy()
    free = solid.free
    solver = splu(matrix[free][:, free])
    response[free] = solver.solve(-(matrix[free][:, ~free] @ trial[~free]))
    reaction = matrix @ response
    gram = response.T @ reaction
    rows = []
    for d in (-OFFSET, 0., OFFSET):
        displacement = -d*response[:, 0]
        elastic = -(trial[~free, 1:].T @ (-d*reaction[~free, 0]))
        volume_force = -solid.stress_work(displacement, trial[:, 1:])
        fd = []
        for m in (1, 2):
            def energy(step):
                moved = displacement+step*response[:, m]
                return .5*moved @ (matrix @ moved)
            h = .0002
            fd.append((energy(-h)-energy(h))/(2*h))
        row = dict(natural_offset=d, generalized_elastic_force=elastic.tolist(),
                   force_from_volume_stress=volume_force.tolist(), from_signed_energy=fd,
                   energy=float(.5*displacement@(matrix@displacement)))
        if max(np.max(abs(elastic-fd)), np.max(abs(elastic-volume_force))) > 1.e-10:
            raise ArithmeticError("Force from boundary reactions, stress and signed energy differs")
        rows.append(row)
    return dict(radius=radius, grid=mesh, volume=float(volume), **source,
                physical_depth=shape.length, radial_fraction=[INSET, 1.],
                reference_area=REST_AREA, stress_readings=rows,
                max_interior_residual=float(np.max(abs(reaction[free]))),
                stiffness_symmetry_error=float(np.max(abs(gram-gram.T))),
                trace_relaxed_stiffness=gram.tolist())


def main():
    manifest = json.loads(RECORD.read_bytes())
    flat = {str(radius): flat_puncture(radius) for radius in (1, 2)}
    annuli = [curved(radius, mesh, manifest) for radius in (1, 2) for mesh in (48, 72)]
    output = dict(ordered_shear=SHEAR, ordered_bulk=BULK, flat=flat, annuli=annuli,
                  convention="Positive material force opens or lengthens. Linear elasticity.",
                  curved_specimen="b/2 <= r <= b. Both radial faces and both ends registered. "
                                  "Order frozen from the archived occupied state, mu=K=chi. "
                                  "Only crystal virtual work; the physical load-bearing region is unknown.")
    print(json.dumps(output, allow_nan=False, indent=2))
    if any(r["max_interior_residual"] > 1.e-11 for r in annuli):
        raise ArithmeticError("Interior elastic equilibrium did not settle")


if __name__ == "__main__":
    main()
