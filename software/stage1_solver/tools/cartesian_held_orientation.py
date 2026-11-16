"""Held-drain force on a fixed Cartesian lattice inside a circular lateral reservoir."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.fft import dct, irfft2, next_fast_len, rfft2
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse import coo_matrix, diags


RATE = 0.25
ANGLES = {"axis": (1, 0), "diagonal": (1, 1), "oblique": (2, 1)}


class CartesianBath:
    """Disk in x,y; mirror at z=0 and bath at z=H. Depth w has two normal faces."""

    def __init__(self, cells, shift, radius=16, height=8, depth_cells=8):
        h = self.h = 1/cells
        self.dz = 1/depth_cells  # step in slab depth w
        self.volume = 2*h**3     # both signs of the out-of-plane brane coordinate z
        indices = np.arange(-radius*cells-2, radius*cells+3)
        x = (indices+.5-shift)*h
        y = (indices+.5)*h
        ix, iy = self.ix, self.iy = np.nonzero(x[:, None]**2+y**2 < radius**2)
        self.x, self.y = x[ix], y[iy]
        count = len(ix)
        ids = self.ids = np.full((len(x), len(y)), -1, dtype=int)
        ids[ix, iy] = np.arange(count)
        self.mirror = None
        if 2*shift == round(2*shift):
            mirror = ids[ix-np.rint(2*self.x/h).astype(int), iy-np.rint(2*self.y/h).astype(int)]
            if np.any(mirror<0) or not np.array_equal(mirror[mirror], np.arange(count)):
                raise ValueError("Central inversion must map the complete disk onto itself")
            self.mirror = mirror
        rows, cols = [], []
        self.wall = np.zeros(count)
        diag = np.zeros(count)
        missing = []
        for ax, ay in ((1,0),(-1,0),(0,1),(0,-1)):
            nb = ids[ix+ax, iy+ay]
            inside = nb >= 0
            rows.extend(np.flatnonzero(inside))
            cols.extend(nb[inside])
            diag[inside] += 0.5/h**2
            # Linear ghost extrapolation to psi=1 at the true circle crossing.
            # The interior edge conductance stays the same in both directions.
            outside = ~inside
            along = (ax*self.x+ay*self.y)[outside]
            across = (ay*self.x-ax*self.y)[outside]
            cut = np.sqrt(radius**2-across**2)-along
            if np.any(cut<=0) or np.any(cut>h*(1+1e-10)):
                raise ValueError("Invalid boundary crossing")
            self.wall[outside] += 0.5/(h*cut)
            missing.append((ax,ay,outside))
        off = np.full(len(rows), -0.5/h**2)
        self.plane = (coo_matrix((off,(rows,cols)),shape=(count,count)).tocsc()
                      + diags(diag+self.wall,format="csc"))
        self.z = (np.arange(height*cells)+.5)*h
        self.weights = np.ones(depth_cells+1)
        self.weights[[0,-1]] = np.sqrt(.5)
        ez = 2*(np.sin(np.pi*(np.arange(len(self.z))+.5)/(2*len(self.z)))/h)**2
        ew = 2*(np.sin(np.pi*np.arange(depth_cells+1)/(2*depth_cells))/self.dz)**2
        self.shifts = ew[:, None]+ez
        self.shape = (*self.shifts.shape, count)
        # The periodic rectangle is only a fast inverse. Replacing its boundary
        # rows removes every edge from the physical disk to the exterior.
        self.size = size = next_fast_len(len(x))
        k = 2*np.sin(np.pi*np.arange(size)/size)**2/h**2
        self.box_eigenvalues = k[:,None]+k[None,:size//2+1]
        self.bx, self.by = bx, by = ix[self.wall>0], iy[self.wall>0]
        self.extra_diagonal = self.wall[self.wall>0].copy()
        self.ghosts = []
        for ax,ay,mask in missing:
            sel = np.flatnonzero(mask[self.wall>0])
            self.extra_diagonal[sel] -= .5/h**2
            self.ghosts.append((sel,bx[sel]+ax,by[sel]+ay))
        self.inverse = {}

    def prepare_inverse(self, shift):
        """Factor only I + V G E on the circle, with G the box inverse.

        E inserts data at the disk's boundary cells. V replaces their exterior
        neighbor terms by the cut-wall diagonal in the physical operator.
        """
        table = 1/(self.shifts[:,:,None,None]+shift+self.box_eigenvalues)
        bx,by = self.bx,self.by
        def offsets(x,y):
            return ((x[:,None]-bx)%self.size,(y[:,None]-by)%self.size)
        boundary_offset = offsets(bx,by)
        ghost_offsets = [(sel,offsets(gx,gy)) for sel,gx,gy in self.ghosts]
        factors = {}
        for k,j in np.ndindex(self.shifts.shape):
            kernel = irfft2(table[k,j],s=(self.size,self.size))
            matrix = self.extra_diagonal[:,None]*kernel[boundary_offset]
            for sel,off in ghost_offsets:
                matrix[sel] += (.5/self.h**2)*kernel[off]
            matrix.flat[::len(bx)+1] += 1
            factors[k,j] = lu_factor(matrix)
        self.inverse[shift] = table,factors

    def elliptic(self, rhs, normal, shift):
        rhs[-1] += normal/self.dz
        rhs *= self.weights[:,None,None]
        modes = dct(dct(rhs, type=1, axis=0, norm="ortho", overwrite_x=True),
                    type=4, axis=1, norm="ortho", overwrite_x=True)
        if shift not in self.inverse:
            self.prepare_inverse(shift)
        table,factors = self.inverse[shift]
        box = np.zeros((*self.shifts.shape,self.size,self.size))
        def box_solve(data):
            freq = rfft2(data)
            freq *= table
            return irfft2(freq,s=(self.size,self.size))
        box[...,self.ix,self.iy] = modes
        v = box_solve(box)
        boundary = self.extra_diagonal*v[...,self.bx,self.by]
        for sel,gx,gy in self.ghosts:
            boundary[...,sel] += (.5/self.h**2)*v[...,gx,gy]
        box.fill(0.)
        for k,j in np.ndindex(self.shifts.shape):
            box[k,j,self.bx,self.by] = lu_solve(factors[k,j],boundary[k,j])
        u = (v-box_solve(box))[...,self.ix,self.iy]
        return dct(dct(u,type=4,axis=1,norm="ortho",overwrite_x=True),
                   type=1,axis=0,norm="ortho",overwrite_x=True)/self.weights[:,None,None]

    def footprint(self, center):
        distance2 = (self.x-center[0])**2+(self.y-center[1])**2+self.z[:,None]**2
        return np.maximum(1-distance2, 0)**3

    def verify(self, f, normal, currents):
        # Independent real-space equation, including the cut wall, endcaps and mouth face.
        u = f-1
        neq = (self.plane @ u.reshape(-1,len(self.x)).T).T.reshape(f.shape)
        end = np.concatenate((u[:,:1],u,-u[:,-1:]),axis=1)
        neq -= np.diff(end,n=2,axis=1)/(2*self.h**2)
        top = u[-2]+2*self.dz*normal
        normal_ghosts = np.concatenate((u[1:2],u,top[None]),axis=0)
        neq -= np.diff(normal_ghosts,n=2,axis=0)/(2*self.dz**2)
        neq += (abs(f)**8-1)*f/4
        error = float(np.max(abs(neq)))
        exits = (f.imag @ self.wall).sum(axis=1)+f[:,-1].imag.sum(axis=1)/self.h**2
        bath = 2*self.volume*self.dz*np.dot(self.weights**2,exits)
        rates = self.volume*currents.sum(axis=(1,2))
        if (not np.isfinite([error,bath,*rates]).all() or error>1e-8 or
                np.max(abs(rates-RATE))>5e-14 or abs(bath-rates.sum())>2e-8):
            raise ArithmeticError(f"Unverified field: {error=} {bath=} {rates=}")
        return dict(equation_residual=error,bath_rate=float(bath),mouth_rates=rates.tolist())

    def solve(self, footprints, calibration=None):
        shapes = footprints if calibration is None else footprints/calibration
        restoring = footprints.sum(axis=0)
        def boundary(face):
            n = abs(face)**2
            raw = shapes if calibration is None else shapes*n
            totals = self.volume*raw.sum(axis=(1,2))
            currents = raw*(RATE/totals)[:,None,None]
            normal = (-restoring*(1-1/abs(face))+1j*currents.sum(axis=0)/n)*face
            return currents,normal
        f = np.ones(self.shape,complex)
        for step in range(1,501):
            _, g = boundary(f[-1])
            p = (abs(f)**8-1)/4
            re = self.elliptic(2*(f.real-1)-p*f.real,g.real,2)
            im = self.elliptic(-p*f.imag,g.imag,0)
            updated = 1+re+1j*im
            delta = np.max(abs(updated-f))
            if not np.isfinite(delta):
                raise ArithmeticError("Nonfinite nonlinear iteration")
            f = updated
            if delta<5e-13:
                currents,g = boundary(f[-1])
                stats = self.verify(f,g,currents)
                return f[-4:].copy(), currents, dict(steps=step,**stats)
        raise ArithmeticError(f"Steady pair failed to converge: {delta}")

    def loads(self, layers, current, weight):
        """Physical x,y components. z vanishes by orbital-plane reflection."""
        iz, ip = np.nonzero(weight>0)
        cut = (iz,ip)
        face = layers[-1]
        n = abs(face)**2
        slope = -2*weight[cut]*(1-1/np.sqrt(n[cut]))
        def normal(a):
            return (11*a[-1]-18*a[-2]+9*a[-3]-2*a[-4])/(6*self.dz)
        fw, nw = normal(layers), normal(abs(layers)**2)
        carried, stress, direct = [], [], []
        for ax,ay in ((1,0),(0,1)):
            idx = [self.ids[self.ix[ip]+ax*k,self.iy[ip]+ay*k] for k in (-2,-1,1,2)]
            if min(i.min() for i in idx)<0:
                raise ValueError("Force stencil reaches the bath")
            def grad(a):
                l2,l1,r1,r2 = (a[iz,i] for i in idx)
                return (l2-8*l1+8*r1-r2)/(12*self.h)
            df = grad(face)
            carried.append(self.volume*np.sum(current[cut]*(face[cut].conj()*df).imag/n[cut]))
            stress.append(self.volume*np.sum(slope*grad(n)/4))
            direct.append(self.volume*np.sum((df.conj()*fw[cut]).real-grad(nw)/4))
        return np.array([carried,stress,direct])


def measure(slab, distance, angle):
    u = np.array(ANGLES[angle],float)
    u /= np.linalg.norm(u)
    perp = np.array([-u[1],u[0]])
    target, source = (slab.footprint(sign*distance*u/2) for sign in (-1,1))
    alone, j0, acheck = slab.solve(target[None])
    baseline = slab.loads(alone,j0[0],target)
    n0 = abs(alone[-1])**2
    # A generic fractional translation breaks central inversion on this grid.
    # In that case calibrate the companion independently at its actual position.
    extra = {}
    if slab.mirror is None:
        source_alone, js, scheck = slab.solve(source[None])
        ns = abs(source_alone[-1])**2
        source_baseline = slab.loads(source_alone,js[0],source)
        extra["source_alone"] = scheck
    else:
        if np.max(abs(source-target[:,slab.mirror]))>3e-14:
            raise ArithmeticError("Companion and target discretizations differ under inversion")
        ns = n0[:,slab.mirror]
        source_baseline = -baseline
    pair, jp, pcheck = slab.solve(np.array([target,source]),np.array([n0,ns]))
    target_load = slab.loads(pair,jp[0],target)-baseline
    other_load = slab.loads(pair,jp[1],source)-source_baseline
    symmetry_error = np.max(abs(target_load+other_load))
    if not np.isfinite([target_load,other_load]).all():
        raise ArithmeticError("Nonfinite load")
    if slab.mirror is not None and symmetry_error>3e-12:
        raise ArithmeticError(f"Paired antisymmetry failed: {symmetry_error}")
    if slab.mirror is None:
        extra["source_force"] = other_load[:2].sum(axis=0).tolist()
        extra["source_direct"] = other_load[2].tolist()
        extra["sum_of_mouth_loads_max_abs"] = float(symmetry_error)
    else:
        extra["paired_antisymmetry_error"] = float(symmetry_error)
    carried,stress,direct = target_load
    total = carried+stress
    return dict(d=distance,angle=angle,direction=u.tolist(),alone=acheck,pair=pcheck,
                **extra,isolated_load=baseline.tolist(),
                carried=carried.tolist(),stress=stress.tolist(),force=total.tolist(),
                parallel=float(total@u),perpendicular=float(total@perp),
                direct_parallel=float(direct@u),direct_perpendicular=float(direct@perp))


def compare(paths):
    rows = {}
    for p in paths:
        for line in p.read_text().splitlines():
            obj = json.loads(line)
            key = obj["cells"],obj["shift"],obj["d"],obj["angle"]
            if key in rows:
                raise ValueError("Duplicate readout")
            rows[key] = obj
    keys = {(c,s,d,a) for c in (2,3,4) for s in (0,.5) for d in (8,16) for a in ANGLES}
    if set(rows) != keys:
        raise ValueError("Missing Cartesian comparison")
    for d in (8,16):
        for c in (2,3,4):
            group = [rows[c,s,d,a] for s in (0,.5) for a in ANGLES]
            for p, t in (("parallel","perpendicular"),("direct_parallel","direct_perpendicular")):
                radial = np.array([r[p] for r in group])
                relative = np.array([abs(r[t]/r[p]) for r in group])
                direction = max(abs(rows[c,s,d,a][p]/rows[c,s,d,"axis"][p]-1)
                                for s in (0,.5) for a in ANGLES)
                translation = max(abs(rows[c,.5,d,a][p]/rows[c,0,d,a][p]-1) for a in ANGLES)
                worst = group[np.argmax(relative)]
                print(f"d={d} cells={c} readout={p} min_radial={min(radial):.12g} "
                      f"max_radial={max(radial):.12g} angular_diff_percent={100*direction:.9g} "
                      f"half_cell_diff_percent={100*translation:.9g} "
                      f"max_abs_transverse_over_radial_percent={100*max(relative):.9g} "
                      f"worst_angle={worst['angle']} worst_shift={worst['shift']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    subs = ap.add_subparsers(dest="action",required=True)
    test = subs.add_parser("measure")
    test.add_argument("--cells",type=int,choices=(2,3,4),required=True)
    test.add_argument("--shift",type=float,choices=(0,.5),required=True)
    test.add_argument("--checkpoint-dir",type=Path,required=True)
    summary = subs.add_parser("compare")
    summary.add_argument("readouts",nargs="+",type=Path)
    args = ap.parse_args()
    if args.action=="measure":
        args.checkpoint_dir.mkdir(parents=True,exist_ok=False)
        slab = CartesianBath(args.cells,args.shift)
        for d in (8,16):
            for name in ANGLES:
                row = dict(cells=args.cells,shift=args.shift,**measure(slab,d,name))
                line = json.dumps(row,allow_nan=False)
                (args.checkpoint_dir/f"d{d}_{name}.json").write_text(line+"\n")
                print(line,flush=True)
    else:
        compare(args.readouts)
