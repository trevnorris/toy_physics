"""Finite-strain R4 crystal with a mobile core–collar bond and angular material flow.

The independent P1(triangle) × P1(angle) offsets are radial fraction, depth
fraction and polar angle in the carrier's six-coordinate wall chart. At the
origin polar angle is a direction, and its limiting value is extended from
the first ring. The center has one axial displacement. Polar reflection
sets meridional tangential motion to zero at the pole and equator. Axisymmetric
azimuthal twist is exactly minimized by zero: replacing phi by phi+omega
leaves J and chi unchanged and adds chi*(R*sin(theta))²*|grad_ref omega|²/2.

The integration follows every material and label triangle. Radial integration
slices split at all crossings of deformed background order and conversion
elements. Depth and angle use composite Gauss integration, with additional
breaks on both undeformed meshes. The same finite energy, including the motion
of quadrature points at these crossings, generates the solid and order forces.
"""
import ctypes as ct
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import os
import subprocess
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator

from mouth_internal_collar_check import read
from mouth_joined_angular_work import bump
from mouth_supercore_shape import preparation


class InvalidPlacement(ValueError):
    """A trial leaves the physical crystal chart or reverses its orientation."""


class UnresolvedActivation(ArithmeticError):
    """Free order can activate at a quadrature site without an oriented label."""

    def __init__(self, info):
        self.info = info
        super().__init__(f"Undefined order deposition at {info['nonpositive_J_points']} quadrature "
                         f"sites, minimum J={info['min_nonpositive_J']} at "
                         f"{info['min_nonpositive_J_site']}.")


Pointer = ct.POINTER(ct.c_double)


class Grid(ct.Structure):
    _fields_ = [('n', ct.c_int), ('v', Pointer)]


class Data(ct.Structure):
    _fields_ = ([(s, Grid) for s in (
        'x t a depth_edges angle_edges base_x base_t field_x field_t field_a cut_x cut_t cut_a reg_e reg_t gauss').split()]
        + [(s, Pointer) for s in 'weights labels base correction offset spline shape'.split()]
        + [('spline_n', ct.c_int), ('diagnose', ct.c_int)]
        + [(s, ct.c_double) for s in 'length ref_length pressure core_ratio epsilon'.split()])


class Results(ct.Structure):
    _fields_ = [(s, Pointer) for s in ('gradient', 'load', 'stats', 'melted', 'undefined_nodes')]


def pointer(values):
    return values.ctypes.data_as(Pointer)


@lru_cache(maxsize=1)
def kernel():
    source = Path(__file__).resolve().with_suffix('.cpp')
    folder = source.parents[3]/'_scratch'/'s11c'/'joint_r4_solid_build'
    folder.mkdir(parents=True, exist_ok=True)
    library = folder/(sha256(source.read_bytes()).hexdigest()+'.so')
    if not library.exists():
        fd, stage = tempfile.mkstemp(suffix='.so', dir=folder)
        os.close(fd)
        try:
            subprocess.run(['g++', '-O3', '-std=c++17', '-fPIC', '-shared',
                            str(source), '-o', stage], check=True)
            os.replace(stage, library)
        finally:
            Path(stage).unlink(missing_ok=True)
    function = ct.CDLL(str(library)).solid_eval
    function.argtypes = [ct.POINTER(Data), ct.POINTER(Results)]
    function.restype = ct.c_int
    return function


def partition(arrays, subdivisions):
    a = np.sort(np.concatenate(arrays))
    a = a[np.r_[True, np.diff(a) > 2e-13]]
    return np.r_[(a[:-1, None]+np.diff(a)[:, None] *
                 np.arange(subdivisions)[None, :]/subdivisions).ravel(), a[-1]]


class Solid:
    def __init__(self, state, shape, epsilon, n, k, n_angle, model, *, degree=2, subdivisions=1, integration_axes=None):
        if min(n, k, n_angle, degree) < 2 or subdivisions < 1:
            raise ValueError('At least two cells/points per direction are required')
        chart, labels, inner, _ = read(state)
        if len(chart.co) != 1 or shape.neck != 4:
            raise ValueError('Requires the registered radius-four preparation')
        self.sh = self.shape = shape
        self.epsilon = float(epsilon)
        self.model = model
        self.melted_labels = False
        self.cut_node = 2*n
        self.cut = np.exp(-.5*np.exp(chart.co[0]))
        self.x = self.radial_mesh(n)
        self.t = np.linspace(0, 1, k+1)
        self.a = np.linspace(0, np.pi/2, n_angle+1)
        self.chart_shape = (k+1, len(self.x), n_angle+1, 3)
        self.ids = np.asarray(model.ids[1, 0])
        self.base_order = inner[2]
        ref, _, registration = preparation(4)
        d = self.data = Data()
        self.buffers = {}
        def use(name, arr, grid=False):
            data = self.buffers[name] = np.ascontiguousarray(arr, dtype=np.float64)
            setattr(d, name, Grid(len(data), pointer(data)) if grid else pointer(data))
        for key, arr in dict(x=self.x,t=self.t,a=self.a,base_t=inner[0],base_x=inner[1],
                             field_t=model.z[1]/shape.length,field_x=model.rr[0],field_a=model.a,
                             reg_e=chart.a,reg_t=chart.t).items():
            use(key, arr, True)
        # A caller comparing field representations can pin a single physical
        # partition containing both meshes. These planes are independent of
        # the field interpolation and shared by both crystal variations.
        physical_axes=(model.rr[0],model.a,model.z[1]) if integration_axes is None else integration_axes
        for name,axis,need in zip(('cut_x','cut_a','cut_t'),physical_axes,
                                  (model.rr[0],model.a,model.z[1])):
            axis=np.asarray(axis)
            if np.any(np.diff(axis)<=0) or max(np.min(abs(axis[:,None]-need),axis=0))>2e-12:
                raise ValueError('Solid integration partition must contain every field element boundary')
            use(name,axis/shape.length if name=='cut_t' else axis,True)
        use('depth_edges', partition([self.t, inner[0], self.buffers['cut_t'],
                                      chart.t, np.array([.06,.62])], subdivisions), True)
        use('angle_edges', partition([self.a, self.buffers['cut_a']], subdivisions), True)
        uz, wz = leggauss(degree)
        use('gauss', (uz+1)/2, True)
        use('weights', wz/2)
        use('base', inner[2])
        use('labels', labels)
        # Rescale scipy's PPoly from its native argument to normalized depth.
        s = shape.reference
        d.spline_n = len(s.x)-1
        scaled_knots = s.x/shape.unit_length
        if np.max(abs(scaled_knots-np.linspace(0,1,len(s.x))))>1e-12:
            raise ValueError('Expected uniform native meridian spline')
        use('spline', s.c.T * shape.unit_length**np.arange(3,-1,-1))
        use('shape', [v for obj in (shape, ref) for v in
                      (obj.exit/4-1.5, *obj.params[2:], obj.extra)])
        d.length, d.ref_length = shape.length, ref.length
        d.pressure = np.exp(chart.co[0])
        d.core_ratio = np.exp(-.5*np.exp(registration[0]))/self.cut
        d.epsilon = epsilon
        self.zero_dc = np.zeros(self.ids.shape)
        self.work = kernel()
        self._build_coordinates()
        self.metric = self._metric()

    def radial_mesh(self, n):
        return np.r_[np.linspace(0,self.cut,2*n+1)[:-1],np.linspace(self.cut,1,n+1)]

    def _build_coordinates(self):
        free = np.ones(self.chart_shape, bool)
        free[[0,-1]] = False
        free[:, -1] = False
        free[:, 0, :, :] = False
        free[1:-1, 0, 0, 1] = True
        free[:, :, [0,-1], 2] = False
        # One shared physical point at the center. The polar parametrization
        # of limiting directions has its regular even radial continuation.
        ix = np.full(self.chart_shape, -1, np.int64)
        ix[free] = np.arange(free.sum())
        ix[:, 0, :, 1] = ix[:, 0, :1, 1]
        ix[:, 0, :, 2] = ix[:, 1, :, 2]
        self.map, self.free = ix, free
        self.dim = int(free.sum())
        self.valid = ix >= 0

    def _metric(self):
        # Approximate diagonal weak elasticity on physical chart displacements.
        t,x,a=np.meshgrid(self.t,self.x,self.a,indexing='ij')
        b=self.sh.radius(t*self.sh.length)[0]
        r=x*b
        dx=(np.gradient(self.x)[None,:,None]*b)
        dz=self.sh.length*np.gradient(self.t)[:,None,None]
        da=np.gradient(self.a)[None,None,:]
        rad=np.maximum(r,dx/3)
        sn=np.maximum(np.sin(a),np.sin(da/3))
        vol=4*np.pi*rad**2*sn*dx*dz*da
        val=vol*(1/dx**2+1/dz**2+1/(rad*da)**2+1)
        scale=np.stack((b**2,np.full_like(t,self.sh.length**2),rad**2),axis=-1)
        return self.reduce(val[...,None]*scale)

    def reduce(self, array):
        return np.bincount(self.map[self.valid],weights=array[self.valid],minlength=self.dim)

    def pack(self, array):
        array = np.asarray(array)
        if array.shape != self.chart_shape or not np.all(np.isfinite(array)):
            raise ValueError('Incompatible solid displacement array')
        result = array[self.free].copy()
        if np.max(abs(self.unpack(result)-array))>5e-12:
            raise ValueError('Displacement violates the shared center or prescribed traces')
        return result

    def unpack(self, q):
        if np.shape(q) != (self.dim,) or not np.all(np.isfinite(q)):
            raise InvalidPlacement('Malformed material coordinates')
        out = np.zeros(self.chart_shape)
        out[self.valid] = np.asarray(q)[self.map[self.valid]]
        return out

    def seed(self, displacement, internal_amplitude=0.):
        """Interpolate a previous physical displacement and invert the wall chart."""
        if displacement.ndim==4:
            return self.pack(displacement)
        nk,nr,nd=displacement.shape
        if nd!=2 or (nr-1)%3:raise ValueError('Unrecognized meridional seed')
        old_t=np.linspace(0,1,nk)
        old_x=self.radial_mesh((nr-1)//3)
        t,x=np.meshgrid(self.t,self.x,indexing='ij')
        sampled=RegularGridInterpolator((old_t,old_x),displacement)((t,x))
        old_radius=x*self.sh.radius(t*self.sh.length)[0]+sampled[...,0]
        z=self.sh.length*t+sampled[...,1]
        width=self.sh.radius(z)[0]
        xp=(old_radius/width)[...,None]
        h=bump(z,self.sh.length)[0][...,None]
        p=(3*np.cos(self.a)**2-1)/2
        dest=old_radius[...,None]+h*xp**3*(self.epsilon+internal_amplitude*(1-xp))*p
        landed=np.broadcast_to(xp,dest.shape).copy()
        for _ in range(5):
            landed -= (landed*width[...,None]+self.epsilon*h*landed**3*p-dest)/(
                       width[...,None]+3*self.epsilon*h*landed**2*p)
        if np.max(abs(landed*width[...,None]+self.epsilon*h*landed**3*p-dest))>2e-12:
            raise InvalidPlacement('Failed to invert physical shape seed')
        out=np.zeros(self.chart_shape)
        out[...,0]=landed-x[...,None]
        out[...,1]=sampled[...,1,None]/self.sh.length
        return self.pack(out)

    def lower_bound_load(self, q, absolute):
        """Assemble the order load with all oriented contributions, including at χ=0.

        On an unoriented, exactly melted site the cost of reordering is
        undefined. Every oriented 4D extension has W >= 0 since each principal
        stretch satisfies (s²-1)/2 - log(s) >= 0. Positive reference weights
        and Q1 functions preserve this lower bound coefficient by coefficient.
        The affected global order nodes are returned in undefined_order_nodes.
        The energy itself remains the exact finite quadrature at the input order.
        """
        if not self.melted_labels:
            raise ValueError('The lower-bound load requires melted_labels=True')
        return self._evaluate(q, absolute, False, True, False, True)

    def evaluate(self, q, correction=None, *, diagnose=False, load=True, gradient=True):
        return self._evaluate(q, correction, diagnose, load, gradient, False)

    def _evaluate(self, q, correction, diagnose, load, gradient, lower_bound):
        offs=self.unpack(q)
        if self.melted_labels and (correction is None or np.any(self.buffers['base'] != 0)):
            raise ValueError('Melted-label evaluation requires the actual absolute order and zero background')
        dc=(self.zero_dc if correction is None else
            np.ascontiguousarray(np.asarray(correction)[self.ids],dtype=float))
        if not np.all(np.isfinite(dc)):raise ValueError('Nonfinite conversion field')
        if self.melted_labels and (np.any(dc<0) or np.any(dc>1)):
            raise ValueError('Absolute order outside [0,1]')
        self.data.offset, self.data.correction = pointer(offs),pointer(dc)
        self.data.diagnose=int(diagnose)
        native=np.zeros_like(offs) if gradient else None
        block=np.zeros(self.ids.shape) if load else None
        stats=np.zeros(30)
        melt=np.zeros(13) if self.melted_labels else None
        undefined=np.zeros(self.ids.shape) if self.melted_labels else None
        result=Results(pointer(native) if gradient else Pointer(),
                       pointer(block) if load else Pointer(),pointer(stats),
                       pointer(melt) if self.melted_labels else Pointer(),
                       pointer(undefined) if self.melted_labels else Pointer())
        code=self.work(ct.byref(self.data),ct.byref(result))
        if code==3:
            raise InvalidPlacement(f'Nonpositive J={melt[9]} at landed {melt[10:13].tolist()} '
                                   'outside a proven exactly zero-order neighborhood')
        if code:raise InvalidPlacement('Material left open order chart' if code==1 else 'Inverted solid')
        g=self.reduce(native) if gradient else None
        deposit=(np.bincount(self.ids.ravel(),weights=block.ravel(),minlength=self.model.size)
                 if load else None)
        info=dict(transfer_energy=stats[1],unit_order_energy=stats[2],
                  reference_volume=stats[3],quadrature_points=int(stats[4]))
        if self.melted_labels:
            count=int(melt[0])
            info.update(nonpositive_J_points=count, min_nonpositive_J=float(melt[1]) if count else None,
                        min_nonpositive_J_site=melt[2:5].tolist() if count else None,
                        min_nonpositive_J_reference_site=melt[5:8].tolist() if count else None,
                        nonpositive_J_reference_volume=float(melt[8]),
                        undefined_order_nodes=np.unique(self.ids[undefined!=0]).tolist(),
                        order_load_is_lower_bound=bool(lower_bound and count))
            if count:
                info['oriented_unit_order_energy']=info['unit_order_energy']
                info['unit_order_energy']=None
                if load and not lower_bound:
                    raise UnresolvedActivation(info)
        if diagnose:
            for name,ix in (('min_J',5),('max_J',9),('minimum_order',13),('min_stretch',17),
                            ('max_stretch',21),('maximum_order',25)):
                info[name]=stats[ix];info[name+'_site']=stats[ix+1:ix+4].tolist()
            info['twist_coefficient']=stats[29]
        if not np.isfinite(stats[0]) or (gradient and not np.all(np.isfinite(g))):
            raise InvalidPlacement('Nonfinite crystal energy or derivative')
        return float(stats[0]),g,deposit,info

    def residual(self,q,g):
        """Convert native generalized forces to spatial radial, axial, polar work."""
        off=self.unpack(q)
        force=np.zeros(self.chart_shape)
        force[self.free]=-g
        t,x,a=np.meshgrid(self.t,self.x,self.a,indexing='ij')
        xp,tp,angle=x+off[...,0], t+off[...,1], a+off[...,2]
        z=tp*self.sh.length
        b,bz=self.sh.radius(z)
        h,hz=bump(z,self.sh.length)
        P=(3*np.cos(angle)**2-1)/2
        R=b*xp+self.epsilon*h*xp**3*P
        Rx=b+3*self.epsilon*h*xp**2*P
        Rz=bz*xp+self.epsilon*hz*xp**3*P
        Ra=-3*self.epsilon*h*xp**3*np.cos(angle)*np.sin(angle)
        force[...,0]/=Rx
        force[...,1]=force[...,1]/self.sh.length-Rz*force[...,0]
        force[...,2]=np.divide(force[...,2]-Ra*force[...,0],R,
                              out=np.zeros_like(R),where=R>1e-12)
        force[~self.free]=0
        index=np.unravel_index(np.argmax(abs(force)),force.shape)
        coherent=force[...,:2].sum(axis=2)
        bulk_index=np.unravel_index(np.argmax(abs(coherent)),coherent.shape)
        last=force[-2,-2]
        return dict(material_force_max=float(abs(force[index])),material_force_site=list(map(int,index)),
                    angular_sum_force_max=float(abs(coherent[bulk_index])),
                    angular_sum_force_site=list(map(int,bulk_index)),
                    bond_force_max=float(abs(force[:,self.cut_node]).max()),
                    last_free_ring_angular_sum=coherent[-2,-2].tolist(),
                    last_free_ring_peak_per_direction=abs(last).max(axis=0).tolist(),
                    last_free_ring_reference=[float(self.t[-2]),float(self.x[-2])],
                    last_free_ring_landed_depth_fraction=[float(tp[-2,-2].min()),float(tp[-2,-2].max())])
