"""One conforming angular edge split of the wound radius-four crystal.

All original coordinates are hierarchical coarse coordinates. Three new
coordinates are the material-chart offset at the inserted site relative to
linear interpolation of its two parents. Fixed parent labels remain identical.
"""
import copy
import ctypes as ct
from functools import lru_cache
from hashlib import sha256
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from mouth_free_shape_solid import Data, Grid, InvalidPlacement, Pointer, Results, pointer
from mouth_joined_angular_work import bump


@lru_cache(maxsize=4)
def refined_kernel(folder):
    source=Path(__file__).with_suffix('.cpp')
    parent=source.with_name('mouth_free_shape_solid.cpp')
    digest=sha256(source.read_bytes()+parent.read_bytes()).hexdigest()
    folder=Path(folder)
    folder.mkdir(exist_ok=True,parents=True)
    library=folder/(digest+'.so')
    if not library.exists():
        fd,stage=tempfile.mkstemp(suffix='.so',dir=folder)
        os.close(fd)
        try:
            subprocess.run(['g++','-O3','-std=c++17','-fPIC','-shared',str(source),'-o',stage],check=True)
            os.replace(stage,library)
        finally:
            Path(stage).unlink(missing_ok=True)
    function=ct.CDLL(str(library)).refined_solid_eval
    function.argtypes=[ct.POINTER(Data),ct.POINTER(Results),Pointer,ct.c_double,
                       *([ct.c_int]*6),Pointer]
    function.restype=ct.c_int
    return function


class RefinedSolid:
    """The existing energy and its gradients integrated on the edge star."""
    vertex=(5,18,0)

    def __init__(self,parent,build_dir=None):
        if parent.chart_shape!=(25,25,7,3) or parent.sh.neck!=4.:
            raise ValueError('This refinement is defined on the R4 coarse 24 x 24 preparation')
        self.parent=parent
        if build_dir is None:
            build_dir=Path.cwd()/'_scratch/s11c/r4_refined_solid_build'
        self.native=refined_kernel(str(Path(build_dir).resolve()))
        edges=parent.buffers['angle_edges']
        a0,a1=parent.a[:2]
        inside=edges[(edges>a0)&(edges<a1)]
        if len(inside)==0:
            raise ValueError('No interior angular integration edge to separate the two child prisms')
        self.seam=float(inside[np.argmin(abs(inside-0.5*(a0+a1)))])
        self.seam_fraction=(self.seam-a0)/(a1-a0)
        self.dim=parent.dim+3
        # A positive preconditioning metric only. Use the adjacent positive
        # angular node since the polar component at the axis is constrained.
        metric_ids=parent.map[5,18,1,:]
        if np.any(metric_ids<0):
            raise ValueError('Missing neighboring unconstrained chart component')
        self.metric=np.r_[parent.metric,parent.metric[metric_ids]]

    def __getattr__(self,name):
        return getattr(self.parent,name)

    def split(self,q):
        q=np.asarray(q,dtype=np.float64)
        if q.shape!=(self.dim,) or not np.all(np.isfinite(q)):
            raise InvalidPlacement('Missing or malformed refined solid coordinates')
        return q[:-3],np.ascontiguousarray(q[-3:])

    def pack(self,offsets,surplus=None):
        return np.r_[self.parent.pack(offsets),np.zeros(3) if surplus is None else surplus]

    def unpack(self,q):
        return self.parent.unpack(self.split(q)[0])

    def _evaluate(self,q,absolute,diagnose,load,gradient,lower_bound,region,refined,edges=None):
        q_parent,surplus=self.split(q)
        if not refined and np.any(surplus):
            raise ValueError('Parent quadrature can only represent zero child surplus')
        extra_grad=np.zeros(3)
        j=i=k=half=-1
        if region is None:
            mode=0
        elif isinstance(region,str) and region=='star':
            mode=1
        else:
            mode=2
            j,i,k,half=map(int,region)
            if not (0<=j<24 and 0<=i<24 and 0<=k<6 and 0<=half<2):
                raise ValueError('Invalid reference prism')
        # Reuse the original, independently pinned Python ABI and all of its
        # chart constraints, absolute-order checks and global load assembly.
        # Give each invocation its own native readout object and data pointers.
        evaluator=copy.copy(self.parent)
        evaluator.data=Data.from_buffer_copy(self.parent.data)
        # If the caller asks for a separate angular integration check, use
        # the native quadrature at each added panel with its original weights.
        if edges is not None:
            evaluator.data.angle_edges=Grid(len(edges),pointer(edges))
        evaluator.work=lambda data,result:self.native(
            data,result,pointer(surplus),self.seam,refined,mode,j,i,k,half,pointer(extra_grad))
        e,g,dep,info=evaluator._evaluate(q_parent,absolute,diagnose,load,gradient,lower_bound)
        return e,np.r_[g,extra_grad] if gradient else None,dep,info

    def evaluate(self,q,absolute,*,diagnose=False,load=True,gradient=True):
        return self._evaluate(q,absolute,diagnose,load,gradient,False,None,True)

    def lower_bound_load(self,q,absolute):
        return self._evaluate(q,absolute,False,True,False,True,None,True)

    def evaluate_region(self,q,absolute,region,*,refined=True,load=False,gradient=True,
                        diagnose=False,angular_subdivisions=1):
        """Integrate an individual parent or six-prism star with its true child panel.

        Optional additional subdivisions re-integrate the angular quadrature
        independently. An undefined order derivative in melted material is
        returned as a marked lower bound, as in the parent's lower_bound_load.
        """
        if not isinstance(angular_subdivisions,int) or angular_subdivisions<1:
            raise ValueError('Expected a positive integer angular subdivision')
        edges=None
        if angular_subdivisions>1:
            edges=self.parent.buffers['angle_edges']
            edges=np.r_[
                (edges[:-1,None]+np.diff(edges)[:,None]*
                 (np.arange(angular_subdivisions)/angular_subdivisions)).ravel(),edges[-1]]
        return self._evaluate(q,absolute,diagnose,load,gradient,load,region,refined,edges)

    def residual(self,q,g):
        qp,surplus=self.split(q)
        g=np.asarray(g)
        if g.shape!=(self.dim,):
            raise ValueError('Incompatible refined solid gradient')
        j,i,k=self.vertex
        # Hierarchical coarse derivatives include hat-interpolated midpoint
        # reactions. Recover independent nodal reactions before conversion.
        gp=g[:-3].copy()
        for end,weight in ((k,1.-self.seam_fraction),(k+1,self.seam_fraction)):
            ids=self.parent.map[j,i,end]
            free=ids>=0
            np.add.at(gp,ids[free],-weight*g[-3:][free])
        out=self.parent.residual(qp,gp)
        ring_offsets=self.parent.unpack(qp)[j,i]
        offsets=ring_offsets[k:k+2]
        offset=(1.-self.seam_fraction)*offsets[0]+self.seam_fraction*offsets[1]+surplus
        locations=np.column_stack((np.full(len(ring_offsets),self.parent.x[i]),
                                   np.full(len(ring_offsets),self.parent.t[j]),self.parent.a))
        locations=np.vstack((locations+ring_offsets,
                             np.array([self.parent.x[i],self.parent.t[j],self.seam])+offset))
        x,t,a=locations.T
        shape=self.parent.sh
        z=shape.length*t
        b,bz=shape.radius(z)
        h,hz=bump(z,shape.length)
        p=(3.*np.cos(a)**2-1.)/2.
        eps=self.parent.epsilon
        r=b*x+eps*h*x**3*p
        rx=b+3.*eps*h*x*x*p
        rz=bz*x+eps*hz*x**3*p
        ra=-3.*eps*h*x**3*np.cos(a)*np.sin(a)
        ring_ids=self.parent.map[j,i]
        valid=ring_ids>=0
        nodal=np.zeros_like(ring_offsets)
        nodal[valid]=-gp[ring_ids[valid]]
        forces=np.vstack((nodal,-g[-3:]))
        forces[:,0]/=rx
        forces[:,1]=forces[:,1]/shape.length-rz*forces[:,0]
        forces[:,2]=(forces[:,2]-ra*forces[:,0])/r
        forces[:-1][~self.parent.free[j,i]]=0.
        force=forces[-1]
        out['child_node_force']=force.tolist()
        out['child_node_reference']=[self.parent.x[i],self.parent.t[j],self.seam]
        component=int(np.argmax(abs(force)))
        if abs(force[component])>out['material_force_max']:
            out['material_force_max']=float(abs(force[component]))
            out['material_force_site']=[j,i,self.seam_fraction,component]
        # Update the coherent (radial, axial) sum. Recomputing the baseline
        # without this ring handles even a previous maximum at the new site.
        if out['angular_sum_force_site'][:2]==[j,i]:
            without=gp.copy()
            without[ring_ids[valid]]=0.
            unaffected=self.parent.residual(qp,without)
            out['angular_sum_force_max']=unaffected['angular_sum_force_max']
            out['angular_sum_force_site']=unaffected['angular_sum_force_site']
        coherent=forces[:,:2].sum(axis=0)
        component=int(np.argmax(abs(coherent)))
        if abs(coherent[component])>out['angular_sum_force_max']:
            out['angular_sum_force_max']=float(abs(coherent[component]))
            out['angular_sum_force_site']=[j,i,component]
        return out
