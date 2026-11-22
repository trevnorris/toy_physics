"""Independent waist-clamped P1 triangular shear cavity; compare actual pole work/flux controls separately."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import spherical_jn, spherical_yn
sys.path.insert(0,'software/stage1_solver/tools')
from mouth_exchange_geometry import ExchangeShape
from mouth_exchange_interior import InteriorShape
from mouth_exchange_shear import OutgoingSlab

def clamped(shape,n):
 z=np.linspace(0,shape.length,2*n+1);x=np.linspace(0,1,n+1); rr=shape.radius(z)[0][:,None]*x;zz=np.broadcast_to(z[:,None],rr.shape)
 rq,wq=np.polynomial.legendre.leggauss(4);rq=(rq+1)/2;wq=wq/2
 pts=np.array([rr.ravel(),zz.ravel()]).T
 ids=np.arange(rr.size).reshape(rr.shape)
 aa,ba,ab,bb=[i.ravel() for i in (ids[:-1,:-1],ids[:-1,1:],ids[1:,:-1],ids[1:,1:])]
 a,b=(v.ravel() for v in np.indices((len(z)-1,len(x)-1)));s=(a+b)%2==0
 tri=np.vstack((np.c_[aa[s],ba[s],bb[s]],np.c_[aa[s],bb[s],ab[s]],np.c_[aa[~s],ba[~s],ab[~s]],np.c_[ab[~s],ba[~s],bb[~s]]))
 P=pts[tri];d=np.roll(P,-1,axis=1)-np.roll(P,-2,axis=1)
 det=(P[:,1,0]-P[:,0,0])*(P[:,2,1]-P[:,0,1])-(P[:,2,0]-P[:,0,0])*(P[:,1,1]-P[:,0,1])
 grad=np.stack((d[:,:,1],-d[:,:,0]),axis=2)/det[:,None,None];assert min(det)>0
 Km=np.zeros((len(tri),3,3));Mm=np.zeros_like(Km);l=4
 for u,wu in zip(rq,wq):
  for t,wt in zip(rq,wq):
   N=np.array([(1-u)*(1-t),u,t*(1-u)])
   r=P[:,:,0]@N
   wt2=det*wu*wt*(1-u)*r*r
   curl=2*N+r[:,None]*grad[:,:,0]; dz=grad[:,:,1]*r[:,None]
   Km+=wt2[:,None,None]*(l*(l+1)*np.outer(N,N)+curl[:,:,None]*curl[:,None,:]+dz[:,:,None]*dz[:,None,:])
   Mm+=(wt2*r*r)[:,None,None]*np.outer(N,N)
 free=(pts[:,0]>1e-12)&(pts[:,1]>1e-12)
 row,col=np.repeat(tri,3,axis=1).ravel(),np.tile(tri,(1,3)).ravel()
 def matrix(a):
  return coo_matrix((a.ravel(),(row,col)),shape=(len(pts),len(pts))).tocsr()[free][:,free].tocsc()
 K,M=matrix(Km),matrix(Mm)
 lam=eigsh(K,M=M,k=2,sigma=0,v0=np.ones(K.shape[0]),return_eigenvectors=False,tol=1e-11)
 return np.sqrt(np.sort(lam))

def dtn():
 R=3.
 port=OutgoingSlab(np.c_[np.full(25,R),np.linspace(-1,0,25)],R)
 om=4.20802-3.9e-8j
 k=np.sqrt(om**2-port.cutoff**2)
 k=np.where((port.cutoff>om.real)&(k.imag<0),-k,k)
 # Real or first weakly evanescent channels only, to avoid loss of relative accuracy in j+i y.
 ss=k[:4]*R;l=4
 val=(spherical_jn(l,ss)+1j*spherical_yn(l,ss));grad=(spherical_jn(l,ss,derivative=True)+1j*spherical_yn(l,ss,derivative=True))
 t=k[:4]*grad/val+1/R
 e=np.abs(t-port.traction(om)[:4]);assert max(e[:3])<1e-8
 print(json.dumps(dict(DtN_disagreement=e.tolist())))
def run(n):
 for key,shape in [('exchange',ExchangeShape),('interior',InteriorShape)]:
  h=.005
  poles=np.array([clamped(shape(q),n) for q in [-h,0.,h]])
  print(json.dumps(dict(n=n,kind=key,clamped_modes=poles.tolist(),slope=(poles[2,0]-poles[0,0])/(2*h),curvature=(poles[2,0]+poles[0,0]-2*poles[1,0])/h**2)),flush=True)
 dtn()
if __name__=='__main__':run(int(sys.argv[1]))
