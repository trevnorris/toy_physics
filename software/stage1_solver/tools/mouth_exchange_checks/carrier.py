"""Independent triangular weak-flow solve with sheet resistance on the side and end cap."""
import sys,json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
sys.path.insert(0,'software/stage1_solver/tools')
from mouth_exchange_geometry import ExchangeShape

def solve(e,n,beta,outer=28.,roof=22.):
 sh=ExchangeShape(e);L=sh.length
 # The artificial radial seams in the feed and far bulk carry no interface cost.
 eta=np.linspace(0,1,n+1);st=np.linspace(0,1,2*n+1);xt=np.expm1(3.5*st)/np.expm1(3.5)
 zf=-eta[::-1];za=L*np.linspace(0,1,2*n+1);zb=L+(roof-L)*xt
 nodes={}; coords=[]; R=[];C=[];K=[]
 def point(r,z,kind):
  if z < -1e-10: kind=0
  elif z > L+1e-10: kind=1
  elif z < 1e-10:
   if kind=='F' or kind=='I' or (beta==0 and r<sh.neck+1e-8):kind=0
   else:kind=1
  elif kind=='I' and beta>0:kind=2
  else:kind=1
  key=round(r,10),round(z,10),kind
  if key not in nodes:
   nodes[key]=len(coords);coords.append((r,z,kind))
  return nodes[key]
 def block(level,section):
  z= {'F':zf,'I':za,'E':za,'B':zb}[level]
  anchors=sh.radius(z)[0]
  rr=anchors[:,None]*eta if section==0 else anchors[:,None]+(outer-anchors[:,None])*xt
  zz=np.broadcast_to(z[:,None],rr.shape)
  idx=np.fromiter((point(r,h,level) for r,h in zip(rr.flat,zz.flat)),int,rr.size).reshape(rr.shape)
  aa,ba,ab,bb=[v.ravel() for v in (idx[:-1,:-1],idx[:-1,1:],idx[1:,:-1],idx[1:,1:])]
  # Alternating triangle orientations on each rectangle
  a,b0=(q.ravel() for q in np.indices((len(z)-1,rr.shape[1]-1)))
  s=(a+b0)%2==0
  tri=np.vstack((np.c_[aa[s],ba[s],bb[s]],np.c_[aa[s],bb[s],ab[s]],
                 np.c_[aa[~s],ba[~s],ab[~s]],np.c_[ab[~s],ba[~s],bb[~s]]))
  P=np.array(coords)[tri,:2];rad=P[:,:,0]
  det=(P[:,1,0]-P[:,0,0])*(P[:,2,1]-P[:,0,1])-(P[:,2,0]-P[:,0,0])*(P[:,1,1]-P[:,0,1]);assert min(det)>0
  area_r2=det*(.5/6)*(np.sum(rad*rad,axis=1)+rad[:,0]*rad[:,1]+rad[:,1]*rad[:,2]+rad[:,2]*rad[:,0])
  delta=np.roll(P,-1,axis=1)-np.roll(P,-2,axis=1)
  grad=np.stack((delta[:,:,1],-delta[:,:,0]),axis=2)/det[:,None,None]
  val=4*np.pi*area_r2[:,None,None]*(grad@grad.transpose(0,2,1))
  R.append(np.repeat(tri,3,axis=1).ravel());C.append(np.tile(tri,(1,3)).ravel());K.append(val.ravel())
  return idx
 block('F',0);block('F',1);inside=block('I',0);outside=block('E',1);above=block('B',0);block('B',1)
 if beta:
  for one,two in ((inside[:,-1],outside[:,0]),(inside[-1,:],above[0,:])):
   pair=np.c_[one[:-1],one[1:],two[:-1],two[1:]]
   edges=np.array(coords)[pair[:,:2],:2]; ds=np.linalg.norm(np.diff(edges,axis=1)[:,0,:],axis=1)
   mat=np.zeros((len(ds),4,4))
   for t,a in zip(*np.polynomial.legendre.leggauss(3)):
    t=(t+1)/2;rad=edges[:,0,0]*(1-t)+edges[:,1,0]*t;u=np.array([1-t,t,t-1,-t])
    mat+=np.pi*2*a*(ds*rad*rad/beta)[:,None,None]*np.outer(u,u)
   R.append(np.repeat(pair,4,axis=1).ravel());C.append(np.tile(pair,(1,4)).ravel());K.append(mat.ravel())
 mat=coo_matrix((np.hstack(K),(np.hstack(R),np.hstack(C))),shape=(len(coords),len(coords))).tocsr()
 rad,z,kind=np.array(coords).T
 side=abs(rad-outer)<1e-9
 feed=side&(kind==0)
 far=(side&(kind==1))|(z>roof-1e-8)
 u=far.astype(float);free=~(feed|far)
 u[free]=spsolve(mat[free][:,free],-(mat@u)[free])
 flux=mat@u
 assert max(abs(flux[free]))<1e-9 and abs(flux[far].sum()+flux[feed].sum())<1e-8
 return u@flux*.5*.4510084282857874**2

def main(n):
 for beta in [0.,.04/6]:
  table={a:solve(a,n,beta) for a in [-.01,-.005,0.,.005,.01]}
  print(json.dumps(dict(n=n,beta=beta,energies=table,first=[(table[a]-table[-a])/(2*a) for a in [.01,.005]],second=[(table[a]+table[-a]-2*table[0])/a**2 for a in [.01,.005]])),flush=True)

if __name__=='__main__':
 main(int(sys.argv[1]))
