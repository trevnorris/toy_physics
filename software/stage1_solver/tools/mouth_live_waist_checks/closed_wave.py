"""Independent triangular curl and mass forms on an ordered, closed 4D cavity."""


import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

sys.path.insert(0, 'software/stage1_solver/tools/mouth_lip_checks')
from shape_check import Shape

def arrays(shape, n):
    z = np.linspace(0., shape.length, 2*n+1)
    x = np.linspace(0., 1., n+1)
    rad = shape.radius(z)[0]
    coords = np.stack(np.broadcast_arrays(rad[:,None]*x,z[:,None]),-1).reshape(-1,2)
    ids = np.arange(len(coords)).reshape(len(z),len(x))
    ll,lr,ul,ur=[g.ravel() for g in (ids[:-1,:-1],ids[:-1,1:],ids[1:,:-1],ids[1:,1:])]
    i,j=np.indices((len(z)-1,len(x)-1)); ev=(i+j).ravel()%2==0
    tri=np.vstack((np.c_[ll[ev],lr[ev],ur[ev]],np.c_[ll[ev],ur[ev],ul[ev]],
                   np.c_[ll[~ev],lr[~ev],ul[~ev]],np.c_[ul[~ev],lr[~ev],ur[~ev]]))
    P=coords[tri]
    v=P[:,1]-P[:,0];w=P[:,2]-P[:,0]
    det=v[:,0]*w[:,1]-v[:,1]*w[:,0]; assert min(det)>0
    d=np.roll(P,-1,axis=1)-np.roll(P,-2,axis=1)
    grad=np.stack((d[:,:,1],-d[:,:,0]),axis=2)/det[:,None,None]
    gauss,wt=np.polynomial.legendre.leggauss(3); gauss=(gauss+1)/2;wt=wt/2
    # Duffy map integrates the weighted radial products on physical straight triangles.
    q=[]
    for u,p in zip(gauss,wt):
        for v,w in zip(gauss,wt):
            N=np.array([(1-u)*(1-v),u,(1-u)*v])
            r=P[:,:,0]@N
            weight=p*w*(1-u)*det*r*r
            q.append((N,r,weight))
    return coords,ids,tri,grad,q
def matrix(local,tri,N):
    row=np.broadcast_to(tri[:,:,None],local.shape)
    col=np.broadcast_to(tri[:,None,:],local.shape)
    return coo_matrix((local.ravel(),(row.ravel(),col.ravel())),shape=(N,N)).tocsr()
def mode(data,chi,source=None,index=0):
    xy,ids,tri,grad,quad=data
    cut=(xy[:,1]>0)&(xy[:,0]>0); K,M,Z,DK,DM=(np.zeros((len(tri),3,3)) for _ in range(5))
    for N,r,vol in quad:
        angle=20*np.outer(N,N)
        curl=2*N+r[:,None]*grad[:,:,0];dw=r[:,None]*grad[:,:,1]
        z=dw[:,:,None]*dw[:,None,:]
        k=angle+curl[:,:,None]*curl[:,None,:]+z
        m=(r*r)[:,None,None]*np.outer(N,N)
        material=(chi[tri]@N)
        w=(vol*material)[:,None,None]
        K+=w*k;M+=w*m;Z+=w*z
        if source is not None:
            dw=(vol*(source[tri]@N))[:,None,None]
            DK+=dw*k;DM+=dw*m
    make=lambda a:matrix(a,tri,len(xy))[cut][:,cut]
    K,M,Z,DK,DM=map(make,(K,M,Z,DK,DM))
    vals,vecs=eigsh(K,M=M,k=3,sigma=0.,v0=np.ones(K.shape[0]),tol=1e-11)
    omega=vals[index]**.5;v=vecs[:,index]
    return omega,float(v@(Z@v)/vals[index]),float(v@(DK@v-vals[index]*(DM@v))/(2*omega)),np.sqrt(vals).tolist()
