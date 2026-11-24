// Frozen native quadrature on the unfolded reference. Order is evaluated on its
// own CP10 grid in both banks; only the depth coordinates of the label samples vary.
#include "mouth_cp10_bucket_work.cpp"

namespace {
struct FrozenShiftWork {
    const SpatialLabels *labels;
    const double *rate;
    const unsigned char *front;
    long double (*sums)[7][2];
    long long &points;
    int operator()(const Data &d,Results &,const Profile &,const TriSlice &s,double t,double a,
                   double x,double w,bool core,double &den,double &denx) const {
        den=denx=0.;
        Vec p=s.X; Mat native=s.M;
        for(int c=0;c<3;c++){p[c]+=s.M[c][0]*(x-s.x0);native[c][2]+=s.Max[c]*(x-s.x0);}
        auto xr=basis(d.field_x,p[0]),zt=basis(d.field_t,p[1]),pa=basis(d.field_a,p[2]);
        double chi=0.,slope=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            double v=xr.w[i]*pa.w[k]*d.correction[field_id(d,xr.i+i,pa.i+k,zt.i+j)];
            chi+=v*zt.w[j];slope+=v*zt.v[j];
        }
        const double work=slope*(zt.w[0]*rate[zt.i]+zt.w[1]*rate[zt.i+1]);
        if(chi==0. && work==0.)return 0;
        if(!(chi>0.))return 1;
        const auto wall=radius(d,p);
        const Mat frame{wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
        const double hs=wall.R*std::sin(p[2]), dj=determinant(native);
        const double dv=w*4*pi*determinant(frame)*hs*dj;
        if(!(dv>0.))return 2;
        int region=front[(xr.i*(d.field_t.n-1)+zt.i)*(d.field_a.n-1)+pa.i];
        for(int k=0;k<7;k++){
            Vec label{x,t,a};Mat m=native;double det=dj;bool side=core;
            if(k){
                label={};Mat g{};labels[k-1].interpolate(p,label,g);
                const double di=determinant(g);
                if(!(di>0.))return 3;
                m=transpose(cofactor(g));for(auto &r:m)for(double &v:r)v/=di;
                det=1./di; side=label[0]<=d.x.v[(d.x.n-1)*2/3];
            }
            auto q=measure(d,label,m,det,frame,hs,dv,chi,work,side);
            if(!std::isfinite(q.energy+q.order))return 4;
            sums[region][k][0]+=q.energy;sums[region][k][1]+=q.order;
        }
        ++points;return 0;
    }
};
}
extern "C" int frozen_shift(Data *d,Results *r,const Grid *axes,const double **nodes,const double *rate,
                            const unsigned char *active,const unsigned char *front,double *out,long long *count){
    SpatialLabels labels[6];
    for(int i=0;i<6;i++)labels[i]={axes[3*i],axes[3*i+1],axes[3*i+2],nodes[i]};
    long double sums[2][7][2]{};
    ActivePrisms selected{{},active,d->x.n-1,d->a.n-1};
    int rc=integrate_solid(d,r,selected,FrozenShiftWork{labels,rate,front,sums,*count});
    for(auto &reg:sums)for(auto &bank:reg)for(auto v:bank)*out++=double(v);
    return rc;
}
