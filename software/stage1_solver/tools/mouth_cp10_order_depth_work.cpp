// At stationary native images, decompose the depth derivative into geometric
// work and work from moving order knots. The additional derivative of the
// resampled spatial labels belongs to the full signed calculation.
#include "mouth_cp10_nested_energy_depth.cpp"

namespace {
struct ActivePrisms : OriginalPrisms {
    const unsigned char *mask;
    int nr, na;
    bool include(int j,int i,int k,int) const {return mask[(j*nr+i)*na+k];}
};
struct OrderWork {
    const SpatialLabels *grids;
    const double *rate;
    long double (*rows)[3];
    long long &count;
    int operator()(const Data &d,Results &,const Profile &,const TriSlice &s,
                   double t,double a,double x,double w,bool core,double &den,double &denx) const {
        den=denx=0.;
        Vec p=s.X;
        Mat m=s.M;
        for(int c=0;c<3;c++){p[c]+=s.M[c][0]*(x-s.x0);m[c][2]+=s.Max[c]*(x-s.x0);}
        const auto r=basis(d.field_x,p[0]),z=basis(d.field_t,p[1]),b=basis(d.field_a,p[2]);
        double chi=0.,slope=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            const double c=d.correction[field_id(d,r.i+i,b.i+k,z.i+j)];
            chi+=r.w[i]*z.w[j]*b.w[k]*c;
            slope+=r.w[i]*z.v[j]*b.w[k]*c;
        }
        const double dchi=-slope*(z.w[0]*rate[z.i]+z.w[1]*rate[z.i+1]);
        if(chi==0. && dchi==0.)return 0;
        if(!(chi>0.))return 1;
        const double mj=determinant(m);
        const auto wall=radius(d,p);
        const Mat frame{wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
        const double hs=wall.R*std::sin(p[2]),dv=w*4*pi*determinant(frame)*hs*mj;
        if(!(dv>0.))return 2;
        const auto old=local_work(d,{x,t,a},m,mj,frame,hs,dv,chi,dchi,core);
        if(!std::isfinite(old.stress+old.order))return 3;
        rows[0][0]+=old.stress; rows[0][1]+=old.order;
        for(int g=0;g<3;g++){
            Vec xi{};
            Mat grad{};
            grids[g].interpolate(p,xi,grad);
            const double det=determinant(grad);
            if(!(det>0.))return 4;
            Mat fwd=transpose(cofactor(grad));
            for(auto &row:fwd)for(double &v:row)v/=det;
            const bool side=xi[0]<=d.x.v[(d.x.n-1)*2/3];
            const auto projected=local_work(d,xi,fwd,1./det,frame,hs,dv,chi,dchi,side);
            if(!std::isfinite(projected.stress+projected.order))return 5;
            rows[g+1][0]+=projected.stress; rows[g+1][1]+=projected.order;
            rows[g+1][2]+=std::abs(projected.order-old.order);
        }
        ++count;
        return 0;
    }
};
}
extern "C" int order_work(Data *d,Results *out,const Grid *axes,const double **labels,
        const double *rate,const unsigned char *mask,double *row,long long *count) {
    SpatialLabels grids[3];
    for(int i=0;i<3;i++)grids[i]={axes[3*i],axes[3*i+1],axes[3*i+2],labels[i]};
    long double totals[4][3]{};
    ActivePrisms selected{{},mask,d->x.n-1,d->a.n-1};
    const int rc=integrate_solid(d,out,selected,OrderWork{grids,rate,totals,*count});
    for(int i=0;i<4;i++)for(int j=0;j<3;j++)row[3*i+j]=double(totals[i][j]);
    return rc;
}
