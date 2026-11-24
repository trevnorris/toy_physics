// Common-point stored energy and the five motions that leave all chart knots fixed.
// The independently integrated signed depth derivative handles the sixth motion.
#include "mouth_cp10_order_depth_work.cpp"

namespace {
struct WallVariation { double radius; Vec grad; };
using Variations=std::array<WallVariation,5>;

Variations variations(const Data &d,const Vec &p) {
    double x=p[0],t=p[1],a=p[2],prod=t*(1-t);
    const double exit=d.shape[0]+1.5;
    const double val[]={exit*t,4*prod,4*prod*(2*t-1),16*prod*prod};
    const double der[]={exit,4*(1-2*t),4*(-6*t*t+6*t-1),32*prod*(1-2*t)};
    Variations out{};
    for(int i=0;i<4;i++)out[i]={x*val[i],{val[i],x*der[i],0.}};
    const double y=(t-.06)/.56, sn=std::sin(a), cs=std::cos(a), p2=(3*cs*cs-1)/2;
    if(y>0 && y<1) {
        const double f=16*y*y*(1-y)*(1-y), ft=32*y*(1-y)*(1-2*y)/.56;
        out[4]={x*x*x*f*p2,{3*x*x*f*p2,x*x*x*ft*p2,-3*x*x*x*f*cs*sn}};
    }
    return out;
}

struct ShapeWork {
    const SpatialLabels *grids;
    const Data *signed_shapes;
    long double (*rows)[6];
    double &worst_fd;
    long long &count;
    int operator()(const Data &d,Results &,const Profile &,const TriSlice &s,double t,double a,
                   double x,double w,bool core,double &den,double &denx) const {
        den=denx=0;
        Vec p=s.X;
        Mat m=s.M;
        for(int c=0;c<3;c++){p[c]+=m[c][0]*(x-s.x0);m[c][2]+=s.Max[c]*(x-s.x0);}
        const auto bx=basis(d.field_x,p[0]),bt=basis(d.field_t,p[1]),ba=basis(d.field_a,p[2]);
        double chi=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++)
            chi+=bx.w[i]*bt.w[j]*ba.w[k]*d.correction[field_id(d,bx.i+i,ba.i+k,bt.i+j)];
        if(chi==0.)return 0;
        const auto wall=radius(d,p);
        const auto v=variations(d,p);
        const Mat frame{wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
        const double hs=wall.R*std::sin(p[2]),mj=determinant(m),dv=w*4*pi*determinant(frame)*hs*mj;
        if(!(chi>0 && dv>0))return 1;
        const bool check=(count++%100003)==0;
        for(int level=0;level<4;level++){
            Vec xi{x,t,a};
            Mat forward=m;
            double dj=mj;
            bool side=core;
            if(level>0){
                xi={}; Mat g{};
                grids[level-1].interpolate(p,xi,g);
                const double det=determinant(g);
                if(!(det>0))return 2;
                dj=1./det;forward=transpose(cofactor(g));
                for(auto &row:forward)for(double &z:row)z/=det;
                side=xi[0]<=d.x.v[(d.x.n-1)*2/3];
            }
            const Ref ref=reference(d,xi[0],xi[1],profile(d,xi[1],true),side);
            const Mat right=product(forward,ref.inv), F=product(frame,right);
            const double hoop=hs/(ref.S*std::sin(xi[2]));
            const double J=determinant(frame)*dj*hoop/(ref.det2*ref.S);
            if(!(J>0 && std::isfinite(J)))return 3;
            const double W=.5*(norm2(F)+hoop*hoop-4.)-std::log(J)+.25*(J-1.)*(J-1.);
            const double weight=chi*dv/J, volume_stress=.5*J*(J-1.)-1.;
            rows[level][0]+=weight*W;
            Vec conjugate;
            for(int c=0;c<3;c++)conjugate[c]=dot(F[0],right[c]);
            conjugate[0]+=volume_stress/wall.grad[0];
            const double radial=(dot(F[2],F[2])+hoop*hoop+2*volume_stress)/wall.R;
            for(int mode=0;mode<5;mode++){
                const double dw=dot(conjugate,v[mode].grad)+radial*v[mode].radius;
                rows[level][mode+1]-=weight*dw;
                if(!check)continue;
                // Central difference at scattered actual Gauss points, using signed production
                // geometry. Divide by the same reference volume before checking the stress.
                double es[2];
                for(int b=0;b<2;b++){
                    const auto &test=signed_shapes[2*mode+b];
                    const auto r=radius(test,p);
                    const Mat C{r.grad,Vec{0.,test.length,0.},Vec{0.,0.,r.R}};
                    const double hh=r.R*std::sin(p[2]), dvol=w*4*pi*determinant(C)*hh*mj;
                    es[b]=constitutive(test,xi,forward,C,dj,hh,dvol,chi,side).weighted/weight;
                }
                const double fd=(es[1]-es[0])/(2.e-5);
                worst_fd=std::max(worst_fd,std::abs(fd-dw)/(1+std::abs(dw)));
            }
        }
        return 0;
    }
};
}

extern "C" int unfolded_shape_work(Data *d,Results *out,const Grid *axes,const double **labels,
        const unsigned char *mask,const Data *probes,double *result,double *fd,long long *count){
    SpatialLabels grids[3];
    for(int k=0;k<3;k++)grids[k]={axes[3*k],axes[3*k+1],axes[3*k+2],labels[k]};
    long double rows[4][6]{};
    ActivePrisms selected{{},mask,d->x.n-1,d->a.n-1};
    const int rc=integrate_solid(d,out,selected,ShapeWork{grids,probes,rows,*fd,*count});
    for(int k=0;k<4;k++)for(int j=0;j<6;j++)result[6*k+j]=double(rows[k][j]);
    return rc;
}
