// Signed log-depth work at the saved physical material points. The reconstructed
// xi(p) function remains fixed in the normalized physical chart p, while the
// original depth-dependent order interpolation knots move. This is a different
// virtual motion from inverting again at every displaced spatial knot.
#include "mouth_cp10_nested_energy.cpp"

namespace {
struct DepthWork {double stress,order,J;};
DepthWork local_work(const Data &d,const Vec &xi,const Mat &forward,double chart_J,const Mat &frame,
                     double azimuth_radius,double dv,double chi,double dchi,bool core){
    auto ref=reference(d,xi[0],xi[1],profile(d,xi[1],true),core);
    Mat F=product(frame,product(forward,ref.inv));
    const double hoop=azimuth_radius/(ref.S*std::sin(xi[2]));
    const double J=determinant(frame)*chart_J*hoop/(ref.det2*ref.S);
    if(!(J>0 && std::isfinite(J)))return {NAN,NAN,J};
    const double W=.5*(norm2(F)+hoop*hoop-4.)-std::log(J)+.25*(J-1)*(J-1);
    const double dW=dot(F[1],F[1])-1.+.5*J*(J-1);
    return {-dv/J*chi*dW,-dv/J*dchi*W,J};
}
struct DepthIntegrand {
    const SpatialLabels &labels;
    const double *tdot;
    long double *v;
    int operator()(const Data &d,Results &out,const Profile &init,const TriSlice &s,
                   double t,double a,double x,double w,bool core,double &den,double &denx)const {
        const int rc=sample(d,out,init,s,t,a,x,w,core,den,denx);
        if(rc)return rc;
        Vec p=s.X;Mat m=s.M;
        for(int c=0;c<3;c++){p[c]+=s.M[c][0]*(x-s.x0);m[c][2]+=s.Max[c]*(x-s.x0);}
        const auto r=basis(d.field_x,p[0]),z=basis(d.field_t,p[1]),b=basis(d.field_a,p[2]);
        double chi=0.,chit=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            const double c=d.correction[field_id(d,r.i+i,b.i+k,z.i+j)];
            chi+=r.w[i]*z.w[j]*b.w[k]*c;
            chit+=r.w[i]*z.v[j]*b.w[k]*c;
        }
        const double dchi=-chit*(z.w[0]*tdot[z.i]+z.w[1]*tdot[z.i+1]);
        if(chi==0.){if(dchi!=0.)return 6;return 0;}
        Vec y{};Mat g{};
        labels.interpolate(p,y,g);
        const double dg=determinant(g),mj=determinant(m);
        if(!(dg>0 && mj>0 && std::isfinite(dg)))return 4;
        Mat inv=transpose(cofactor(g)); for(auto &row:inv)for(double &q:row)q/=dg;
        const auto wall=radius(d,p);
        const Mat frame{wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
        const double hs=wall.R*std::sin(p[2]),dv=w*4*pi*determinant(frame)*hs*mj;
        const auto native=local_work(d,{x,t,a},m,mj,frame,hs,dv,chi,dchi,core);
        const bool reconstructed_core=y[0]<=d.x.v[(d.x.n-1)*2/3];
        const auto projected=local_work(d,y,inv,1./dg,frame,hs,dv,chi,dchi,reconstructed_core);
        const auto gradient=local_work(d,{x,t,a},inv,1./dg,frame,hs,dv,chi,dchi,core);
        const double work[]={native.stress,native.order,projected.stress,projected.order,gradient.stress,gradient.order};
        for(int i=0;i<6;i++){if(!std::isfinite(work[i]))return 5;v[i]+=work[i];}
        v[6]+=1;
        return 0;
    }
};
}
extern "C" int nested_depth_work(const Data *d,Results *out,const double *x,int nx,const double *t,int nt,
                                const double *a,int na,const double *xi,const double *tdot,double *report){
    SpatialLabels labels{{nx,x},{nt,t},{na,a},xi};
    long double rows[7]{};
    const int rc=integrate_solid(d,out,OriginalPrisms{},DepthIntegrand{labels,tdot,rows});
    for(int i=0;i<7;i++)report[i]=double(rows[i]);
    return rc;
}
