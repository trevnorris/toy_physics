// Shared physical material points of the checkpoint-10 order-work integral.
// Accumulate by native spatial cell and original material band/side. The
// reference metric and where it is evaluated can be held fixed while only
// the inverse-label gradient is replaced, providing the hybrid control.
#include "mouth_cp10_order_depth_work.cpp"

namespace {
constexpr int ncol=10, levels=3, npoint=17;
struct State {double J,W,energy,order;};

State measure(const Data &d,const Vec &xi,const Mat &forward,double dj,const Mat &frame,
              double hoop_radius,double dv,double chi,double minus_dchi,bool core) {
    const auto ref=reference(d,xi[0],xi[1],profile(d,xi[1],true),core);
    const Mat F=product(frame,product(forward,ref.inv));
    const double hoop=hoop_radius/(ref.S*std::sin(xi[2]));
    const double J=determinant(frame)*dj*hoop/(ref.det2*ref.S);
    if(!(J>0. && std::isfinite(J)))return {J,NAN,NAN,NAN};
    const double W=.5*(norm2(F)+hoop*hoop-4.)-std::log(J)+.25*(J-1.)*(J-1.);
    return {J,W,chi*dv*W/J,dv/J*minus_dchi*W};
}
struct Bins {
    const SpatialLabels *grids;
    const double *rate;
    long double *total,*spatial,*band;
    long double **per_grid;
    double *witness;
    long long &count;
    int operator()(const Data &d,Results &,const Profile &,const TriSlice &s,
                   double t,double a,double x,double weight,bool core,double &den,double &denx)const {
        den=denx=0.;
        Vec p=s.X;
        Mat m=s.M;
        for(int c=0;c<3;c++){
            p[c]+=s.M[c][0]*(x-s.x0);
            m[c][2]+=s.Max[c]*(x-s.x0);
        }
        const auto r=basis(d.field_x,p[0]),z=basis(d.field_t,p[1]),b=basis(d.field_a,p[2]);
        double chi=0.,slope=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            const double c=d.correction[field_id(d,r.i+i,b.i+k,z.i+j)];
            chi+=r.w[i]*z.w[j]*b.w[k]*c;
            slope+=r.w[i]*z.v[j]*b.w[k]*c;
        }
        const double minus_dchi=slope*(z.w[0]*rate[z.i]+z.w[1]*rate[z.i+1]);
        if(chi==0. && minus_dchi==0.)return 0;
        if(!(chi>0.))return 1;
        const double mj=determinant(m);
        const auto wall=radius(d,p);
        const Mat frame{wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
        const double hs=wall.R*std::sin(p[2]),dv=weight*4*pi*determinant(frame)*hs*mj;
        if(!(dv>0.))return 2;
        const Vec original{x,t,a};
        const auto old=measure(d,original,m,mj,frame,hs,dv,chi,minus_dchi,core);
        const int home=(r.i*(d.field_t.n-1)+z.i)*(d.field_a.n-1)+b.i;
        const int region=2*cell(d.t,t)+(core?0:1);
        for(int level=0;level<levels;level++){
            const auto &grid=grids[level];
            Vec xi{};Mat grad{};
            grid.interpolate(p,xi,grad);
            const double dg=determinant(grad);
            if(!(dg>0.))return 3;
            Mat inv=transpose(cofactor(grad));
            for(auto &row:inv)for(double &v:row)v/=dg;
            const bool side=xi[0]<=d.x.v[(d.x.n-1)*2/3];
            const auto projected=measure(d,xi,inv,1./dg,frame,hs,dv,chi,minus_dchi,side);
            const auto gradient=measure(d,original,inv,1./dg,frame,hs,dv,chi,minus_dchi,core);
            const double de=projected.energy-old.energy,df=projected.order-old.order;
            const double values[ncol]={old.energy,projected.energy,old.order,projected.order,de,df,
                                       gradient.energy-old.energy,gradient.order-old.order,std::abs(de),std::abs(df)};
            for(int c=0;c<ncol;c++){
                if(!std::isfinite(values[c]))return 4;
                total[level*ncol+c]+=values[c];
                spatial[(home*levels+level)*ncol+c]+=values[c];
                band[(region*levels+level)*ncol+c]+=values[c];
            }
            const int i=cell(grid.x,p[0]),j=cell(grid.t,p[1]),k=cell(grid.a,p[2]);
            const int dest=(i*(grid.t.n-1)+j)*(grid.a.n-1)+k;
            per_grid[level][2*dest]+=de;
            per_grid[level][2*dest+1]+=df;
            if(level==2){
                double *v=witness+npoint*dest;
                if(std::abs(df)>std::abs(v[0])){
                    const double at_point[npoint]={df,wall.R,d.length*p[1],p[2],chi,minus_dchi,
                                                   old.J,projected.J,old.W,projected.W,dv,t,core?0.:1.,
                                                   x,a,p[0],p[1]};
                    std::copy_n(at_point,npoint,v);
                }
            }
        }
        ++count;
        return 0;
    }
};
}

extern "C" int bucket_work(Data *d,Results *out,const Grid *axes,const double **labels,
                            const double *rate,const unsigned char *active,long double *total,
                            long double *spatial,long double *band,long double **per_grid,
                            double *witness,long long *count) {
    SpatialLabels grids[levels];
    for(int i=0;i<levels;i++)grids[i]={axes[3*i],axes[3*i+1],axes[3*i+2],labels[i]};
    ActivePrisms selected{{},active,d->x.n-1,d->a.n-1};
    return integrate_solid(d,out,selected,Bins{grids,rate,total,spatial,band,per_grid,witness,*count});
}
