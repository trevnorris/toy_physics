// Crystal shape-work control on each smooth material chart, using the original
// quadrature and the Eulerian Cauchy stress. This does not glue the charts.
#include "mouth_r4_refined_solid.cpp"

struct SixChartWorks {
    long double *out;
    double h;
    int operator()(const Data &d,Results &,const Profile &initial,const TriSlice &s,
                   double t,double a,double x,double weight,bool core,double &den,double &denx) const {
        den=denx=0.;
        Vec X=s.X; Mat M=s.M;
        for(int i=0;i<3;++i) { X[i]+=M[i][0]*(x-s.x0); M[i][2]+=s.Max[i]*(x-s.x0); }
        const Basis bx=basis(d.field_x,X[0]),bt=basis(d.field_t,X[1]),ba=basis(d.field_a,X[2]);
        double chi=0.;
        for(int i=0;i<2;++i)for(int k=0;k<2;++k)for(int j=0;j<2;++j)
            chi+=bx.w[i]*ba.w[k]*bt.w[j]*d.correction[field_id(d,bx.i+i,ba.i+k,bt.i+j)];
        if(chi==0.) return 0;
        if(chi<0.) return 3;
        const Ref ref=reference(d,x,t,initial,core);
        const Map m=radius(d,X);
        const double z=X[1],r=X[0],th=X[2];
        const double sn=std::sin(th),sn0=std::sin(a),hoop=m.R*sn/(ref.S*sn0);
        const Mat C={m.grad,Vec{0,d.length,0},Vec{0,0,m.R}};
        const Mat B=product(M,ref.inv),F=product(C,B);
        const double J=m.grad[0]*d.length*m.R*determinant(M)*hoop/(ref.det2*ref.S);
        if(!(J>0.)) return 2;
        const double jinv=1./J, wt=4*pi*ref.det2*ref.S*ref.S*sn0*weight;
        const double phys=wt*J, pressure=-1.+.5*J*(J-1.);
        const auto energy=[](const Mat &f,double ring,double vol){
            return .5*(norm2(f)+ring*ring-4.)-std::log(vol)+.25*(vol-1.)*(vol-1.);
        };
        out[0]+=phys*jinv*chi*energy(F,hoop,J);
        const double T=z*(1-z),Tz=1-2*z;
        const double b[5]={0.,4*(d.shape[0]+1.5)*.25*z,4*T,4*T*(2*z-1),16*T*T};
        const double bz[5]={0.,4*(d.shape[0]+1.5)*.25,4*Tz,4*(Tz*(2*z-1)+2*T),32*T*Tz};
        for(int n=0;n<6;++n){
            double dr=0.; Vec grad{};
            if(n==5) {
                double y=(z-.06)/.56;
                const double f=(y>0. && y<1.)?16*y*y*(1-y)*(1-y):0.;
                const double fz=(y>0. && y<1.)?32*y*(1-y)*(1-2*y)/.56:0.;
                const double co=std::cos(th),P=(3*co*co-1)/2.;
                dr=r*r*r*f*P; grad={3*r*r*f*P,r*r*r*fz*P,-3*r*r*r*f*co*sn};
            } else {
                dr=r*b[n]; grad={b[n],r*bz[n],0.};
            }
            const double depth_change=n==0?d.length:0.;
            Mat dC={grad,Vec{0,depth_change,0},Vec{0,0,dr}};
            Mat dF=product(dC,B);
            const double div=grad[0]/m.grad[0]+depth_change/d.length+2*dr/m.R;
            const double dh=hoop*dr/m.R;
            // J sigma : grad(v) = F:dF + hoop*d(hoop) + pressure*div(v).
            double work=hoop*dh+pressure*div;
            for(int k=0;k<3;++k)work+=dot(F[k],dF[k]);
            out[1+n]+=phys*jinv*chi*work;
            for(int probe=0;probe<4;++probe){
                const double delta=(probe<2?h:h/2)*(probe%2?1:-1);
                const double eta=n==1?std::expm1(.25*delta)/.25:(n==0?std::expm1(delta):delta);
                Mat cf=C;
                for(int i=0;i<3;++i)for(int j=0;j<3;++j)cf[i][j]+=eta*dC[i][j];
                const double ratio=(cf[0][0]/m.grad[0])*(cf[1][1]/d.length)*std::pow(cf[2][2]/m.R,2);
                const double jac=J*ratio;
                if(!(jac>0.))return 2;
                const double euler_volume=phys*ratio;
                const double euler_density=chi/jac*energy(product(cf,B),hoop+eta*dh,jac);
                out[7+4*n+probe]+=euler_volume*euler_density;
            }
        }
        return 0;
    }
};

extern "C" int inverse_chart_work(const Data *data,const double *surplus,double seam,
                                   double h,double *output){
    long double accumulation[31]={};
    double stats[30]={}; Results results{};
    results.stats=stats;
    AngularEdgeStar star{surplus,seam,true,0,-1,-1,-1,-1};
    int rc=integrate_solid(data,&results,star,SixChartWorks{accumulation,h});
    for(int n=0;n<31;++n)output[n]=static_cast<double>(accumulation[n]);
    return rc;
}
