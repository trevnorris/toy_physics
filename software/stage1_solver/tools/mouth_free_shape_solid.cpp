// Reference-volume weak work of a bonded four-dimensional axisymmetric solid.
// Orthogonal physical rows are (radial, depth, polar arc, azimuthal arc).
// Logical columns are (radius fraction, depth fraction, polar angle, azimuth).
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using Vec = std::array<double,3>;
using Mat = std::array<Vec,3>;
constexpr double pi=3.14159265358979323846;
struct Grid { int n; const double *v; };
struct Data {
    Grid x,t,a,depth_edges,angle_edges,base_x,base_t,field_x,field_t,field_a,cut_x,cut_t,cut_a,reg_e,reg_t,gauss;
    const double *weights,*labels,*base,*correction,*offset,*spline,*shape;
    int spline_n, diagnose;
    double length, ref_length, pressure, core_ratio, epsilon;
};
struct Results { double *gradient,*load,*stats; };
int cell(Grid g,double p) {
    return std::clamp(int(std::upper_bound(g.v,g.v+g.n,p)-g.v)-1,0,g.n-2);
}
double dot(const Vec &v,const Vec &w) {return v[0]*w[0]+v[1]*w[1]+v[2]*w[2];}
Vec cross(const Vec &v,const Vec &w) {return {v[1]*w[2]-v[2]*w[1],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0]};}
Mat transpose(const Mat &b) {Mat c;for(int i=0;i<3;i++)for(int j=0;j<3;j++)c[i][j]=b[j][i];return c;}
Mat product(const Mat &a,const Mat &b) {
    Mat bt=transpose(b),c;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)c[i][j]=dot(a[i],bt[j]);
    return c;
}
double determinant(const Mat &a) {return dot(a[0],cross(a[1],a[2]));}
Mat cofactor(const Mat &a) {return {cross(a[1],a[2]),cross(a[2],a[0]),cross(a[0],a[1])};}
double norm2(const Mat &a) {return dot(a[0],a[0])+dot(a[1],a[1])+dot(a[2],a[2]);}
struct Profile {double b,bt,btt;};
Profile profile(const Data &d,double t,bool ref=false) {
    // scipy CubicHermiteSpline coefficients are in normalized t on uniform knots.
    const int n=d.spline_n, i=std::clamp(int(t*n),0,n-1);
    const double *c=d.spline+4*i;
    double h=t-double(i)/n;
    const double *v=d.shape+(ref?4:0);
    const double e=v[0], m=v[1], a=v[2], extra=v[3], prod=t*(1-t);
    return {4*((((c[0]*h+c[1])*h+c[2])*h+c[3])+e*t+4*prod*(m+a*(2*t-1))+16*extra*prod*prod),
            4*((3*c[0]*h+2*c[1])*h+c[2]+e+4*m*(1-2*t)+4*a*(-6*t*t+6*t-1)+32*extra*prod*(1-2*t)),
            4*(6*c[0]*h+2*c[1]-8*m+4*a*(6-12*t)+32*extra*(1-6*prod))};
}
struct Map {double R; Vec grad; Mat hess;};
Map radius(const Data &d,const Vec &X) {
    double x=X[0],t=X[1],a=X[2];
    auto b=profile(d,t);
    double y=(t-.06)/.56, f=0,ft=0,ftt=0;
    if(y>0 && y<1) {
        f=16*y*y*(1-y)*(1-y)*d.epsilon;
        ft=32*y*(1-y)*(1-2*y)/.56*d.epsilon;
        ftt=32*(1-6*y+6*y*y)/(.56*.56)*d.epsilon;
    }
    double co=std::cos(a),sn=std::sin(a), p=(3*co*co-1)/2, pa=-3*sn*co, paa=3*(sn*sn-co*co);
    double x2=x*x,x3=x2*x;
    Map r{};r.R=x*b.b+x3*f*p;
    r.grad={b.b+3*x2*f*p,x*b.bt+x3*ft*p,x3*f*pa};
    r.hess={Vec{6*x*f*p,b.bt+3*x2*ft*p,3*x2*f*pa},
            Vec{b.bt+3*x2*ft*p,x*b.btt+x3*ftt*p,x3*ft*pa},
            Vec{3*x2*f*pa,x3*ft*pa,x3*f*paa}};
    return r;
}
struct Ref { double S,Sx,det2,logVx; Mat inv,H; };
Ref reference(const Data &d,double x,double t,const Profile &initial,bool core) {
    double sx,st,tx,tt,S,sxx,stx,txx,ttx;
    if(core) {
        sx=initial.b*d.core_ratio;st=x*initial.bt*d.core_ratio;
        S=x*sx;tx=0;tt=d.ref_length;sxx=txx=ttx=0;stx=initial.bt*d.core_ratio;
    } else {
        double e=-std::log(x)/d.pressure;
        int i=cell(d.reg_e,e), j=cell(d.reg_t,t), ne=d.reg_e.n;
        double de=d.reg_e.v[i+1]-d.reg_e.v[i],dt=d.reg_t.v[j+1]-d.reg_t.v[j];
        double u=(e-d.reg_e.v[i])/de,v=(t-d.reg_t.v[j])/dt;
        bool up=u+v>1.;int base=j*ne+i;
        const double *p=d.labels+2*(base+(up?ne+1:0)),
                     *r=d.labels+2*(base+(up?ne:1)),*z=d.labels+2*(base+(up?1:ne));
        double sgn=up?-1.:1.;if(up){u=1-u;v=1-v;}
        S=4+p[0]+u*(r[0]-p[0])+v*(z[0]-p[0]);
        sx=(r[0]-p[0])*sgn/(-d.pressure*x*de);tx=(r[1]-p[1])*sgn/(-d.pressure*x*de);
        st=(z[0]-p[0])*sgn/dt;tt=(z[1]-p[1])*sgn/dt;
        sxx=-sx/x;txx=-tx/x;stx=ttx=0;
    }
    Ref out{};out.S=S;out.Sx=sx;
    double det=sx*tt-st*tx;
    out.det2=std::abs(det);
    out.inv={Vec{tt/det,-st/det,0},Vec{-tx/det,sx/det,0},Vec{0,0,1/S}};
    // H = (d Gref/dx) inv(Gref). It acts on the right of F.
    out.H=product(Mat{Vec{sxx,stx,0},Vec{txx,ttx,0},Vec{0,0,sx}},out.inv);
    out.logVx=out.H[0][0]+out.H[1][1]+2*sx/S;
    return out;
}
struct TriSlice {
    int id[6]; double n[6],nx[6],nt[6],na[6],nax[6],x0;
    Vec X; Mat M; Vec Max{};
    TriSlice(const Data &d,int i,int j,int k,double t,double a,double mid) {
        double dx=d.x.v[i+1]-d.x.v[i],dt=d.t.v[j+1]-d.t.v[j],da=d.a.v[k+1]-d.a.v[k];
        double u=(mid-d.x.v[i])/dx,v=(t-d.t.v[j])/dt,y=(a-d.a.v[k])/da;
        bool up=u+v>1;int base=j*d.x.n+i;int ids[]={base+(up?d.x.n+1:0),base+(up?d.x.n:1),base+(up?1:d.x.n)};
        double sign=up?-1.:1.;if(up){u=1-u;v=1-v;}
        double bary[]={1-u-v,u,v},bx[]={-sign/dx,sign/dx,0},bz[]={-sign/dt,0,sign/dt};
        x0=mid;X={mid,t,a};M={Vec{1,0,0},Vec{0,1,0},Vec{0,0,1}};
        for(int p=0;p<3;p++)for(int b=0;b<2;b++){
            int l=2*p+b;double w=b?y:1-y,dw=b?1/da:-1/da;
            id[l]=ids[p]*d.a.n+k+b;
            n[l]=bary[p]*w;nx[l]=bx[p]*w;nt[l]=bz[p]*w;na[l]=bary[p]*dw;nax[l]=bx[p]*dw;
            const double *off=d.offset+3*id[l];
            for(int c=0;c<3;c++) {X[c]+=n[l]*off[c];M[c][0]+=nx[l]*off[c];
                M[c][1]+=nt[l]*off[c];M[c][2]+=na[l]*off[c];Max[c]+=nax[l]*off[c];}
        }
    }
};
struct Edge {double x;int component;};
void cut(const TriSlice &s,double l,double r,Grid grid,int c,std::vector<Edge> &edges) {
    double slope=s.M[c][0]; if(std::abs(slope)<1e-13)return;
    double v0=s.X[c]+slope*(l-s.x0),v1=s.X[c]+slope*(r-s.x0);
    auto it=std::upper_bound(grid.v+1,grid.v+grid.n-1,std::min(v0,v1));
    for(;it!=grid.v+grid.n-1 && *it<std::max(v0,v1);++it){
        double x=s.x0+(*it-s.X[c])/slope;
        if(x>l+1e-14 && x<r-1e-14)edges.push_back({x,c});
    }
}
struct Basis {int i;double w[2],v[2];};
Basis basis(Grid g,double x) {
    int i=cell(g,x);double h=g.v[i+1]-g.v[i],u=(x-g.v[i])/h;
    return {i,{1-u,u},{-1/h,1/h}};
}
void extremes(const double val,int slot,bool lower,double *stats,const Vec &X){
    if(lower?val<stats[slot]:val>stats[slot]) {stats[slot]=val;for(int i=0;i<3;i++)stats[slot+1+i]=X[i];}
}
void stretches(const Mat &F,double h,double *stats,const Vec &X){
    Mat B=product(F,transpose(F)); double tr=(B[0][0]+B[1][1]+B[2][2])/3;
    for(int i=0;i<3;i++)B[i][i]-=tr;
    double p=std::sqrt(norm2(B)/6),lo=tr,hi=tr;
    if(p>1e-14){
        double angle=std::acos(std::clamp(determinant(B)/(2*p*p*p),-1.,1.))/3;
        hi+=2*p*std::cos(angle); lo+=2*p*std::cos(angle+2*pi/3);
    }
    extremes(std::min(h,std::sqrt(std::max(0.,lo))),17,true,stats,X);
    extremes(std::max(h,std::sqrt(std::max(0.,hi))),21,false,stats,X);
}
// Derivatives include reference integration-point and weight motion at each
// crossing of an Eulerian order knot, so the residual differentiates precisely
// the finite quadrature used by the minimizer.
int sample(const Data &d,Results &out,const Profile &initial,const TriSlice &s,double t,double a,
           double x,double sample_weight,bool core,double &density,double &density_x) {
    double dx=x-s.x0, sn0=std::sin(a);Vec X=s.X;Mat M=s.M;
    for(int i=0;i<3;i++){X[i]+=M[i][0]*dx; M[i][2]+=s.Max[i]*dx;}
    if(X[0]<=0||X[0]>=1||X[1]<=0||X[1]>=1||X[2]<=0||X[2]>=pi/2)return 1;
    Ref ref=reference(d,x,t,initial,core);
    Map wall=radius(d,X);
    Mat C={wall.grad,Vec{0,d.length,0},Vec{0,0,wall.R}};
    Mat B=product(M,ref.inv),F=product(C,B);
    double sn=std::sin(X[2]),h=wall.R*sn/(ref.S*sn0);
    double J=wall.grad[0]*d.length*wall.R*determinant(M)*h/(ref.det2*ref.S);
    if(!(J>0 && wall.R>0 && wall.grad[0]>0 && ref.S>0))return 2;
    double press=-1+.5*J*(J-1),W=.5*(norm2(F)+h*h-4)-std::log(J)+.25*(J-1)*(J-1);
    const Basis bx=basis(d.base_x,X[0]),bt=basis(d.base_t,X[1]);
    double chi=0;Vec dc{};
    for(int j=0;j<2;j++)for(int i=0;i<2;i++){
        double v=d.base[(bt.i+j)*d.base_x.n+bx.i+i];
        chi+=bx.w[i]*bt.w[j]*v;dc[0]+=bx.v[i]*bt.w[j]*v;dc[1]+=bx.w[i]*bt.v[j]*v;
    }
    const Basis fx=basis(d.field_x,X[0]),ft=basis(d.field_t,X[1]),fa=basis(d.field_a,X[2]);
    double V=4*pi*ref.det2*ref.S*ref.S*sn0, wt=V*sample_weight, delta=0.;
    for(int i=0;i<2;i++)for(int k=0;k<2;k++)for(int j=0;j<2;j++){
        int id=((fx.i+i)*d.field_a.n+fa.i+k)*d.field_t.n+ft.i+j;
        double w=fx.w[i]*fa.w[k]*ft.w[j],v=d.correction[id];
        delta+=v*w;dc[0]+=v*fx.v[i]*fa.w[k]*ft.w[j];dc[1]+=v*fx.w[i]*fa.w[k]*ft.v[j];
        dc[2]+=v*fx.w[i]*fa.v[k]*ft.w[j];
        if(out.load)out.load[id]+=wt*W*w;
    }
    chi+=delta; density=V*W*chi;
    out.stats[0]+=sample_weight*density;
    out.stats[1]+=wt*W*delta;out.stats[2]+=wt*W;out.stats[3]+=wt;
    ++out.stats[4];
    if(d.diagnose){
        extremes(J,5,true,out.stats,X);extremes(J,9,false,out.stats,X);
        extremes(chi,13,true,out.stats,X);stretches(F,h,out.stats,X);
        extremes(chi,25,false,out.stats,X);
        // Positive extra energy for trial twist omega=sin(pi*t) x(1-x) cos(a).
        Vec twist={(1-2*x)*std::sin(pi*t)*std::cos(a),pi*x*(1-x)*std::cos(pi*t)*std::cos(a),
                   -x*(1-x)*std::sin(pi*t)*sn0},gw{};
        for(int j=0;j<3;j++)for(int i=0;i<3;i++)gw[j]+=twist[i]*ref.inv[i][j];
        out.stats[29]+=.5*wt*chi*wall.R*wall.R*sn*sn*dot(gw,gw);
    }
    if(!out.gradient)return 0;
    // P = F + press F^{-T}. The signed determinant of F may be negative on
    // the collar because its original inverse labels use inward-increasing e.
    Mat P=cofactor(F);double p=press/determinant(F);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)P[i][j]=F[i][j]+p*P[i][j];
    Mat gradC=product(P,transpose(B));
    Mat gradM=product(product(transpose(C),P),transpose(ref.inv));
    Vec pos{}; double hoop=h*h+press;
    for(int c=0;c<3;c++){
        double derivative=dot(gradC[0],wall.hess[c])+gradC[2][2]*wall.grad[c]
                            +hoop*(wall.grad[c]/wall.R+(c==2?std::cos(X[2])/sn:0));
        pos[c]=W*dc[c]+chi*derivative;
    }
    for(int n=0;n<6;n++)for(int c=0;c<3;c++)out.gradient[3*s.id[n]+c]+=wt*(
        pos[c]*(s.n[n]+dx*s.nx[n])+
        chi*(gradM[c][0]*s.nx[n]+gradM[c][1]*s.nt[n]+gradM[c][2]*(s.na[n]+dx*s.nax[n])));
    // d/dx at fixed material t and angle, within one smooth integration piece.
    Mat FtP=product(transpose(F),P);
    double refpart=0.;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)refpart+=FtP[i][j]*ref.H[i][j];
    double eprime=chi*(W*ref.logVx-refpart-hoop*ref.Sx/ref.S);
    for(int c=0;c<3;c++) eprime+=pos[c]*M[c][0]+chi*gradM[c][2]*s.Max[c];
    density_x=V*eprime;
    return 0;
}
extern "C" int solid_eval(const Data *data,Results *res) {
    const Data &d=*data;Results &out=*res;
    std::fill(out.stats,out.stats+30,0.);
    for(int slot:{5,13,17})out.stats[slot]=std::numeric_limits<double>::infinity();
    for(int slot:{9,21,25})out.stats[slot]=-std::numeric_limits<double>::infinity();
    std::vector<Edge> edges;std::vector<double> fixed;
    int cutnode=(d.x.n-1)*2/3;
    for(int ti=0;ti<d.depth_edges.n-1;ti++)for(int tz=0;tz<d.gauss.n;tz++){
        double dt=d.depth_edges.v[ti+1]-d.depth_edges.v[ti],t=d.depth_edges.v[ti]+dt*d.gauss.v[tz];
        int j=cell(d.t,t);double v=(t-d.t.v[j])/(d.t.v[j+1]-d.t.v[j]);
        int jreg=cell(d.reg_t,t); double vr=(t-d.reg_t.v[jreg])/(d.reg_t.v[jreg+1]-d.reg_t.v[jreg]);
        auto initial=profile(d,t,true); fixed.clear();
        for(int l=0;l<d.reg_e.n-1;l++){
            if(d.reg_e.v[l]<-1e-15)continue;
            fixed.push_back(std::exp(-d.pressure*d.reg_e.v[l]));
            fixed.push_back(std::exp(-d.pressure*(d.reg_e.v[l]+(1-vr)*(d.reg_e.v[l+1]-d.reg_e.v[l]))));
        }
        std::sort(fixed.begin(),fixed.end());
        for(int ai=0;ai<d.angle_edges.n-1;ai++)for(int az=0;az<d.gauss.n;az++){
            double da=d.angle_edges.v[ai+1]-d.angle_edges.v[ai],a=d.angle_edges.v[ai]+da*d.gauss.v[az];
            int k=cell(d.a,a);double weight=dt*d.weights[tz]*da*d.weights[az];
            for(int i=0;i<d.x.n-1;i++){
                double lo=d.x.v[i],hi=d.x.v[i+1],diagonal=lo+(1-v)*(hi-lo);
                for(int half=0;half<2;half++){
                    double l=half?diagonal:lo,r=half?hi:diagonal;
                    if(r-l<1e-15)continue;
                    TriSlice s(d,i,j,k,t,a,(l+r)/2);
                    edges.clear();edges.push_back({l,-1});edges.push_back({r,-1});
                    if(i>=cutnode)for(auto it=std::upper_bound(fixed.begin(),fixed.end(),l);it!=fixed.end() && *it<r;it++)
                        edges.push_back({*it,-1});
                    cut(s,l,r,d.base_x,0,edges);cut(s,l,r,d.base_t,1,edges);
                    cut(s,l,r,d.cut_x,0,edges);cut(s,l,r,d.cut_t,1,edges);cut(s,l,r,d.cut_a,2,edges);
                    std::sort(edges.begin(),edges.end(),[](Edge p,Edge q){return p.x<q.x;});
                    for(size_t z=0;z+1<edges.size();z++) {
                        double xl=edges[z].x,xr=edges[z+1].x,delta=xr-xl,workl=0,workr=0;
                        if(delta<1e-15)continue;
                        for(int q=0;q<d.gauss.n;q++){
                            double g=d.gauss.v[q],w=d.weights[q]*weight,den=0,denx=0;
                            int err=sample(d,out,initial,s,t,a,xl+g*delta,w*delta,i<cutnode,den,denx);
                            if(err)return err;
                            workl+=w*((1-g)*delta*denx-den);workr+=w*(g*delta*denx+den);
                        }
                        if(out.gradient)for(int e=0;e<2;e++) {
                            const Edge &cutedge=edges[z+e];int c=cutedge.component;if(c<0)continue;
                            double adj=-(e?workr:workl)/s.M[c][0];
                            for(int n=0;n<6;n++)out.gradient[3*s.id[n]+c]+=adj*(s.n[n]+(cutedge.x-s.x0)*s.nx[n]);
                        }
                    }
                }
            }
        }
    }
    return 0;
}
