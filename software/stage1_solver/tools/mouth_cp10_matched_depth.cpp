// Frozen CP10 material map. At fixed native polar angle each physical prism
// slice is an affine triangle. Sweep this triangle at fixed physical depth:
// moving spatial depth faces then delimit finite ribbons on BOTH signed maps.
#include "mouth_cp10_nested_energy.cpp"
#include <stdexcept>
namespace matched {
using List=std::vector<double>;
Vec mix(const Vec &a,const Vec &b,double f){
    Vec v;for(int i=0;i<3;i++)v[i]=a[i]+f*(b[i]-a[i]);return v;
}
void append_knots(List &v,Grid g,double l,double h){
    for(auto p=std::upper_bound(g.v,g.v+g.n,l);p<g.v+g.n&&*p<h;p++)v.push_back(*p);
}
void tidy(List &c){std::sort(c.begin(),c.end());c.erase(std::unique(c.begin(),c.end()),c.end());}
void inner(List &c,double x,double lo=0.,double hi=1.){if(x>lo&&x<hi)c.push_back(x);}
void linear_faces(List &v,Grid g,double l,double h){
    if(l==h)return;
    for(auto p=std::upper_bound(g.v,g.v+g.n,std::min(l,h));p<g.v+g.n&&*p<std::max(l,h);p++)inner(v,(*p-l)/(h-l));
}
struct Poly {
    double a,b,c;
    double operator()(double x)const{return (c*x+b)*x+a;}
    double deriv(double x)const{return 2*c*x+b;}
    static Poly quadratic(double l,double m,double h){double c=2*((h-m)-(m-l));return {l,h-l-c,c};}
    static Poly line(double l,double h){return {l,h-l,0};}
};
void roots(List &out,double A,double B,double C,double lo,double hi){
    // Stable quadratic discriminant/root pair, including near-linear cases.
    if(A==0.){if(B!=0.)inner(out,-C/B,lo,hi);return;}
    double D=B*B-4*A*C;if(D<0)return;
    double q=-.5*(B+std::copysign(std::sqrt(D),B));
    if(q==0.){inner(out,-B/(2*A),lo,hi);return;}
    inner(out,q/A,lo,hi);inner(out,C/q,lo,hi);
}
template<class Fn> void crossings(List &out,const List &extrema,Fn value){
    for(size_t k=1;k<extrema.size();k++){
        double l=extrema[k-1],r=extrema[k],fl=value(l),fr=value(r);
        if(fl==0.){inner(out,l);continue;}
        if(fr==0.){inner(out,r);continue;}
        if(std::signbit(fl)==std::signbit(fr))continue;
        for(int it=0;it<54;it++){double m=(l+r)*.5;if(m==l||m==r)break;
            if(std::signbit(value(m))==std::signbit(fl))l=m;else r=m;}
        inner(out,(l+r)*.5);
    }
}
void poly_faces(List &v,const Poly &p,Grid nodes){
    double l=p(0),r=p(1),mn=std::min(l,r),mx=std::max(l,r);
    if(p.c!=0.){double e=-p.b/(2*p.c);if(e>0&&e<1){mn=std::min(mn,p(e));mx=std::max(mx,p(e));}}
    for(auto z=std::lower_bound(nodes.v,nodes.v+nodes.n,mn);z<nodes.v+nodes.n&&*z<=mx;z++)
        roots(v,p.c,p.b,p.a-*z,0.,1.);
}
void metric_faces(List &all,const Data &d,const Poly &r,const Poly &z,Grid radial,double weld){
    List c{0.,1.};
    roots(c,r.c,r.b,r.a-weld,0.,1.);
    double mx=std::max(r(0),r(1));
    if(r.c<0){double v=-r.b/(2*r.c);if(v>0&&v<1)mx=std::max(mx,r(v));}
    if(mx>=weld){
        poly_faces(c,r,radial);poly_faces(c,z,d.reg_t);tidy(c);
        for(size_t i=1;i<c.size();i++){
            double l=c[i-1],u=c[i],mid=(l+u)*.5;if(r(mid)<=weld)continue;
            const int I=cell(d.reg_e,-std::log(r(mid))/d.pressure),T=cell(d.reg_t,z(mid));
            const double de=d.reg_e.v[I+1]-d.reg_e.v[I],dt=d.reg_t.v[T+1]-d.reg_t.v[T];
            const auto f=[&](double x){return (-std::log(r(x))/d.pressure-d.reg_e.v[I])/de+(z(x)-d.reg_t.v[T])/dt-1;};
            // The derivative numerator g = p*de*r*z' - dt*r' is cubic.
            const double b=d.pressure*de;
            const double g3=b*2*r.c*z.c,g2=b*(r.c*z.b+2*r.b*z.c),
                         g1=b*(r.b*z.b+2*r.a*z.c)-2*dt*r.c,g0=b*r.a*z.b-dt*r.b;
            List crit{l,u};roots(crit,3*g3,2*g2,g1,l,u);tidy(crit);
            List ex{l,u};crossings(ex,crit,[&](double x){return ((g3*x+g2)*x+g1)*x+g0;});tidy(ex);
            crossings(all,ex,f);
        }
    }
    all.insert(all.end(),c.begin(),c.end());
}
struct Node{Vec m,p;};
struct Line{Vec m,p,dy;};
struct Bank {
    SpatialLabels mesh;
    const double *chi;
    double L;
};
struct Fields {
    Vec label[2][2]; Vec depth_der[2][2]; double order[2][2];
    double x0,dx,a0,da;
    bool active=false;
    Fields(const Bank &bank,const Vec &p){
        const auto &g=bank.mesh;
        int i=cell(g.x,p[0]),k=cell(g.a,p[2]); auto t=basis(g.t,p[1]);
        x0=g.x.v[i];dx=g.x.v[i+1]-x0;a0=g.a.v[k];da=g.a.v[k+1]-a0;
        for(int X=0;X<2;X++)for(int A=0;A<2;A++){
            auto n=((i+X)*g.t.n+t.i)*g.a.n+k+A;
            double lm=bank.chi[n],lp=bank.chi[n+g.a.n];
            order[X][A]=t.w[0]*lm+t.w[1]*lp;
            active|=order[X][A]>0.;
            for(int c=0;c<3;c++){
                double l=g.nodes[3*n+c],r=g.nodes[3*(n+g.a.n)+c];
                label[X][A][c]=t.w[0]*l+t.w[1]*r;
                depth_der[X][A][c]=(r-l)*t.v[1];
            }
        }
    }
    Vec value(const Vec &p)const{
        double x=(p[0]-x0)/dx,a=(p[2]-a0)/da;return mix(mix(label[0][0],label[1][0],x),mix(label[0][1],label[1][1],x),a);
    }
    double evaluate(const Vec &p,Vec &xi,Mat &D)const {
        const double x=(p[0]-x0)/dx,a=(p[2]-a0)/da;
        const double wx[]={1-x,x},wa[]={1-a,a}; double chi=0;
        for(int X=0;X<2;X++)for(int A=0;A<2;A++){
            const double w=wx[X]*wa[A]; chi+=w*order[X][A];
            for(int c=0;c<3;c++){
                xi[c]+=w*label[X][A][c];D[c][0]+=(2*X-1)*wa[A]/dx*label[X][A][c];
                D[c][1]+=w*depth_der[X][A][c];D[c][2]+=(2*A-1)*wx[X]/da*label[X][A][c];
            }
        }
        return chi;
    }
};
struct Rows { long double E[2][2]{},force[2]{},vol=0,count=0,active=0;
    double min_j=std::numeric_limits<double>::infinity();};
struct Audit{
    const Data &d; Grid gauss,radial;const double *weights;
    const Bank banks[2];double h,weld;Rows rows;
    void section(int i,int j,int k,int half,const Line &lo,const Line &hi,double weight){
        const Vec &l=lo.p,&r=hi.p; const auto &g=banks[0].mesh;
        const double det=std::abs((hi.m[0]-lo.m[0])*lo.dy[1]-(hi.m[1]-lo.m[1])*lo.dy[0]);
        if(!det)return;
        rows.vol+=det*weight;
        List cuts{0.,1.};linear_faces(cuts,g.x,l[0],r[0]);linear_faces(cuts,g.a,l[2],r[2]);tidy(cuts);
        const bool native_core=i<(d.x.n-1)*2/3;
        for(size_t idx=1;idx<cuts.size();idx++){
            const double bl=cuts[idx-1],br=cuts[idx]; if(br<=bl)continue;
            const Vec PL=mix(l,r,bl),PR=mix(l,r,br),PM=mix(PL,PR,.5),ML=mix(lo.m,hi.m,bl),MR=mix(lo.m,hi.m,br);
            Fields f[2]={Fields(banks[0],PM),Fields(banks[1],PM)};
            if(!f[0].active&&!f[1].active)continue;
            List parts{0.,1.};
            if(!native_core)metric_faces(parts,d,Poly::line(ML[0],MR[0]),Poly::line(ML[1],MR[1]),radial,weld);
            for(int b=0;b<2;b++)if(f[b].active){
                auto yl=f[b].value(PL),ym=f[b].value(PM),yr=f[b].value(PR);
                metric_faces(parts,d,Poly::quadratic(yl[0],ym[0],yr[0]),Poly::quadratic(yl[1],ym[1],yr[1]),radial,weld);
            }
            tidy(parts);
            for(size_t v=1;v<parts.size();v++)for(int z=0;z<gauss.n;z++){
                double dl=parts[v]-parts[v-1];if(dl<=0)continue;
                double s=parts[v-1]+dl*gauss.v[z];
                const Vec p=mix(PL,PR,s),m=mix(ML,MR,s);
                const double w=weight*det*(br-bl)*dl*weights[z];
                // Use an interior seed strictly on this triangle. The affine
                // element derivative has a unique extension to its closure.
                double vj=(m[1]-d.t.v[j])/(d.t.v[j+1]-d.t.v[j]);
                double xm=d.x.v[i]+(d.x.v[i+1]-d.x.v[i])*(half?1-.5*vj:.5*(1-vj));
                TriSlice elem(d,i,j,k,m[1],m[2],xm);
                for(int c=0;c<3;c++)elem.M[c][2]+=elem.Max[c]*(m[0]-xm);
                const double jm=determinant(elem.M);
                auto wall=radius(d,p); double hs=wall.R*std::sin(p[2]);
                double energies[2][2]{};
                for(int b=0;b<2;b++)if(f[b].active){
                    Vec y{};Mat G{};double chi=f[b].evaluate(p,y,G);if(chi<=0.)continue;
                    double dg=determinant(G);if(!(dg>0&&jm>0))throw std::runtime_error("ordered map orientation");
                    Mat inv=transpose(cofactor(G));for(auto &row:inv)for(auto &c:row)c/=dg;
                    Mat C={wall.grad,Vec{0.,banks[b].L,0.},Vec{0.,0.,wall.R}};
                    double dv=4*pi*determinant(C)*hs*jm;
                    auto orig=constitutive(d,m,elem.M,C,jm,hs,dv,chi,native_core);
                    auto proj=constitutive(d,y,inv,C,1/dg,hs,dv,chi,y[0]<=weld);
                    if(!std::isfinite(orig.weighted+proj.weighted))throw std::runtime_error("nonfinite strain");
                    rows.min_j=std::min({rows.min_j,1/orig.J,1/proj.J});
                    energies[b][0]=orig.weighted;energies[b][1]=proj.weighted;rows.active++;
                }
                for(int s=0;s<2;s++){
                    for(int b=0;b<2;b++)rows.E[b][s]+=w*energies[b][s];
                    rows.force[s]+=w*((energies[0][s]-energies[1][s])/(2*h));
                }rows.count++;
            }
        }
    }
    void prism(int i,int j,int k,int half,double a,double aw,double highest){
        Node vertices[3];
        for(int v=0;v<3;v++){
            const int ix=i+(half?(v!=1):(v==1)), it=j+(half?(v!=2):(v==2));
            const double y=(a-d.a.v[k])/(d.a.v[k+1]-d.a.v[k]);
            auto &m=vertices[v].m;auto &p=vertices[v].p;
            m={d.x.v[ix],d.t.v[it],a}; p=m; const int off=((it*d.x.n+ix)*d.a.n+k)*3;
            for(int c=0;c<3;c++)p[c]+=(1-y)*d.offset[off+c]+y*d.offset[off+3+c];
        }
        List zz{vertices[0].p[1],vertices[1].p[1],vertices[2].p[1]};tidy(zz);
        const double l=zz.front(),u=std::min(zz.back(),highest);if(u<=l)return;
        for(const Bank &b:banks)append_knots(zz,b.mesh.t,l,u);
        zz.push_back(u);tidy(zz);
        for(size_t q=1;q<zz.size();q++){
            const double dz=zz[q]-zz[q-1];if(zz[q]>u||dz<=0)continue;
            double zm=(zz[q]+zz[q-1])*.5;
            std::array<std::pair<int,int>,3> sides;int count=0;
            for(int edge=0;edge<3;edge++){
                int e=(edge+1)%3;double za=vertices[edge].p[1],zb=vertices[e].p[1];
                if(std::min(za,zb)<zm&&zm<std::max(za,zb))sides[count++]={edge,e};
            }
            if(count!=2){
                if(dz<3e-15)continue;
                throw std::runtime_error("triangle-depth intersection");
            }
            for(int T=0;T<gauss.n;T++){
                double z=zz[q-1]+dz*gauss.v[T];Line end[2];
                for(int o=0;o<2;o++){
                    auto [l,r]=sides[o];auto &v=vertices[l];auto &y=vertices[r];double dz=y.p[1]-v.p[1],s=(z-v.p[1])/dz;
                    end[o].p=mix(v.p,y.p,s);end[o].p[1]=z;end[o].m=mix(v.m,y.m,s);
                    for(int c=0;c<3;c++)end[o].dy[c]=(y.m[c]-v.m[c])/dz;
                }
                section(i,j,k,half,end[0],end[1],aw*dz*weights[T]);
            }
        }
    }
};
}
extern "C" int matched_depth(const Data *d,const double *x,int nx,const double *a,int na,int nt,
                 const double *tm,const double *tp,const double *ym,const double *yp,const double *chi,
                 double lm,double lp,double h,const double *points,const double *w,int n,int first,int past,double *out){
    using namespace matched;
    List radial;for(int i=0;i<d->reg_e.n;i++)radial.push_back(std::exp(-d->pressure*d->reg_e.v[i]));tidy(radial);
    double cut=d->x.v[(d->x.n-1)*2/3];
    Audit job{*d,Grid{n,points},Grid{int(radial.size()),radial.data()},w,
              {Bank{SpatialLabels{{nx,x},{nt,tm},{na,a},ym},chi,lm},Bank{SpatialLabels{{nx,x},{nt,tp},{na,a},yp},chi,lp}},h,cut,{}};
    int final=0;
    for(int i=0;i<nx;i++)for(int j=0;j<nt;j++)for(int k=0;k<na;k++)if(chi[(i*nt+j)*na+k]>0.)final=std::max(final,j);
    const int active_end=std::min(final+1,nt-1);
    double high=std::max(tm[active_end],tp[active_end]);
    try{
        for(int k=0;k<d->a.n-1;k++){
            if(k<first||k>=past)continue;
            List angles{d->a.v[k],d->a.v[k+1]};append_knots(angles,Grid{na,a},angles.front(),angles.back());tidy(angles);
            for(size_t u=1;u<angles.size();u++)for(int A=0;A<n;A++){
                double da=angles[u]-angles[u-1],alpha=angles[u-1]+da*points[A];
                for(int j=0;j<d->t.n-1;j++)for(int i=0;i<d->x.n-1;i++)for(int v=0;v<2;v++)
                    job.prism(i,j,k,v,alpha,da*w[A],high);
            }
        }
    }catch(std::exception &){return 1;}
    auto &r=job.rows;int c=0;for(auto &e:r.E)for(auto v:e)out[c++]=v;
    for(auto v:r.force)out[c++]=v;out[c++]=r.vol;out[c++]=r.count;out[c++]=r.active;out[c++]=r.min_j;
    return 0;
}
