// Physical-cell integration of one Q1 inverse label chart with two stored natural metrics.
#include "mouth_free_shape_solid.cpp"

namespace {
Mat invert(const Mat &m) {
    Mat q=transpose(cofactor(m));
    double d=determinant(m);
    for (auto &row:q) for(double &v:row) v/=d;
    return q;
}
Mat physical_basis(const Map &r,double length,const Vec &widths) {
    Mat c={r.grad,Vec{0,length,0},Vec{0,0,r.R}};
    for(auto &row:c) for(int i=0;i<3;i++)row[i]*=widths[i];
    return c;
}
struct Landed {
    Vec xi{}; Mat D{}; double chi=0.;
};
struct Q1 {
    Vec labels[2][2][2]; double order[2][2][2];
    Landed at(double u,double v,double w) const {
        Landed p{};
        for(int i=0;i<2;i++) for(int j=0;j<2;j++) for(int k=0;k<2;k++){
            const double x=i?u:1-u, t=j?v:1-v, a=k?w:1-w, h=x*t*a;
            const Vec g={(2*i-1)*t*a,(2*j-1)*x*a,(2*k-1)*x*t};
            p.chi+=h*order[i][j][k];
            for(int c=0;c<3;c++){
                p.xi[c]+=h*labels[i][j][k][c];
                for(int dir=0;dir<3;dir++)p.D[c][dir]+=g[dir]*labels[i][j][k][c];
            }
        }
        return p;
    }
};
struct Geometry {
    Data d; std::array<double,3> lo,width;
};
struct Evaluated {
    Map wall; Mat C; double R,hoop,J,volume,W; Vec coord;
};
Evaluated spatial(const Geometry &g,const Vec &frac,const Mat &right,double det_right,double S,double an,double chi) {
    // right = (d xi / d local_cell)^(-1) (d X_s / d xi)^(-1).
    Vec pos; for(int j=0;j<3;j++)pos[j]=g.lo[j]+g.width[j]*frac[j];
    Map wall=radius(g.d,pos);
    Mat C=physical_basis(wall,g.d.length,g.width);
    const double h=wall.R*std::sin(pos[2])/(S*an), detC=determinant(C);
    const double J=detC*det_right*h;
    const double volume=4*pi*detC*wall.R*std::sin(pos[2]);
    if(!(J>0 && std::isfinite(J)))return {wall,C,wall.R,h,J,volume,0.,pos};
    const Mat F=product(C,right);
    double W=.5*(norm2(F)+h*h-4)-std::log(J)+.25*(J-1)*(J-1);
    return {wall,C,wall.R,h,J,volume,W,pos};
}
} // namespace

// arrays gs contain 1+4*6 independent shape samplings with physical knot locations.
// out [E, six forces, 24 signed energies]; bad [flags, xi, chi, signed inv J,
// physical chart, parent cell, minimum Jinv active, global min chart determinant].
// Angular hoop and the dependent azimuthal coordinate are retained exactly.
extern "C" int cp10_spatial_integral(const Data *ds,int ns,const double *t_deriv,const double *xi,const double *chi,
              const double *nodes,const double *weights,int n,double *out,double *bad,double *checks,double *el) {
    // The nested-grid audit only integrates positive-order cells. The legacy
    // values 1 and 25 additionally survey pure-melt orientation, and 25 also
    // carries 24 signed perturbations. Both use exactly the same integrand.
    const bool active_only=ns==0 || ns==-1;
    const bool derivatives=ns==25 || ns==-1;
    const bool signed_probes=ns==25;
    if(!active_only && ns!=1 && !signed_probes)return 1;
    const Data &d=ds[0]; int nx=d.field_x.n,nt=d.field_t.n,na=d.field_a.n;
    std::vector<long double> result(31,0.), per((d.x.n-1)*(d.t.n-1)*(d.a.n-1)*2,0.);
    bad[15]=bad[16]=std::numeric_limits<double>::infinity();
    auto node=[=](int i,int j,int k){return (i*nt+j)*na+k;};
    const double cut=d.x.v[(d.x.n-1)*2/3];
    for(int i=0;i<nx-1;i++)for(int j=0;j<nt-1;j++)for(int k=0;k<na-1;k++){
        if(active_only){
            bool has_order=false;
            for(int I=0;I<2;I++)for(int T=0;T<2;T++)for(int A=0;A<2;A++)
                if(chi[node(i+I,j+T,k+A)]>0)has_order=true;
            if(!has_order)continue;
        }
        Geometry gs[25];Q1 q;bool melted=true;
        for(int l=0;l<(active_only?1:ns);l++){
            const auto &h=ds[l];
            gs[l]={h,{h.field_x.v[i],h.field_t.v[j],h.field_a.v[k]},
                       {h.field_x.v[i+1]-h.field_x.v[i],h.field_t.v[j+1]-h.field_t.v[j],h.field_a.v[k+1]-h.field_a.v[k]}};
        }
        for(int I=0;I<2;I++)for(int T=0;T<2;T++)for(int A=0;A<2;A++){
            int id=node(i+I,j+T,k+A);
            q.order[I][T][A]=chi[id]; if(chi[id]>0)melted=false;
            for(int c=0;c<3;c++)q.labels[I][T][A][c]=xi[id*3+c];
        }
        const Geometry &g=gs[0];
        for(int T=0;T<n;T++)for(int A=0;A<n;A++){
            double v=nodes[T],w=nodes[A];
            Landed l=q.at(0.,v,w),r=q.at(1.,v,w);
            std::vector<double> cuts={0.,1.};
            if(!melted){
                // At fixed spatial depth and angle each label is affine along this radial line.
                // Split the true weld at xi_r=cut and the material-element coordinate switches.
                const Grid axes[]={d.x,d.t,d.a,d.reg_t};
                const int components[]={0,1,2,1};
                for(int s=0;s<4;s++){
                    int c=components[s];double dx=r.xi[c]-l.xi[c];
                    if(dx==0)continue;
                    for(int h=1;h<axes[s].n-1;h++){
                        double f=(axes[s].v[h]-l.xi[c])/dx;
                        if(f>0 && f<1)cuts.push_back(f);
                    }
                }
                if(r.xi[0]!=l.xi[0])for(int p=0;p<d.reg_e.n;p++){
                    double h=std::exp(-d.pressure*d.reg_e.v[p]);
                    double f=(h-l.xi[0])/(r.xi[0]-l.xi[0]);
                    if(h>=cut && f>0 && f<1)cuts.push_back(f);
                }
                std::sort(cuts.begin(),cuts.end());
#ifdef CP10_SPLIT_REGISTRY_DIAGONALS
                // The derivatives of the collar rest map jump at each triangular-cell
                // diagonal. With moving labels, sampling that jump would contaminate
                // a signed force even though the metric location is known exactly.
                std::vector<double> diagonals;
                const double dr=r.xi[0]-l.xi[0],dt=r.xi[1]-l.xi[1];
                for(size_t z=1;z<cuts.size();z++){
                    const double lower=cuts[z-1],upper=cuts[z],mid=(lower+upper)*.5;
                    const double xm=l.xi[0]+mid*dr,tm=l.xi[1]+mid*dt;
                    if(!(upper>lower && xm>cut))continue;
                    const int ireg=cell(d.reg_e,-std::log(xm)/d.pressure),jreg=cell(d.reg_t,tm);
                    const double de=d.reg_e.v[ireg+1]-d.reg_e.v[ireg],ds=d.reg_t.v[jreg+1]-d.reg_t.v[jreg];
                    auto coordinate=[&](double x){return (-std::log(l.xi[0]+dr*x)/d.pressure-d.reg_e.v[ireg])/de
                                         +(l.xi[1]+dt*x-d.reg_t.v[jreg])/ds-1.;};
                    std::vector<double> sides={lower,upper};
                    if(dr!=0. && dt!=0.){
                        const double turning=(dr*ds/(d.pressure*de*dt)-l.xi[0])/dr;
                        if(turning>lower && turning<upper)sides.insert(sides.begin()+1,turning);
                    }
                    for(size_t k=1;k<sides.size();k++){
                        double a=sides[k-1],b=sides[k],fa=coordinate(a),fb=coordinate(b);
                        if((fa>0.)==(fb>0.))continue;
                        for(int nroot=0;nroot<56;nroot++){
                            double m=(a+b)*.5;
                            if(m==a || m==b)break;
                            if((coordinate(m)>0.)==(fa>0.))a=m;else b=m;
                        }
                        diagonals.push_back((a+b)*.5);
                    }
                }
                cuts.insert(cuts.end(),diagonals.begin(),diagonals.end());
                std::sort(cuts.begin(),cuts.end());
#endif
            }
            for(size_t z=1;z<cuts.size();z++)for(int X=0;X<n;X++){
                double du=cuts[z]-cuts[z-1]; if(du<=0)continue;
                double u=cuts[z-1]+du*nodes[X], wt=du*weights[X]*weights[T]*weights[A];
                Vec frac={u,v,w}; auto p=q.at(u,v,w);
                double detD=determinant(p.D);
                const double dchart=detD/(g.width[0]*g.width[1]*g.width[2]);
                bad[16]=std::min(bad[16],dchart);
                const bool core=p.xi[0]<=cut;
                Ref ref=reference(d,p.xi[0],p.xi[1],profile(d,p.xi[1],true),core);
                const double an=std::sin(p.xi[2]);
                if(!(ref.S>0 && an>0 && ref.det2>0))return 2;
                Vec pos;for(int cc=0;cc<3;cc++)pos[cc]=g.lo[cc]+g.width[cc]*frac[cc];
                auto map=radius(d,pos);
                double physical_den=map.grad[0]*d.length*map.R*map.R*std::sin(pos[2]);
                double invJ=dchart*ref.det2*ref.S*ref.S*an/physical_den;
                if(p.chi>0)bad[15]=std::min(bad[15],invJ);
                if(!(invJ>0)){
                    // Retain the worst active violation (or the worst pure-melt violation while none is active).
                    int severity=p.chi>0?2:1;
                    if(severity>bad[0] || (severity==bad[0] && invJ<bad[5])){
                        bad[0]=severity; for(int c=0;c<3;c++){bad[1+c]=p.xi[c];bad[6+c]=pos[c];}
                        bad[4]=p.chi;bad[5]=invJ;bad[9]=i;bad[10]=j;bad[11]=k;bad[12]=core;
                    }
                    if(p.chi>0){bad[13]++;continue;}
                    bad[14]++;continue;
                }
                if(p.chi==0)continue;
                // Signed control at the first real positive-order point on each side.
                // In a single cell, reflect the first label component about its current value.
                // This affine map of all eight first components leaves the point's labels
                // and material metric fixed and flips exactly one gradient row.
                int side=core?0:1;
                if(checks[4*side]==0){
                    Q1 test=q;
                    for(auto &ii:test.labels)for(auto &jj:ii)for(auto &vertex:jj)vertex[0]=2*p.xi[0]-vertex[0];
                    auto neg=test.at(u,v,w);
                    double neg_inv=invJ*determinant(neg.D)/detD;
                    if(neg_inv>0 || std::abs(neg.xi[0]-p.xi[0])>2e-14)return 3;
                    checks[4*side]=p.chi;checks[4*side+1]=invJ;checks[4*side+2]=neg_inv;
                    checks[4*side+3]=std::abs(neg.xi[0]-p.xi[0]);
                }
                Mat right=product(invert(p.D),ref.inv);
                double detright=1/(detD*ref.S*ref.det2);
                auto b=spatial(g,frac,right,detright,ref.S,an,p.chi);
                if(!(b.J>0))return 4;
                int I=cell(d.x,p.xi[0]),Tref=cell(d.t,p.xi[1]),K=cell(d.a,p.xi[2]);
                int tri=(p.xi[0]-d.x.v[I])/(d.x.v[I+1]-d.x.v[I])+(p.xi[1]-d.t.v[Tref])/(d.t.v[Tref+1]-d.t.v[Tref])>1;
                int ie=((Tref*(d.x.n-1)+I)*(d.a.n-1)+K)*2+tri;
                const double pref=wt*b.volume*p.chi/b.J;
                result[0]+=pref*b.W;per[ie]+=pref*b.W;
                if(!derivatives)continue;
                // Because xi and chi are nodal transported material/scalar coordinates in ALE cell fractions,
                // vref*dchi is stationary. All current-grid terms appear in F and dJ.
                const Mat F=product(b.C,right);
                double tx=b.coord[1],rr=b.coord[0],aa=b.coord[2],prod=tx*(1-tx);
                double tdot=t_deriv[j]+v*(t_deriv[j+1]-t_deriv[j]);
                double dwidth=t_deriv[j+1]-t_deriv[j];
                double bw[]={0.,(d.shape[0]+1.5)*tx,4*prod,4*prod*(2*tx-1),16*prod*prod};
                double bp[]={0.,d.shape[0]+1.5,4*(1-2*tx),4*(-6*tx*tx+6*tx-1),32*prod*(1-2*tx)};
                for(int mode=0;mode<6;mode++){
                    Vec dg{}; double dr=0., depth=mode==0?d.length:0.;
                    if(mode==0){
                        dr=b.wall.grad[1]*tdot;
                        for(int c=0;c<3;c++)dg[c]=b.wall.hess[c][1]*tdot;
                    }else if(mode<5){dr=rr*bw[mode];dg={bw[mode],rr*bp[mode],0.};}
                    else{
                        const double h=(tx-.06)/.56;
                        const double f=h>0 && h<1?16*h*h*(1-h)*(1-h):0.;
                        const double ft=h>0 && h<1?32*h*(1-h)*(1-2*h)/.56:0.;
                        const double c=std::cos(aa),P=(3*c*c-1)/2;
                        dr=rr*rr*rr*f*P;
                        dg={3*rr*rr*f*P,rr*rr*rr*ft*P,-3*rr*rr*rr*f*std::sin(aa)*c};
                    }
                    Mat dm={dg,Vec{0.,depth,0.},Vec{0.,0.,dr}};
                    for(auto &row:dm)for(int c=0;c<3;c++)row[c]*=g.width[c];
                    if(mode==0){dm[0][1]+=b.wall.grad[1]*dwidth;dm[1][1]+=d.length*dwidth;}
                    Mat df=product(dm,right);
                    double work=b.hoop*b.hoop*dr/b.R;
                    for(int c=0;c<3;c++)work+=dot(F[c],df[c]);
                    double div=dg[0]/b.wall.grad[0]+depth/d.length+2*dr/b.R+(mode==0?dwidth/g.width[1]:0.);
                    work+=(-1+.5*b.J*(b.J-1))*div;
                    result[1+mode]-=pref*work;
                    if(!signed_probes)continue;
                    for(int sig=0;sig<4;sig++){
                        auto ss=spatial(gs[1+mode*4+sig],frac,right,detright,ref.S,an,p.chi);
                        if(!(ss.J>0))return 5;
                        result[7+mode*4+sig]+=wt*ss.volume*p.chi/ss.J*ss.W;
                    }
                }
            }
        }
    }
    for(size_t i=0;i<per.size();i++)el[i]=double(per[i]);
    for(int i=0;i<31;i++)out[i]=double(result[i]);
    return 0;
}
