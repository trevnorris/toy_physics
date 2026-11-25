// Physical-volume crystal energy with one welded Q1 label field, two unchanged
// reference embeddings and exact variations of the quadrature's radial cuts.
#include "mouth_cp10_spatial_crystal.cpp"
#include <stdexcept>

namespace euler {
struct Cut {double u; Vec normal;};  // boundary in material-label space
using Cuts=std::vector<Cut>;
constexpr Vec zero{};

void add_face(Cuts& cuts,double left,double right,double target,int component){
    if(left==right)return;
    const double u=(target-left)/(right-left);
    if(u>0&&u<1){Vec normal{};normal[component]=1;cuts.push_back({u,normal});}
}
void grid_faces(Cuts& cuts,Grid axis,double left,double right,int c){
    if(left==right)return;
    for(auto s=std::upper_bound(axis.v,axis.v+axis.n,std::min(left,right));
             s<axis.v+axis.n&&*s<std::max(left,right);s++)add_face(cuts,left,right,*s,c);
}
void sorted(Cuts &cuts){
    std::sort(cuts.begin(),cuts.end(),[](const Cut&a,const Cut&b){return a.u<b.u;});
    cuts.erase(std::unique(cuts.begin(),cuts.end(),[](const Cut&a,const Cut&b){return a.u==b.u;}),cuts.end());
}
void radial_cuts(Cuts& cuts,const Data&d,const Vec &l,const Vec&r,double weld){
    cuts.clear();cuts.push_back({0.,zero});cuts.push_back({1.,zero});
    if(d.collar_coeff){
        add_face(cuts,l[0],r[0],weld,0);
        if(cuts.size()>2)std::swap(cuts[1],cuts[2]);
        return;
    }
    // Native cell planes reproduce the production evaluator's quadrature. Registry
    // radial, depth and diagonal faces are genuine jumps of the reference metric.
    grid_faces(cuts,d.x,l[0],r[0],0);grid_faces(cuts,d.t,l[1],r[1],1);
    grid_faces(cuts,d.a,l[2],r[2],2);grid_faces(cuts,d.reg_t,l[1],r[1],1);
    if(l[0]!=r[0] && std::max(l[0],r[0])>weld){
        for(int k=0;k<d.reg_e.n;k++){
            const double s=std::exp(-d.pressure*d.reg_e.v[k]);
            if(s>=weld)add_face(cuts,l[0],r[0],s,0);
        }
    }
    sorted(cuts);
    size_t end=cuts.size();const double dr=r[0]-l[0],dt=r[1]-l[1];
    for(size_t k=1;k<end;k++){
        const double lower=cuts[k-1].u,upper=cuts[k].u,mid=.5*(lower+upper);
        const double xm=l[0]+mid*dr,tm=l[1]+mid*dt;
        if(xm<=weld)continue;
        const int i=cell(d.reg_e,-std::log(xm)/d.pressure),j=cell(d.reg_t,tm);
        const double de=d.reg_e.v[i+1]-d.reg_e.v[i],dz=d.reg_t.v[j+1]-d.reg_t.v[j];
        auto face=[&](double u){return (-std::log(l[0]+u*dr)/d.pressure-d.reg_e.v[i])/de+(l[1]+u*dt-d.reg_t.v[j])/dz-1.;};
        double bounds[3]={lower,upper,upper};int pieces=1;
        if(dt!=0&&dr!=0){const double turn=(dr*dz/(d.pressure*de*dt)-l[0])/dr;
            if(turn>lower&&turn<upper){bounds[1]=turn;pieces=2;}}
        for(int p=0;p<pieces;p++){
            double a=bounds[p],b=bounds[p+1];const double fa=face(a),fb=face(b);
            if(fa==0||fb==0||std::signbit(fa)==std::signbit(fb))continue;
            for(int n=0;n<52;n++){double m=.5*(a+b);if(m==a||m==b)break;
                if(std::signbit(face(m))==std::signbit(fa))a=m;else b=m;}
            const double m=.5*(a+b);
            cuts.push_back({m,{-1./(d.pressure*de*(l[0]+m*dr)),1./dz,0.}});
        }
    }
    if(cuts.size()!=end)sorted(cuts);
}

// det(dxi/d physical chart) is degree <=2 in each local coordinate. Use
// tri-quadratic Bernstein coefficients to certify orientation on the WHOLE cell.
// Actual nonpositive samples and inconclusive positive samples have distinct codes.
using Net = std::array<double,27>;
struct Orientation {
    double lower=std::numeric_limits<double>::infinity(), witness=0., margin=0.;
    Vec where{};
    int depth=0;
};

// De Casteljau preserves the convex-hull enclosure on both half intervals.
void split_net(const Net &in, int axis, Net &left, Net &right){
    const int stride=axis==0?9:axis==1?3:1;
    for(int n=0;n<27;n++){
        if((n/stride)%3!=0)continue;
        const double a=in[n], b=in[n+stride], c=in[n+2*stride];
        const double ab=(a+b)*.5, bc=(b+c)*.5, mid=(ab+bc)*.5;
        left[n]=a; left[n+stride]=ab; left[n+2*stride]=mid;
        right[n]=mid; right[n+stride]=bc; right[n+2*stride]=c;
    }
}

int certify_net(const Net &net,const Vec &lo,double span,int depth,Orientation &proof){
    const double lower=*std::min_element(net.begin(),net.end())-proof.margin;
    if(lower>0.){
        proof.lower=std::min(proof.lower,lower);
        return 0;
    }
    int corner=0;
    for(int i=0;i<3;i+=2)for(int j=0;j<3;j+=2)for(int k=0;k<3;k+=2){
        const int n=9*i+3*j+k;
        if(net[n]<net[corner])corner=n;
    }
    if(net[corner]<=proof.margin || depth==12){
        proof.witness=net[corner]; proof.lower=lower; proof.depth=depth;
        proof.where={lo[0]+(corner/9)*span*.5, lo[1]+((corner/3)%3)*span*.5, lo[2]+(corner%3)*span*.5};
        return proof.witness<=0.?2:3;
    }
    const double half=span*.5;
    Net xs[2],ys[2],zs[2];
    split_net(net,0,xs[0],xs[1]);
    for(int i=1;i>=0;i--){
        split_net(xs[i],1,ys[0],ys[1]);
        for(int j=1;j>=0;j--){
            split_net(ys[j],2,zs[0],zs[1]);
            for(int k=1;k>=0;k--){
                int rc=certify_net(zs[k],{lo[0]+i*half,lo[1]+j*half,lo[2]+k*half},half,depth+1,proof);
                if(rc)return rc;
            }
        }
    }
    return 0;
}

int orientation(const Q1 &q,const Vec &width,Orientation &proof){
    // Form columns from edge differences. This removes the constant label
    // before any interpolation and keeps the margin invariant under translation.
    Vec columns[3][3][3];
    double scale=1.;
    for(int axis=0;axis<3;axis++){
        const int d1=(axis==0)?1:0, d2=(axis==2)?1:2;
        Vec edge[2][2];
        double largest=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++){
            int b[3]={0,0,0}; b[d1]=i; b[d2]=j;
            const Vec &a=q.labels[b[0]][b[1]][b[2]];
            b[axis]=1;
            const Vec &z=q.labels[b[0]][b[1]][b[2]];
            for(int c=0;c<3;c++)edge[i][j][c]=(z[c]-a[c])/width[axis];
            largest=std::max(largest,std::sqrt(dot(edge[i][j],edge[i][j])));
        }
        scale*=largest;
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
            const double u=i*.5,v=j*.5;
            for(int c=0;c<3;c++)columns[axis][i][j][c]=
                (1-u)*((1-v)*edge[0][0][c]+v*edge[0][1][c])+
                    u*((1-v)*edge[1][0][c]+v*edge[1][1][c]);
        }
    }
    // Inputs are binary64 spatial labels. Keep the same guard in promoted
    // derivative readers evaluating these saved labels.
    proof.margin=128*0x1p-52*scale;
    Net net;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++){
        const double value=dot(columns[0][j][k],cross(columns[1][i][k],columns[2][i][j]));
        if(!std::isfinite(value)){
            proof.lower=value; proof.witness=value; proof.where={.5*i,.5*j,.5*k};
            return 3;
        }
        net[9*i+3*j+k]=value;
    }
    for(int stride: {9,3,1})for(int n=0;n<27;n++)if((n/stride)%3==0){
        net[n+stride]=2*net[n+stride]-.5*(net[n]+net[n+2*stride]);
    }
    return certify_net(net,{0.,0.,0.},1.,0,proof);
}

struct PointWork {double unit=0.,flux=0.,Jinv=0.,grad_norm=0.,W=0.; Vec value{}; Mat tangent{};};
PointWork material(const Geometry& g,const Vec&local,const Landed&p,bool need_gradient){
    PointWork out{};Vec pos{};for(int k=0;k<3;k++)pos[k]=g.lo[k]+g.width[k]*local[k];
    const auto &d=g.d;const bool core=p.xi[0]<=d.x.v[(d.x.n-1)*2/3];
    const auto prof=profile(d,p.xi[1],true);const Ref ref=reference(d,p.xi[0],p.xi[1],prof,core);
    const double sn=std::sin(p.xi[2]);
    if(!(ref.S>0 && ref.det2>0 && sn>0))throw std::domain_error("reference chart");
    const auto wall=radius(d,pos);const Mat C=physical_basis(wall,d.length,g.width);
    const double detD=determinant(p.D),volume=4*pi*detD*ref.det2*ref.S*ref.S*sn;
    const Mat di=invert(p.D),bi=product(di,ref.inv),F=product(C,bi);
    const double hoop=wall.R*std::sin(pos[2])/(ref.S*sn);
    const double J=determinant(C)*hoop/(detD*ref.S*ref.det2);
    if(!(volume>0 && J>0 && std::isfinite(J)))throw std::domain_error("oriented volume");
    out.Jinv=1./J;out.W=.5*(norm2(F)+hoop*hoop-4.)-std::log(J)+.25*(J-1.)*(J-1.);
    out.unit=volume*out.W;
    // Gradient norm in the orthonormal meridional physical frame of the current state.
    const Mat grad=product(p.D,invert(C));out.grad_norm=std::sqrt(norm2(grad));
    if(!need_gradient)return out;
    const double load=p.chi*volume,vcoef=out.W+1.-.5*J*(J-1.);
    Mat s=product(transpose(F),F);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)s[i][j]=(i==j?vcoef:0.)-s[i][j];
    const Mat gradB=product(s,transpose(bi)),gradH=product(gradB,transpose(p.D));
    const Mat H=invert(ref.inv);
    out.tangent=product(transpose(H),gradB);
    for(auto &r:out.tangent)for(double &z:r)z*=load;
    const double ss=ref.S,sx=H[0][0],st=H[0][1],polar=gradH[2][2]+(vcoef-hoop*hoop)/ss;
    out.value={polar*sx,polar*st,(vcoef-hoop*hoop)*std::cos(p.xi[2])/sn};
    if(core){
        const double slope=d.core_ratio*prof.bt;
        out.value[0]+=gradH[0][1]*slope;
        out.value[1]+=gradH[0][0]*slope+gradH[0][1]*p.xi[0]*d.core_ratio*prof.btt;
    }else if(d.collar_coeff){
        for(int i=0;i<2;i++){
            out.value[i]+=gradH[0][0]*ref.s_hess[i]+gradH[0][1]*ref.s_hess[i+1]
                         +gradH[1][0]*ref.z_hess[i]+gradH[1][1]*ref.z_hess[i+1];
        }
    }else out.value[0]-=(gradH[0][0]*sx+gradH[1][0]*H[1][0])/p.xi[0];
    for(double &v:out.value)v*=load;
    const double rateR=wall.grad[0]*g.width[0]/wall.R;
    const Mat gradC=product(F,transpose(bi));
    double moving=hoop*hoop*rateR+(-1.+.5*J*(J-1.))*(wall.hess[0][0]/wall.grad[0]*g.width[0]+2.*rateR);
    for(int k=0;k<3;k++)moving+=gradC[0][k]*wall.hess[k][0]*g.width[0]*g.width[k];
    moving+=gradC[2][2]*wall.grad[0]*g.width[0]*g.width[2];
    out.flux=load*moving;
    return out;
}
}

// Standalone whole-chart entry. Physics and guard call the same cell predicate.
extern "C" int r4_chart_orientation(Grid x,Grid t,Grid a,const double *labels,
                                      const unsigned char *active,double *out){
    const Grid axes[]={x,t,a};
    out[0]=out[1]=std::numeric_limits<double>::infinity();
    for(int i=0;i<x.n-1;i++)for(int j=0;j<t.n-1;j++)for(int k=0;k<a.n-1;k++){
        const int cell=(i*(t.n-1)+j)*(a.n-1)+k;
        Q1 q{};
        for(int u=0;u<2;u++)for(int v=0;v<2;v++)for(int w=0;w<2;w++){
            const int n=((i+u)*t.n+j+v)*a.n+k+w;
            for(int c=0;c<3;c++)q.labels[u][v][w][c]=labels[3*n+c];
        }
        const int index[]={i,j,k};
        Vec width;
        for(int c=0;c<3;c++)width[c]=axes[c].v[index[c]+1]-axes[c].v[index[c]];
        euler::Orientation proof;
        const int rc=euler::orientation(q,width,proof);
        if(rc || proof.lower<out[0]){
            out[0]=proof.lower; out[2]=i; out[3]=j; out[4]=k;
            for(int c=0;c<3;c++)out[5+c]=proof.where[c];
            out[8]=proof.witness; out[9]=proof.margin; out[10]=proof.depth;
        }
        if(rc)return rc;
        if(active[cell])out[1]=std::min(out[1],proof.lower);
    }
    return 0;
}

extern "C" int r4_eulerian_material(const Data*dp,const double *labels,const double *order,
           Grid rule,const double *weights,double *output,double *grad,double *loads,double *diagnostics,double *intact){
    using namespace euler;const Data &d=*dp;const int nx=d.field_x.n,nt=d.field_t.n,na=d.field_a.n;
    auto id=[&](int x,int z,int a){return (x*nt+z)*na+a;};
    long double energy=0.;std::fill(diagnostics,diagnostics+5,std::numeric_limits<double>::infinity());
    Cuts cuts;cuts.reserve(60);std::vector<double> cut_work;
    const double weld=d.x.v[(d.x.n-1)*2/3];
    for(int i=0;i<nx-1;i++)for(int j=0;j<nt-1;j++)for(int k=0;k<na-1;k++){
        Q1 q;bool has_order=false;int ids[2][2][2];
        for(int X=0;X<2;X++)for(int Z=0;Z<2;Z++)for(int A=0;A<2;A++){
            int n=ids[X][Z][A]=id(i+X,j+Z,k+A);q.order[X][Z][A]=order[n];has_order|=order[n]>0;
            for(int c=0;c<3;c++)q.labels[X][Z][A][c]=labels[3*n+c];
        }
        Geometry g{d,{d.field_x.v[i],d.field_t.v[j],d.field_a.v[k]},
                   {d.field_x.v[i+1]-d.field_x.v[i],d.field_t.v[j+1]-d.field_t.v[j],d.field_a.v[k+1]-d.field_a.v[k]}};
        Orientation proof;
        int rc=orientation(q,g.width,proof);
        if(rc){diagnostics[2]=proof.lower;diagnostics[5]=i;diagnostics[6]=j;diagnostics[7]=k;
            for(int a=0;a<3;a++)diagnostics[8+a]=proof.where[a];
            diagnostics[11]=q.at(proof.where[0],proof.where[1],proof.where[2]).chi;
            return rc;}
        diagnostics[2]=std::min(diagnostics[2],proof.lower);
        if(has_order)diagnostics[3]=std::min(diagnostics[3],proof.lower);
        const bool need=has_order&&grad;
        for(int T=0;T<rule.n;T++)for(int A=0;A<rule.n;A++){
            double v=rule.v[T],a=rule.v[A],wa=weights[T]*weights[A];
            double bta[2][2],bt[2][2],ba[2][2];
            for(int Z=0;Z<2;Z++)for(int K=0;K<2;K++){
                bta[Z][K]=(Z?v:1-v)*(K?a:1-a);bt[Z][K]=(2*Z-1)*(K?a:1-a);ba[Z][K]=(Z?v:1-v)*(2*K-1);
            }
            const Landed l=q.at(0,v,a),r=q.at(1,v,a);Vec delta{};Mat dD{};
            for(int c=0;c<3;c++){delta[c]=r.xi[c]-l.xi[c];for(int h=0;h<3;h++)dD[c][h]=r.D[c][h]-l.D[c][h];}
            const double dc=r.chi-l.chi;
            radial_cuts(cuts,d,l.xi,r.xi,weld);cut_work.assign(cuts.size(),0.);
            for(size_t z=1;z<cuts.size();z++){
                const double du=cuts[z].u-cuts[z-1].u;if(!(du>0))continue;
                for(int XG=0;XG<rule.n;XG++){
                    const double f=rule.v[XG],u=cuts[z-1].u+du*f,w=wa*weights[XG],wt=w*du;
                    auto p=l;p.chi+=u*dc;
                    for(int c=0;c<3;c++){p.xi[c]+=u*delta[c];for(int h=0;h<3;h++)p.D[c][h]+=u*dD[c][h];}
                    PointWork state;
                    try{state=material(g,{u,v,a},p,need);}catch(std::domain_error&){return 4;}
                    if(!std::isfinite(state.unit+state.grad_norm))return 5;
                    if(state.Jinv<diagnostics[0]){diagnostics[0]=state.Jinv;diagnostics[5]=i;diagnostics[6]=j;diagnostics[7]=k;
                        diagnostics[8]=u;diagnostics[9]=v;diagnostics[10]=a;diagnostics[11]=p.chi;diagnostics[12]=state.grad_norm;diagnostics[13]=state.W;}
                    if(p.chi>0&&state.Jinv<diagnostics[1]){diagnostics[1]=state.Jinv;diagnostics[14]=i;diagnostics[15]=j;diagnostics[16]=k;
                        diagnostics[17]=u;diagnostics[18]=v;diagnostics[19]=a;diagnostics[20]=p.chi;diagnostics[21]=state.grad_norm;diagnostics[22]=state.W;}
                    if(p.chi>=.5 && state.W>intact[4*j]){
                        intact[4*j]=state.W;intact[4*j+1]=d.length*(g.lo[1]+v*g.width[1]);
                        intact[4*j+2]=p.chi;intact[4*j+3]=state.Jinv;
                    }
                    diagnostics[4]=std::min(diagnostics[4],1./state.Jinv);
                    energy+=wt*p.chi*state.unit;
                    double tangent_u=state.flux+dot(state.value,delta)+state.unit*dc;
                    if(need){
                        for(int c=0;c<3;c++)tangent_u+=dot(state.tangent[c],dD[c]);
                        cut_work[z-1]+=w*(-state.unit*p.chi+du*(1-f)*tangent_u);
                        cut_work[z]+=w*(state.unit*p.chi+du*f*tangent_u);
                    }
                    for(int XX=0;XX<2;XX++)for(int Z=0;Z<2;Z++)for(int K=0;K<2;K++){
                        int idx=ids[XX][Z][K];const double bx=XX?u:1-u,sh=bx*bta[Z][K];
                        loads[idx]+=wt*sh*state.unit;
                        if(need){Vec ds{(2*XX-1)*bta[Z][K],bx*bt[Z][K],bx*ba[Z][K]};
                            for(int c=0;c<3;c++)grad[3*idx+c]+=wt*(dot(state.tangent[c],ds)+state.value[c]*sh);}
                    }
                }
            }
            if(need)for(size_t n=1;n+1<cuts.size();n++){
                const Cut &cut=cuts[n];const double denom=dot(cut.normal,delta);
                if(denom==0){if(cut_work[n]!=0)return 6;continue;}
                const double adj=-cut_work[n]/denom;
                for(int X=0;X<2;X++)for(int Z=0;Z<2;Z++)for(int K=0;K<2;K++){
                    double factor=adj*bta[Z][K]*(X?cut.u:1-cut.u);int idx=ids[X][Z][K];
                    for(int c=0;c<3;c++)grad[3*idx+c]+=factor*cut.normal[c];}
            }
        }
    }
    *output=double(energy);return 0;
}
