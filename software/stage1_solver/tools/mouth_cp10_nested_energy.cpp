// Nested spatial label interpolation evaluated on unchanged checkpoint-10 material
// quadrature. Only the sampling lattice of the inverse labels changes.
#include "mouth_cp10_imaged_labels.cpp"

namespace {
enum Column {point_count,native_energy,projected_energy,gradient_only_energy,absolute_error,columns};
enum Correspondence {same_prism,same_side_different_prism,different_side,groups};

struct SpatialLabels {
    Grid x,t,a;
    const double *nodes;
    void interpolate(const Vec &p,Vec &xi,Mat &derivative) const {
        const auto r=basis(x,p[0]),z=basis(t,p[1]),b=basis(a,p[2]);
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            const double *node=nodes+3*(((r.i+i)*t.n+z.i+j)*a.n+b.i+k);
            const double w=r.w[i]*z.w[j]*b.w[k];
            const Vec g{r.v[i]*z.w[j]*b.w[k],r.w[i]*z.v[j]*b.w[k],r.w[i]*z.w[j]*b.v[k]};
            for(int c=0;c<3;c++){
                xi[c]+=node[c]*w;
                for(int l=0;l<3;l++)derivative[c][l]+=node[c]*g[l];
            }
        }
    }
};
struct Energy { double J, weighted; };
Energy constitutive(const Data &d,const Vec &xi,const Mat &forward_chart,const Mat &frame,
                    double chart_J,double hoop_radius,double dv,double chi,bool core) {
    const Ref ref=reference(d,xi[0],xi[1],profile(d,xi[1],true),core);
    const Mat F=product(frame,product(forward_chart,ref.inv));
    const double h=hoop_radius/(ref.S*std::sin(xi[2]));
    const double J=determinant(frame)*chart_J*h/(ref.det2*ref.S);
    if(!(J>0 && std::isfinite(J)))return {J,NAN};
    const double W=.5*(norm2(F)+h*h-4)-std::log(J)+.25*(J-1)*(J-1);
    return {J,chi*dv*W/J};
}
int material_element(const Data &d,const Vec &q) {
    const int i=cell(d.x,q[0]),j=cell(d.t,q[1]),k=cell(d.a,q[2]);
    const bool tri=(q[0]-d.x.v[i])/(d.x.v[i+1]-d.x.v[i])+
                   (q[1]-d.t.v[j])/(d.t.v[j+1]-d.t.v[j])>1.;
    return ((j*(d.x.n-1)+i)*(d.a.n-1)+k)*2+tri;
}
struct Readout {
    const MaterialBin &source;
    const SpatialLabels &labels;
    long double *by_element;
    double *stat,*witness;
    int operator()(const Data &d, Results &out, const Profile &initial, const TriSlice &s,
                   double t,double a,double x,double w,bool core,double &den,double &denx) const {
        int rc=sample(d,out,initial,s,t,a,x,w,core,den,denx);
        if(rc)return rc;
        const Vec original{x,t,a};
        Vec p=s.X;
        Mat m=s.M;
        for(int c=0;c<3;c++){
            p[c]+=s.M[c][0]*(x-s.x0);
            m[c][2]+=s.Max[c]*(x-s.x0);
        }
        const auto rr=basis(d.field_x,p[0]),tt=basis(d.field_t,p[1]),aa=basis(d.field_a,p[2]);
        double chi=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++)
            chi+=rr.w[i]*tt.w[j]*aa.w[k]*d.correction[field_id(d,rr.i+i,aa.i+k,tt.i+j)];
        if(chi==0.)return 0;
        Vec xi{};
        Mat grad{};
        labels.interpolate(p,xi,grad);
        const double dg=determinant(grad),mj=determinant(m);
        for(int c=0;c<3;c++){witness[c]=p[c];witness[c+3]=original[c];witness[c+6]=xi[c];}
        witness[9]=chi;witness[10]=dg;witness[11]=mj;
        if(!(dg>0 && mj>0 && std::isfinite(dg)))return 4;
        Mat inv=transpose(cofactor(grad));
        for(auto &row:inv)for(double &val:row)val/=dg;
        const Map wall=radius(d,p);
        const Mat frame{wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
        const double hs=wall.R*std::sin(p[2]),dv=w*4*pi*determinant(frame)*hs*mj;
        const bool projected_core=xi[0]<=d.x.v[(d.x.n-1)*2/3];
        const Energy mapped=constitutive(d,xi,inv,frame,1./dg,hs,dv,chi,projected_core);
        const Energy hybrid=constitutive(d,original,inv,frame,1./dg,hs,dv,chi,core);
        const Energy exact=constitutive(d,original,m,frame,mj,hs,dv,chi,core);
        const double old=den*w;
        if(!std::isfinite(mapped.weighted+hybrid.weighted+exact.weighted))return 5;
        const int g=(projected_core!=core)?different_side:(material_element(d,xi)==source.index?same_prism:same_side_different_prism);
        long double *bin=by_element+(source.index*groups+g)*columns;
        bin[point_count]++;
        bin[native_energy]+=old;
        bin[projected_energy]+=mapped.weighted;
        bin[gradient_only_energy]+=hybrid.weighted;
        bin[absolute_error]+=std::abs(mapped.weighted-old);
        stat[0]=std::min(stat[0],1./mapped.J);
        stat[1]=std::max(stat[1],1./mapped.J);
        stat[2]=std::max(stat[2],std::abs(exact.weighted-old));
        stat[3]+=std::abs(exact.weighted-old);
        return 0;
    }
};
}

extern "C" int nested_image_energy(const Data *d,Results *out,const double *x,int nx,const double *t,int nt,
                                  const double *a,int na,const double *xi,double *table,double *stat,double *witness){
    MaterialBin bin{{},-1,d->x.n-1,d->a.n-1};
    const int count=(d->t.n-1)*bin.nr*bin.na*2*groups*columns;
    SpatialLabels grid{{nx,x},{nt,t},{na,a},xi};
    std::vector<long double> values(count,0.);
    stat[0]=std::numeric_limits<double>::infinity();
    const int rc=integrate_solid(d,out,bin,Readout{bin,grid,values.data(),stat,witness});
    for(int k=0;k<count;k++)table[k]=double(values[k]);
    return rc;
}
