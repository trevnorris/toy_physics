// Compare the two strains at the very same physical points, weighted by the
// image of the original material quadrature. No inverse cells are re-fitted.
#include "mouth_free_shape_solid.cpp"

struct MaterialBin : OriginalPrisms {
    mutable int index=-1;
    int nr,na;
    bool include(int j,int i,int k,int triangle) const {
        index=((j*nr+i)*na+k)*2+triangle;
        return true;
    }
};

struct CompareImages {
    const MaterialBin &bin;
    const double *labels;
    long double *energy, *groups;
    double *order_load;
    int operator()(const Data &d,Results &out,const Profile &initial,const TriSlice &s,
                   double t,double a,double x,double w,bool source_core,double &den,double &denx) const {
        const int rc=sample(d,out,initial,s,t,a,x,w,source_core,den,denx);
        if(rc)return rc;
        Vec p=s.X;
        Mat m=s.M;
        for(int c=0;c<3;c++){
            p[c]+=s.M[c][0]*(x-s.x0);
            m[c][2]+=s.Max[c]*(x-s.x0);
        }
        const auto r=basis(d.field_x,p[0]),z=basis(d.field_t,p[1]),b=basis(d.field_a,p[2]);
        Vec y{};
        Mat grad{};
        double chi=0.;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
            const double h=r.w[i]*z.w[j]*b.w[k];
            chi+=h*d.correction[field_id(d,r.i+i,b.i+k,z.i+j)];
            const int n=((r.i+i)*d.field_t.n+z.i+j)*d.field_a.n+b.i+k;
            for(int c=0;c<3;c++){
                const double label=labels[3*n+c];
                y[c]+=h*label;
                grad[c][0]+=label*r.v[i]*z.w[j]*b.w[k];
                grad[c][1]+=label*r.w[i]*z.v[j]*b.w[k];
                grad[c][2]+=label*r.w[i]*z.w[j]*b.v[k];
            }
        }
        if(chi==0.)return 0;
        const double dg=determinant(grad);
        if(!(dg>0 && determinant(m)>0))return 4;
        const bool label_core=y[0]<=d.x.v[(d.x.n-1)*2/3];
        auto ref=reference(d,y[0],y[1],profile(d,y[1],true),label_core);
        auto wall=radius(d,p);
        Mat inv=transpose(cofactor(grad));
        for(auto &row:inv) for(double &v:row)v/=dg;
        Mat frame={wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
        Mat F=product(frame,product(inv,ref.inv));
        const double hoop=wall.R*std::sin(p[2])/(ref.S*std::sin(y[2]));
        const double J=determinant(frame)*hoop/(ref.det2*ref.S*dg);
        if(!(J>0 && std::isfinite(J)))return 5;
        const double W=.5*(norm2(F)+hoop*hoop-4)-std::log(J)+.25*(J-1)*(J-1);
        const double dv=w*4*pi*determinant(frame)*wall.R*std::sin(p[2])*determinant(m);
        const double projected=chi*dv/J*W, original=w*den;
        energy[2*bin.index]+=original;
        energy[2*bin.index+1]+=projected;
        int idx=(label_core!=source_core)?3:0;
        groups[idx]+=1;
        groups[idx+1]+=original;
        groups[idx+2]+=projected;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++)
            order_load[field_id(d,r.i+i,b.i+k,z.i+j)]+=dv/J*W*r.w[i]*z.w[j]*b.w[k];
        return 0;
    }
};

extern "C" int compare_images(const Data *data,Results *out,const double *labels,double *elements,double *report,double *order_load){
    MaterialBin bin{{},-1,data->x.n-1,data->a.n-1};
    const int count=2*(data->t.n-1)*bin.nr*bin.na*2;
    std::vector<long double> energy(count,0.);
    long double groups[6]{};
    int rc=integrate_solid(data,out,bin,CompareImages{bin,labels,energy.data(),groups,order_load});
    if(rc)return rc;
    for(int i=0;i<count;i++)elements[i]=double(energy[i]);
    for(int i=0;i<6;i++)report[i]=double(groups[i]);
    return 0;
}
