// Physical-volume change-of-variables control on EACH registered reference patch.
// This evaluates the inverse tensor on the imaged material quadrature. It does
// not assert that the two patches define one continuous Euclidean reference map.
#include "mouth_r4_refined_solid.cpp"

struct TensorEnergy {
    double inverse_volume;
    double energy;
    Mat deformation;
};

Mat inverse_matrix(const Mat &m) {
    Mat inverse = transpose(cofactor(m));
    const double det = determinant(m);
    for (auto &row : inverse) for (double &v : row) v /= det;
    return inverse;
}

bool inverse_energy(const Mat &G, double hoop_inverse, double chi, TensorEnergy &out) {
    const double d = determinant(G) * hoop_inverse;
    if (!(d > 0 && std::isfinite(d) && chi > 0)) return false;
    const Mat F = inverse_matrix(G);
    const double h = 1. / hoop_inverse;
    const double J = 1. / d;
    out = {d, chi*d*(.5*(norm2(F)+h*h-4) + std::log(d) + .25*(J-1)*(J-1)), F};
    return std::isfinite(out.energy);
}

struct ImageCells : AngularEdgeStar {
    int radial_count, polar_count;
    mutable int parent = -1, child = -1;
    bool include(int j,int i,int k,int tri) const {
        parent = ((j*radial_count+i)*polar_count+k)*2+tri;
        return true;
    }
    void enrich(TriSlice &s,const Data &d,int j,int i,int k,int tri,double t,double a) const {
        AngularEdgeStar::enrich(s,d,j,i,k,tri,t,a);
        child = refined && incident(j,i,k,tri) && a>=seam ? 1:0;
    }
};

struct ImageEnergy {
    const ImageCells &cells;
    double *elements, *diagnostics;
    int operator()(const Data &d,Results &out,const Profile &initial,const TriSlice &s,
                   double t,double a,double x,double weight,bool core,double &den,double &denx) const {
        const int code=sample(d,out,initial,s,t,a,x,weight,core,den,denx);
        if (code) return code;
        Mat grad=s.M;
        Vec landed=s.X;
        for (int i=0;i<3;i++) {
            landed[i]+=grad[i][0]*(x-s.x0);
            grad[i][2]+=s.Max[i]*(x-s.x0);
        }
        const auto bx=basis(d.field_x,landed[0]),bt=basis(d.field_t,landed[1]),ba=basis(d.field_a,landed[2]);
        double chi=0;
        for(int i=0;i<2;i++)for(int k=0;k<2;k++)for(int j=0;j<2;j++)
            chi+=d.correction[field_id(d,bx.i+i,ba.i+k,bt.i+j)]*bx.w[i]*ba.w[k]*bt.w[j];
        if (chi==0) return 0;
        const auto ref=reference(d,x,t,initial,core);
        const auto wall=radius(d,landed);
        Mat forward=product(Mat{wall.grad,Vec{0,d.length,0},Vec{0,0,wall.R}},product(grad,ref.inv));
        // Select a right-handed local reference frame. A constant reflection
        // preserves all singular values and energy on that patch.
        if(determinant(forward)<0) for(auto &row:forward) row[1]=-row[1];
        const double hoop=wall.R*std::sin(landed[2])/(ref.S*std::sin(a));
        const Mat inverse=inverse_matrix(forward);
        TensorEnergy e;
        if(!inverse_energy(inverse,1/hoop,chi,e))return 4;
        // Physical volume directly from the deformed chart; independent of
        // reference-volume weight and of the computed inverse determinant.
        const double physical=4*pi*d.length*wall.grad[0]*wall.R*wall.R*std::sin(landed[2])*determinant(grad)*weight;
        if (!(physical>0))return 5;
        const int idx=2*(2*cells.parent+cells.child);
        elements[idx]+=weight*den;
        elements[idx+1]+=physical*e.energy;
        const double nativeJ=wall.grad[0]*d.length*wall.R*determinant(grad)*hoop/(ref.det2*ref.S);
        diagnostics[0]+=1;
        diagnostics[1]=std::max(diagnostics[1],std::abs(1/e.inverse_volume-nativeJ));
        Mat delta=e.deformation;
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)delta[i][j]-=forward[i][j];
        diagnostics[2]=std::max(diagnostics[2],std::sqrt(norm2(delta)+std::pow(1/(1/hoop)-hoop,2)));
        if(diagnostics[3]==0){
            Mat negative=inverse;
            for(double &v:negative[0])v=-v;
            if(inverse_energy(negative,1/hoop,chi,e))return 6;
            diagnostics[3]=chi;diagnostics[4]=determinant(negative)/hoop;
            for(int i=0;i<3;i++)diagnostics[5+i]=landed[i];
        }
        return 0;
    }
};

extern "C" int image_energy(const Data *d,Results *out,const double *surplus,double seam,
                            int refined,double *elements,double *diagnostics) {
    ImageCells cells{{surplus,seam,refined!=0,0,0,0,0,0},d->x.n-1,d->a.n-1};
    return integrate_solid(d,out,cells,ImageEnergy{cells,elements,diagnostics});
}
