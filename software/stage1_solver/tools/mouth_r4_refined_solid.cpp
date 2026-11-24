// Conforming refinement of a single angular edge and its triangle star.
// The dependent nodes on the patch boundary retain the exact coarse trace.
#include "mouth_free_shape_solid.cpp"

struct AngularEdgeStar {
    static constexpr int depth=5, radial=18, polar=0;
    const double *surplus;
    double seam;
    bool refined;
    int region;  // 0 all, 1 six-prism star, 2 one parent prism
    int target_j,target_i,target_k,target_half;

    static bool incident(int j,int i,int k,int half) {
        if(k!=polar)return false;
        int dz=depth-j,dx=radial-i;
        if(dz<0 || dx<0 || dz>1 || dx>1)return false;
        return half ? dx+dz>=1 : dx+dz<=1;
    }
    bool include(int j,int i,int k,int half) const {
        if(region==1)return incident(j,i,k,half);
        if(region==2)return j==target_j && i==target_i && k==target_k && half==target_half;
        return true;
    }
    void enrich(TriSlice &s,const Data &d,int j,int i,int k,int half,double t,double a) const {
        if(!refined || !incident(j,i,k,half))return;
        // Read meridional basis functions from their known triangle rather
        // than dividing a tensor-product basis by a vanishing polar weight.
        const double dx=d.x.v[i+1]-d.x.v[i],dt=d.t.v[j+1]-d.t.v[j];
        const double u=(s.x0-d.x.v[i])/dx,v=(t-d.t.v[j])/dt;
        double lambda,lambda_x,lambda_t;
        if((i==radial && j==depth && half==0) ||
           (i+1==radial && j+1==depth && half==1)) {
            const double sign=half ? 1. : -1.;
            lambda=half ? u+v-1. : 1.-u-v;
            lambda_x=sign/dx;
            lambda_t=sign/dt;
        } else if((half==0 && i+1==radial) || (half==1 && i==radial)) {
            lambda=half ? 1.-u : u;
            lambda_x=(half ? -1. : 1.)/dx;
            lambda_t=0.;
        } else {
            lambda=half ? 1.-v : v;
            lambda_x=0.;
            lambda_t=(half ? -1. : 1.)/dt;
        }
        const double a0=d.a.v[k],a1=d.a.v[k+1];
        const double h=a<seam ? (a-a0)/(seam-a0) : (a1-a)/(a1-seam);
        const double hp=a<seam ? 1./(seam-a0) : -1./(a1-seam);
        const int n=s.count++;
        s.id[n]=d.x.n*d.t.n*d.a.n;
        s.n[n]=lambda*h; s.nx[n]=lambda_x*h; s.nt[n]=lambda_t*h;
        s.na[n]=lambda*hp; s.nax[n]=lambda_x*hp;
        for(int c=0;c<3;c++) {
            s.X[c]+=s.n[n]*surplus[c];
            s.M[c][0]+=s.nx[n]*surplus[c];
            s.M[c][1]+=s.nt[n]*surplus[c];
            s.M[c][2]+=s.na[n]*surplus[c];
            s.Max[c]+=s.nax[n]*surplus[c];
        }
    }
};

extern "C" int refined_solid_eval(const Data *data,Results *res,const double *surplus,
                                  double seam,int refined,int region,int j,int i,int k,int half,
                                  double *surplus_gradient) {
    AngularEdgeStar star{surplus,seam,refined!=0,region,j,i,k,half};
    if(!res->gradient)return integrate_solid(data,res,star);
    // The existing binding holds a six-node coarse chart. A seventh local
    // basis at the shared new node contributes one global 3-vector.
    const int size=3*data->x.n*data->t.n*data->a.n;
    std::vector<double> gradient(size+3,0.);
    Results result=*res;
    result.gradient=gradient.data();
    int err=integrate_solid(data,&result,star);
    std::copy_n(gradient.data(),size,res->gradient);
    std::copy_n(gradient.data()+size,3,surplus_gradient);
    return err;
}
