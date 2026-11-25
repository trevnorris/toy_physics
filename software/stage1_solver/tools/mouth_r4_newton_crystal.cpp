// Exact fixed-order directional derivative of the verified smooth corner gauge.
// Reuse the pointwise barrier's extended-precision determinant and signed work.
#include "mouth_r4_melt_barrier.cpp"

extern "C" int corner_tangent(int count,const int* ids,const double* width,
        const double* weight,const double* labels,const double* step,
        double* action,double* diagonal) {
    for (int c=0;c<count;++c) for (int v=0;v<8;++v) {
        const int node=ids[8*c+v];
        std::array<Edge,3> a,da,cof;
        std::array<int,3> next;
        Edge dx,center{};
        for (int d=0;d<3;++d) {
            const int bit=4>>d;
            next[d]=ids[8*c+(v^bit)];
            dx[d]=((v&bit)?-1.L:1.L)/width[3*c+d];
            for (int k=0;k<3;++k) {
                a[d][k]=(static_cast<long double>(labels[3*next[d]+k])-labels[3*node+k])*dx[d];
                if(step) da[d][k]=(static_cast<long double>(step[3*next[d]+k])-step[3*node+k])*dx[d];
            }
        }
        const auto det=triple(a[0],a[1],a[2]);
        if (!(det>0) || !std::isfinite(det)) return 1;
        long double ddet=0;
        for (int d=0;d<3;++d) {
            cof[d]=cross3(a[(d+1)%3],a[(d+2)%3]);
            if(step) for (int k=0;k<3;++k) ddet+=cof[d][k]*da[d][k];
        }
        for (int d=0;d<3;++d) {
            const int i=(d+1)%3,j=(d+2)%3;
            Edge dc1{},dc2{};
            if(step) { dc1=cross3(da[i],a[j]); dc2=cross3(a[i],da[j]); }
            for (int k=0;k<3;++k) {
                if(step) {
                    const auto h=weight[c]*dx[d]*((1-1/det)*(dc1[k]+dc2[k])+cof[d][k]*ddet/(det*det));
                    action[3*next[d]+k]+=h;
                    action[3*node+k]-=h;
                }
                if(diagonal) {
                    const auto z=dx[d]*cof[d][k]/det;
                    diagonal[3*next[d]+k]+=weight[c]*z*z;
                    center[k]-=z;
                }
            }
        }
        if(diagonal) for(int k=0;k<3;++k) diagonal[3*node+k]+=weight[c]*center[k]*center[k];
    }
    return 0;
}
