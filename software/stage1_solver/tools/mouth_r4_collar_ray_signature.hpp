// Enumerate exact reference triangle sequences with the promoted production
// ray cuts and production Gauss sites before and after a directional probe.
#include "mouth_r4_segmented_endpoint_forces.hpp"
extern "C" int direction_changed_rays(const Data *pd, const Real *labels, const Real *delta,
        Grid rule, const int *active, int *flags) {
    const Data &d=*pd;
    const int nx=d.field_x.n, nt=d.field_t.n, na=d.field_a.n;
    int count=0;
    std::vector<endpoint_depth::Sample> base_seq, moved_seq;
    for (int i=0;i+1<nx;++i) for (int j=0;j+1<nt;++j) for(int k=0;k+1<na;++k) {
        const int cell=(i*(nt-1)+j)*(na-1)+k;
        if (!active[cell]) continue;
        Q1 q{}, moved{};
        for (int x=0;x<2;++x) for(int z=0;z<2;++z) for(int a=0;a<2;++a) {
            const int n=((i+x)*nt+j+z)*na+k+a;
            for(int c=0;c<3;++c){
                q.labels[x][z][a][c]=labels[3*n+c];
                moved.labels[x][z][a][c]=labels[3*n+c]+delta[3*n+c];
            }
        }
        for(int t=0;t<rule.n;++t) for(int a=0;a<rule.n;++a) {
            endpoint_depth::branches(d,q,rule,t,a,base_seq);
            endpoint_depth::branches(d,moved,rule,t,a,moved_seq);
            if (base_seq!=moved_seq) {flags[(cell*rule.n+t)*rule.n+a]=1;++count;}
        }
    }
    return count;
}
