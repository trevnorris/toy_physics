// Reuse the exact archived release and finite-contact kernel, including its original
// translation unit, so the preceding experiment remains reproducible by its hash.
#define main archived_ripple_experiment
#include "mouth_middepth_rearrangement.cpp"
#undef main

struct WallReading {
    double upstream=0, downstream=0, fresh=0, total_weight=0;
    double energy_force=0, direct_force=0;
};

struct BendHold {
    const Patch& patch;
    Vec shape;
    bool follow_feed;
    double end;
    const Vec amplitudes;

    BendHold(const Patch& p, const char* table, const std::string& feed, double edge):
        patch(p), follow_feed(feed=="following"), end(edge), amplitudes{0,p.h,-p.h,p.h/2,-p.h/2} {
        std::ifstream input(table);
        double v;
        while(input>>v) shape.push_back(v);
        if(shape.size()!=static_cast<size_t>((p.length+2)*2*p.sub+1) || !(-6<end && end<0))
            throw std::runtime_error("Incomplete bend shape or edge outside central gauge");
    }
    double wall_shape(int i, int phase) const {
        return shape[(i/patch.lanes)*2*patch.sub+((i%patch.lanes)%2)*patch.sub+phase];
    }
    Vec velocity(const Vec& sites, double amplitude, int phase) const {
        auto f=patch.forces(sites);
        for(int i=0;i<patch.count();++i) {
            if(amplitude!=0 && patch.anchor(i)!=0)
                f.wall[i]=-patch.anchor(i)*patch.wall_law(sites[i]-amplitude*wall_shape(i,phase))[1];
            f.internal[i]+=f.wall[i];
        }
        return f.internal;
    }
    void advance(Vec& sites, double amplitude, int j) const {
        double dt=patch.period/patch.sub;
        Vec mid=sites;
        const auto first=velocity(sites,amplitude,2*j);
        for(int i=0;i<patch.count();++i) mid[i]+=.5*dt*first[i];
        const auto slope=velocity(mid,amplitude,2*j+1);
        for(int i=0;i<patch.count();++i) sites[i]+=dt*slope[i];
    }
    void inject(Vec& sites, double amp) const {
        std::move_backward(sites.begin(),sites.end()-patch.lanes,sites.end());
        for(int i=0;i<patch.lanes;++i) sites[i]=follow_feed ? amp*wall_shape(i,0) : 0.;
    }
    WallReading read(const Vec& sites, const Vec& zero, double amp, int phase, int arrivals,
                     bool verify_energy=false) const {
        WallReading r;
        const double trial=.00001;
        for(int i=0;i<patch.count();++i) {
            const double pos=patch.reference_x(i)+phase/(2.*patch.sub)-patch.length*.5;
            if(pos< -6 || pos>=6 || !patch.anchor(i)) continue;
            const double factor=patch.weight()*patch.anchor(i);
            const double slip=sites[i]-amp*wall_shape(i,phase);
            const double inc=factor*(patch.wall_law(slip)[1]-patch.wall_law(zero[i])[1]);
            if(pos<end) r.upstream+=inc;
            else r.downstream+=inc;
            r.total_weight+=1;
            if(i/patch.lanes<arrivals) r.fresh+=1;
            if(verify_energy) {
                // Force on the wall is minus the derivative with respect to a uniform
                // virtual displacement of these wall ends, with sites frozen.
                const double ep=patch.wall_law(slip-trial)[0]-patch.wall_law(zero[i]-trial)[0];
                const double em=patch.wall_law(slip+trial)[0]-patch.wall_law(zero[i]+trial)[0];
                r.energy_force-=factor*(ep-em)/(2*trial);
                r.direct_force+=inc;
            }
        }
        if(r.total_weight==0) throw std::runtime_error("Empty wall gauge");
        return r;
    }
    std::array<double,4> slopes(const Copies& states, int phase, int n,
                                double& half_error, double& symmetry, double& fresh) const {
        std::array<WallReading,5> r;
        for(int k=0;k<5;++k) r[k]=read(states[k],states[0],amplitudes[k],phase,n);
        const auto net=[](const WallReading& a){return a.upstream+a.downstream;};
        const double scale=std::max(1e-16,std::abs(net(r[1]))+std::abs(net(r[2])));
        if(net(r[0])!=0.) throw std::runtime_error("Nonzero subtracted baseline");
        half_error=std::abs((net(r[1])-net(r[2]))*.5-(net(r[3])-net(r[4])));
        symmetry=std::abs(net(r[1])+net(r[2]))/scale;
        fresh=r[0].fresh/r[0].total_weight;
        return {net(r[1]),net(r[2]),(r[1].upstream-r[2].upstream)*.5,
                (r[1].downstream-r[2].downstream)*.5};
    }
    void run() const {
        std::mt19937 gen(std::abs(patch.seed)+1234);
        std::uniform_real_distribution<double> rand(-1.,1.);
        Vec s(patch.count());
        for(double& v:s) v=1e-5*rand(gen);
        Vec prior, two_back;
        for(int n=0;n<6*patch.length;++n) {
            two_back=prior;prior=patch.center(s);
            for(int j=0;j<patch.sub;++j) advance(s,0.,j);
            inject(s,0.);
        }
        Vec center=patch.center(s), diff=center;
        for(size_t i=0;i<center.size();++i) diff[i]-=two_back[i];
        const double lag_two=rms(diff);
        for(size_t i=0;i<center.size();++i) diff[i]=center[i]-prior[i];
        std::cout<<"{\"test\":\"release\",\"replacements\":6,\"center\":";array(center);
        std::cout<<",\"lag_one\":"<<rms(diff)<<",\"lag_two\":"<<lag_two<<"}\n";
        if(!std::isfinite(lag_two) || lag_two>1e-5) throw std::runtime_error("No recurrent release");
        Copies states(amplitudes.size(),s);
        double energy_error=0,half_error=0,symmetry=0;
        auto sample=[&](int phase,int n,double& fresh) {
            double he,sy;
            auto a=slopes(states,phase,n,he,sy,fresh);
            for(double v:a) if(!std::isfinite(v)) throw std::runtime_error("Nonfinite wall readout");
            half_error=std::max(half_error,he);symmetry=std::max(symmetry,sy);
            return a;
        };
        const auto verify=[&](int n) {
            const auto r=read(states[1],states[0],patch.h,0,n,true);
            energy_error=std::max(energy_error,std::abs(r.direct_force-r.energy_force));
        };
        const auto emit=[](const char* test,int n,double t,const std::array<double,4>& a,double fresh) {
            std::cout<<"{\"test\":\""<<test<<"\",\"arrivals\":"<<n<<",\"time\":"<<t
                     <<",\"positive_wall\":"<<a[0]<<",\"negative_wall\":"<<a[1]
                     <<",\"odd_upstream\":"<<a[2]<<",\"odd_downstream\":"<<a[3]
                     <<",\"fresh\":"<<fresh<<"}\n";
        };
        double fresh;
        auto initial=sample(0,0,fresh);emit("step",0,0.,initial,fresh);verify(0);
        std::array<double,4> sums{};double fresh_sum=0;int samples=0;
        for(int n=0;n<3*patch.length;++n) {
            for(int j=0;j<patch.sub;++j) {
                if(j%(patch.sub/16)==0) {
                    auto a=sample(2*j,n,fresh);
                    for(int k=0;k<4;++k) sums[k]+=a[k];
                    fresh_sum+=fresh;++samples;
                }
                for(size_t k=0;k<states.size();++k) advance(states[k],amplitudes[k],j);
            }
            for(size_t k=0;k<states.size();++k) inject(states[k],amplitudes[k]);
            if((n+1)%2==0) {
                for(double& v:sums) v/=samples;
                emit("two_arrival_mean",n+1,(n+1)*patch.period,sums,fresh_sum/samples);
                sums={};fresh_sum=0;samples=0;
            }
        }
        verify(3*patch.length);
        std::cout<<"{\"test\":\"checks\",\"max_energy_force_difference\":"<<energy_error
                 <<",\"max_full_half_force_difference\":"<<half_error
                 <<",\"max_even_fraction\":"<<symmetry<<"}\n";
        if(energy_error>1e-7 || half_error>1e-7) throw std::runtime_error("Failed energy or weak-drive check");
    }
};

int main(int argc,char** argv) {
    try {
        if(argc!=25) throw std::runtime_error("Use the first-bend Python driver");
        std::cout<<std::setprecision(16);
        Patch patch(22,argv);
        BendHold(patch,argv[22],argv[23],std::stod(argv[24])).run();
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
