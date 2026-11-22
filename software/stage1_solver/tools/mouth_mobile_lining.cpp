// Two-way momentum exchange between mobile transverse rows and fresh crossing sites.
// Coordinates are in measured lattice pitches and velocities are pitches per carrier time.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Vec = std::vector<double>;
constexpr double spacing = .86602540378443864676;
constexpr double tau = 6.2831853071795864769;

Vec read_table(int size) {
    Vec v(size);
    if (!std::cin.read(reinterpret_cast<char*>(v.data()), size*sizeof(double)))
        throw std::runtime_error("Incomplete material table");
    return v;
}

struct Packet {
    int birth;
    Vec y, v, a;
};

struct Reading {
    double cosine=0, sine=0, speed=0, square=0, count=0;
    void add(const Packet& p, double wall_y, double wall_v) {
        for (size_t k=0; k<p.y.size(); ++k) {
            double angle=tau*(p.y[k]-wall_y);
            double velocity=p.v[k]-wall_v;
            cosine+=std::cos(angle);
            sine+=std::sin(angle);
            speed+=velocity;
            square+=velocity*velocity;
            ++count;
        }
    }
    void print() const {
        if (count==0) throw std::runtime_error("Empty crossing readout");
        std::cout << "{\"tags\":" << count << ",\"phase_amplitude\":" << std::hypot(cosine,sine)/count
                  << ",\"mean_relative_speed\":" << speed/count
                  << ",\"rms_relative_speed\":" << std::sqrt(square/count) << "}";
    }
};

struct RowContact {
    int row, index;
    double fraction;
};

int main(int argc, char** av) {
    if (argc!=13) throw std::invalid_argument("Expected strip parameters and case");
    const int n=std::stoi(av[1]), stride=std::stoi(av[2]), samples=std::stoi(av[3]);
    const int ny=std::stoi(av[4]), nx=std::stoi(av[5]), cutoff=std::stoi(av[6]), neighbors=std::stoi(av[7]);
    const double dt=std::stod(av[8]), scale=std::stod(av[9]);
    const int gauge=std::stoi(av[10]), steps=n*stride;
    const std::string name=av[12];
    const double phase=std::stod(av[11])+(name=="shift" ? .273 : 0.);
    const bool edge=name.rfind("edge",0)==0;
    const double boost=name=="boost" ? .013 : 0.;
    double drift=.002;
    if (name=="reverse" || name=="edge-reverse") drift=-drift;
    if (name=="half" || name=="edge-half") drift*=.5;
    if (name=="rest") drift=0;
    const double weight=1./samples;
    Vec x=read_table(steps+1), order=read_table(steps+1);
    Vec cross=read_table((nx*cutoff+1)*(ny+1)), within=read_table(neighbors*(ny+1));
    std::vector<std::vector<RowContact>> contacts(steps+1);
    for (int t=0; t<=steps; ++t) {
        for (int j=0; j<n; ++j) {
            double dx=std::abs(x[t]-j*spacing)*nx;
            if (dx>=cutoff*nx) continue;
            int k=static_cast<int>(dx);
            contacts[t].push_back({j,k*(ny+1),dx-k});
        }
    }
    Vec row(n), wall_v(n,boost), wall_a(n);
    for (int j=0;j<n;++j) row[j]=phase+.5*(j%2);
    const Vec start=row;
    std::vector<Packet> packets;
    Reading at_gauge, at_exit;
    double injected=0,exited=0,support=0,coasting_time=0;
    int births=0,departures=0;
    auto forces=[&](int t){
        std::fill(wall_a.begin(),wall_a.end(),0.);
        for (int jump=1; jump<=neighbors; ++jump) {
            int base=(jump-1)*(ny+1);
            for (int j=0; j+jump<n; ++j) {
                double dy=row[j+jump]-row[j];
                double u=(dy-std::floor(dy))*ny;
                int k=std::min(ny-1,static_cast<int>(u));
                double f=scale*(within[base+k]+(u-k)*(within[base+k+1]-within[base+k]));
                wall_a[j]+=f;
                wall_a[j+jump]-=f;
            }
        }
        for (auto& p : packets) {
            int age=t-p.birth;
            double amplitude=(name=="smooth" ? 0. : -scale*order[age]);
            std::fill(p.a.begin(),p.a.end(),0.);
            for (const auto& contact : contacts[age]) {
                double reaction=0;
                for (int i=0; i<samples; ++i) {
                    double dy=p.y[i]-row[contact.row];
                    double u=(dy-std::floor(dy))*ny;
                    int k=std::min(ny-1,static_cast<int>(u));
                    int index=contact.index+k;
                    double lo=cross[index]+(u-k)*(cross[index+1]-cross[index]);
                    index+=ny+1;
                    double hi=cross[index]+(u-k)*(cross[index+1]-cross[index]);
                    double f=amplitude*(lo+contact.fraction*(hi-lo));
                    p.a[i]+=f;
                    reaction-=f*weight;
                }
                wall_a[contact.row]+=reaction;
            }
        }
    };
    auto kick=[&](){
        if (edge) {
            support-=.5*dt*wall_a[0];
            wall_a[0]=0;
        }
        for (int j=0; j<n; ++j) wall_v[j]+=.5*dt*wall_a[j];
        for (auto& p:packets)
            for (int k=0;k<samples;++k) p.v[k]+=.5*dt*p.a[k];
    };
    for (int t=0; t<=2*steps; ++t) {
        for (const auto& p : packets) {
            int age=t-p.birth;
            if (p.birth<steps/2 || p.birth>=steps || (age!=gauge && age!=steps)) continue;
            int j=std::clamp(static_cast<int>(std::lround(x[age]/spacing)),0,n-1);
            (age==gauge ? at_gauge : at_exit).add(p,row[j],wall_v[j]);
        }
        if (!packets.empty() && t-packets.front().birth==steps) {
            const auto& p=packets.front();
            exited+=weight*std::accumulate(p.v.begin(),p.v.end(),0.);
            packets.erase(packets.begin());
            ++departures;
        }
        if (t==2*steps) break;
        if (t%stride==0) {
            Packet packet{t,Vec(samples),Vec(samples,drift+boost),Vec(samples)};
            for (int k=0;k<samples;++k) packet.y[k]=phase+boost*t*dt+(k+.5)*weight;
            packets.push_back(std::move(packet));
            injected+=drift+boost;
            ++births;
        }
        forces(t);
        kick();
        for (int j=0;j<n;++j) row[j]+=dt*wall_v[j];
        for (auto& p:packets)
            for (int k=0;k<samples;++k) p.y[k]+=dt*p.v[k];
        forces(t+1);
        kick();
        coasting_time+=dt;
    }
    double inside=0;
    for (const auto& p:packets) inside+=weight*std::accumulate(p.v.begin(),p.v.end(),0.);
    double wall=std::accumulate(wall_v.begin(),wall_v.end(),0.)-n*boost;
    double residual=wall+inside+exited-injected-support;
    double drift_sum=0,c=0,s=0,max_drift=0;
    for (int j=0;j<n;++j) {
        double offset=row[j]-start[j]-boost*coasting_time;
        drift_sum+=offset;
        c+=std::cos(tau*offset);
        s+=std::sin(tau*offset);
        max_drift=std::max(max_drift,std::abs(offset));
    }
    std::cout<<std::setprecision(14)<<"{\"test\":\"strip_passage\",\"case\":\""<<name
             <<"\",\"fed_packets\":"<<births<<",\"departed_packets\":"<<departures
             <<",\"feed_drift\":"<<drift<<",\"wall_momentum_change\":"<<wall
             <<",\"inside_momentum\":"<<inside<<",\"exported_momentum\":"<<exited
             <<",\"injected_momentum\":"<<injected<<",\"edge_holder_impulse\":"<<support
             <<",\"momentum_residual\":"<<residual<<",\"mean_wall_shift\":"<<drift_sum/n
             <<",\"max_wall_shift\":"<<max_drift<<",\"wall_row_coherence\":"<<std::hypot(c,s)/n
             <<",\"depth_gauge\":"; at_gauge.print();
    std::cout<<",\"outlet\":"; at_exit.print(); std::cout<<"}\n";
    if (!std::isfinite(residual) || std::abs(residual)>1e-9 || births!=2*n || departures!=n+1)
        throw std::runtime_error("Open transverse momentum or packet ledger failed");
}
