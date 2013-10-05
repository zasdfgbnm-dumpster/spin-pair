#include "spin.hpp" 
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace Eigen;
using namespace Tensor;
using namespace spin;
using namespace quantum;

typedef Tensor::vector<double> vec;
typedef Tensor::vector<Operator> veco;

class constants {
public:
	static constexpr double omegab = 2*spin::pi*9.7_GHz;
	
	static constexpr double g = 2.002319304373782;
	static constexpr double g_a = g;
	static constexpr double g_b = g;
	static constexpr double beta = 9.2740089937E-24;
	static constexpr double beta_a = beta;
	static constexpr double beta_b = beta;
	static constexpr double D = 9.4825228E26; /* mu0/(4*pi*hbar) */
};

class spin_pair {
public:
	int index; // start from 0
	
	double r_ab;
	double theta_ab;
	double phi_ab;

	double r_a;
	double theta_a;
	double phi_a;
	
	double r_b() const { return abs(rb()); 	}
	double theta_b() const { return acos(rb()(2)/r_b()); }
	double phi_b() const {
		double c = (rb())(0)/(r_b()*sin(theta_b()));
		double s = (rb())(1)/(r_b()*sin(theta_b()));
		double phi_b = acos(c);
		return phi_b+static_cast<int>(s<0)*pi;
	}
	
	double D_ab() const { return constants::D * constants::g_a * constants::beta_a * constants::g_b * constants::beta_b / pow(r_ab,3); }
	double D_a() const { return constants::D * constants::g_a * constants::beta_a * constants::g * constants::beta / pow(r_a,3); }
	double D_b() const { return constants::D * constants::g_b * constants::beta_b * constants::g * constants::beta / pow(r_b(),3); 	}
	
	vec rab() const { return { r_ab*sin(theta_ab)*cos(phi_ab),r_ab*sin(theta_ab)*sin(phi_ab),r_ab*cos(theta_ab) }; }
	vec ra() const { return { r_a*sin(theta_a)*cos(phi_a),r_a*sin(theta_a)*sin(phi_a),r_a*cos(theta_a) }; }
	vec rb() const { return ra()+rab(); }
	
	/* the value of omega1 when resonance */
	double omegar() const { return -0.5*D_ab()*(1-3*cos(theta_ab)*cos(theta_ab)); }
	/* the angular frequency of oscillation */
	double omegao() const { return 0.5*( D_b()*(1-3*cos(theta_b())*cos(theta_b())) - D_a()*(1-3*cos(theta_a)*cos(theta_a)) ); }
	
	tensor<double> Dab() const {
		tensor<double> I = Tensor::I<double>();
		return D_ab() *(I-3*prod(rab() ,rab() )/pow(r_ab ,2));
	}
	tensor<double> Da () const {
		tensor<double> I = Tensor::I<double>();
		return D_a() *(I-3*prod(ra(),ra())/pow(r_a ,2));
	}
	tensor<double> Db() const {
		tensor<double> I = Tensor::I<double>();
		return D_b() *(I-3*prod(rb() ,rb() )/pow(r_b() ,2));
	}
	
	int subspace() const {
		return index+1;
	}
	int subspace1() const {
		return 2*index+1;
	}
	int subspace2() const {
		return 2*index+2;
	}
	
	double coeff_Sxpes() const {
		return 0.5 * ( Dab()(0,0)+Dab()(1,1) );
	}
	double coeff_Sz0_Szpes() const {
		return Da()(2,2)-Db()(2,2);
	}
	
	Operator H_pair() const {
		veco Sz0 = { O(0),O(0),Sz(0) };
		Operator Hb = constants::omegab*Sz(subspace1())+constants::omegab*Sz(subspace2())+dot(S(subspace1()),Dab(),S(subspace2()));
		Operator Hi = dot(Sz0,Da(),S(subspace1()))+dot(Sz0,Db(),S(subspace2()));
		return Hb + Hi;
	}
	Operator H_pes() const {
		return coeff_Sxpes() * Sx(subspace()) + coeff_Sz0_Szpes() * Sz(0) * Sz(subspace());
	}
	Operator rho0_pair() const {
		return spin::I(subspace1())*spin::I(subspace2())/4;
	}
	Operator rho0_pes() const {
		return spin::I(subspace())/2;
	}
	
};

ostream &operator<<(ostream &out,spin_pair p) {
	out << "index    = " << p.index << endl;
	out << "r_ab     = " << p.r_ab/1_nm << "nm" << endl;
	out << "theta_ab = " << p.theta_ab/1_deg << " degree" << endl;
	out << "phi_ab   = " << p.phi_ab/1_deg << " degree" << endl;
	out << "r_a      = " << p.r_a/1_nm << "nm" << endl;
	out << "theta_a  = " << p.theta_a/1_deg << " degree" << endl;
	out << "phi_a    = " << p.phi_a/1_deg << " degree" << endl;
	out << "Dab      = " << p.D_ab()/1_MHz << "MHz" << endl;
	out << "Da       = " << p.D_a()/1_kHz << "kHz" << endl;
	out << "Db       = " << p.D_b()/1_kHz << "kHz" << endl;
	out << "omegar   = " << p.omegar()/1_MHz << "MHz" << endl;
	out << "omegao   = " << p.omegao()/1_kHz << "kHz";
	return out;
}

const double t_end = 20_ms;
const int n_points = 2000;

spin_pair random_pair() {
	/* random numbers */
	random_device device;
	default_random_engine engine(device());
	uniform_real_distribution<> dist_phi(0,2*pi);
	uniform_real_distribution<> dist_theta(-1,1);
	/* generate */
	spin_pair p;
	p.r_ab = 2.2_nm;
	p.r_a = 20_nm;
	p.phi_a = dist_phi(engine);
	p.phi_ab = dist_phi(engine);
	p.theta_a = acos(dist_theta(engine));
	p.theta_ab = acos(dist_theta(engine));
	return p;
}

std::vector<spin_pair> random_pairs(int n) {
	std::vector<spin_pair> pairs(n);
	generate(pairs.begin(),pairs.end(),random_pair);
	for(int i=0;i<n;i++)
		pairs[i].index = i;
	return pairs;
}

void pairs_simulation(double omega1,const std::vector<spin_pair> &pairs,ostream &out){
	for_each(pairs.begin(),pairs.end(),[&out](spin_pair p){ out<<p<<endl<<endl;} );
	out << endl;
	/* Hamiltonian */
	std::vector<Operator> H_pairs(pairs.size());
	transform(pairs.begin(),pairs.end(),H_pairs.begin(),[](spin_pair p){return p.H_pair();});
	Operator Hs = omega1*Sx(0);
	Operator H = accumulate(H_pairs.begin(),H_pairs.end(),Hs);
	/* initial state */
	std::vector<Operator> rho0_pairs(pairs.size());
	transform(pairs.begin(),pairs.end(),rho0_pairs.begin(),[](spin_pair p){return p.rho0_pair();});
	Operator rho0e = Op<2>(0,0.5,0.5,0.5,0.5);
	Operator rho0 = accumulate(rho0_pairs.begin(),rho0_pairs.end(),rho0e,multiplies<Operator>());
	/* simulation */
	auto U = H.U();
	double step = t_end/n_points;
	for(int i=0;i<=n_points;i++) {
		double t = i*step;
		Operator rhot = U(t)*rho0*U(-t);
		Operator rhote = rhot;
		for(auto &i:pairs){
			rhote = rhote.tr(i.subspace1());
			rhote = rhote.tr(i.subspace2());
		}
		double sx = real(tr(rhote*Sx(0)));
		out << t/1_us << "\t" << sx << endl;
	}
}

void pes_simulation(double omega1,const std::vector<spin_pair> &pairs,ostream &out){
	for_each(pairs.begin(),pairs.end(),[&out](spin_pair p){ out<<p<<endl<<endl;} );
	out << endl;
	/* Hamiltonian */
	std::vector<Operator> H_pes(pairs.size());
	transform(pairs.begin(),pairs.end(),H_pes.begin(),[](spin_pair p){return p.H_pes();});
	Operator Hs = omega1*Sx(0);
	Operator H = accumulate(H_pes.begin(),H_pes.end(),Hs);
	/* initial state */
	std::vector<Operator> rho0_pes(pairs.size());
	transform(pairs.begin(),pairs.end(),rho0_pes.begin(),[](spin_pair p){return p.rho0_pes();});
	Operator rho0e = Op<2>(0,0.5,0.5,0.5,0.5);
	Operator rho0 = accumulate(rho0_pes.begin(),rho0_pes.end(),rho0e,multiplies<Operator>());
	/* simulation */
	auto U = H.U();
	double step = t_end/n_points;
	for(int i=0;i<=n_points;i++) {
		double t = i*step;
		Operator rhot = U(t)*rho0*U(-t);
		Operator rhote = rhot;
		for(auto &i:pairs)
			rhote = rhote.tr(i.subspace());
		double sx = real(tr(rhote*Sx(0)));
		out << t/1_us << "\t" << sx << "\t" << sx*0.5+0.25 << endl;
	}
}
