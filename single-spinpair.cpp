#include "spin.hpp" 
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;
using namespace Eigen;
using namespace Tensor;
using namespace spin;
using namespace quantum;

typedef Tensor::vector<double> vec;
typedef Tensor::vector<Operator> veco;

class constants {
public:
	static constexpr double omegab = 2*spin::pi*9.7E9;
	
	static constexpr double g = 2.002319304373782;
	static constexpr double g_a = g;
	static constexpr double g_b = g;
	static constexpr double beta = 9.2740089937E-24;
	static constexpr double beta_a = beta;
	static constexpr double beta_b = beta;
	static constexpr double D = 9.4825228E26; /* mu0/(4*pi*hbar) */
};

class parameters {
public:
	double r_ab;
	double theta_ab;
	double phi_ab;

	double r_a;
	double theta_a;
	double phi_a;
	
	double r_b() const { return abs(rb()); 	}
	double theta_b() const { return acos(rb()(2)); }
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
	
	double omega1;
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
	

};

ostream &operator<<(ostream &out,parameters p) {
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
	out << "omegao   = " << p.omegao()/1_kHz << "kHz" << endl;
	out << "omega1   = " << p.omega1/1_MHz << "MHz";
	return out;
}

double omega1;
double omegar;


const double t_end = 20_ms;
const int n_points = 2000;

void simulation(parameters p,ostream &out){
	out << p << endl;
	/* Hamiltonian */
	double omegab = constants::omegab;
	double omega1 = p.omega1;
	tensor<double> Da  = p.Da() ;
	tensor<double> Db  = p.Db() ;
	tensor<double> Dab = p.Dab();
	veco Sz0 = { O(0),O(0),Sz(0) };
	Operator Hs = omega1*Sx(0);
	Operator Hb = omegab*Sz(1)+omegab*Sz(2)+dot(S(1),Dab,S(2));
	Operator Hi = dot(Sz0,Da,S(1))+dot(Sz0,Db,S(2));
	Operator H = Hs + Hb + Hi;
	/* simulation */
	Operator rho0 = Op<2>(0,0.5,0.5,0.5,0.5)*spin::I(1)*spin::I(2)/4;
	auto U = H.U();
	double step = t_end/n_points;
	for(int i=0;i<=n_points;i++) {
		double t = i*step;
		Operator rhot = U(t)*rho0*U(-t);
		Operator rhote = rhot.tr(1,2);
		double sx = real(tr(rhote*Sx(0)));
		out << t/1_us << "\t" << sx << endl;
	}
}

int main(){
	parameters p;
	p.r_ab = 2.2_nm;
	p.theta_ab = 30_deg;
	p.theta_a = 90_deg;
	p.r_a = 20_nm;
 	p.phi_a = 0;
	p.phi_ab = 0;
	// varied omega1
	for(double r=0.6;r<=1.4;r+=0.1){
		p.omega1 = r*p.omegar();
		stringstream fnstream;
		fnstream << "w1=" << p.omega1/1_MHz << "MHz.txt";
		string fn = fnstream.str();
		ofstream out(fn);
		simulation(p,out);
		out.close();
	}
	// varied phi_ab
	p.omega1 = p.omegar();
	for(p.phi_ab=0;p.phi_ab<360_deg;p.phi_ab+=60_deg) {
		stringstream fnstream;
		fnstream << "phi_ab=" << p.phi_ab/1_deg << "deg.txt";
		string fn = fnstream.str();
		ofstream out(fn);
		simulation(p,out);
		out.close();
	}
}