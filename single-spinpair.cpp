#include "pair.hpp"

int main(){
	spin_pair p;
	p.r_ab = 2.2_nm;
	p.theta_ab = 90_deg;
	p.theta_a = 90_deg;
	p.r_a = 20_nm;
 	p.phi_a = 0;
	p.phi_ab = 0;
	std::vector<spin_pair> pairs = { p };
	// varied omega1
	for(double r=0.99;r<=1.01;r+=0.001){
		double w1 = r*p.omegar();
		stringstream fnstream;
		fnstream << "w1=" << w1/1_MHz << "MHz.txt";
		string fn = fnstream.str();
		ofstream out(fn);
		pairs_simulation(w1,pairs,out);
		out.close();
	}
	// varied phi_ab
	double w1 = p.omegar();
	for(p.phi_ab=0;p.phi_ab<360_deg;p.phi_ab+=60_deg) {
		stringstream fnstream;
		fnstream << "phi_ab=" << p.phi_ab/1_deg << "deg.txt";
		string fn = fnstream.str();
		ofstream out(fn);
		pairs_simulation(w1,pairs,out);
		out.close();
	}
}